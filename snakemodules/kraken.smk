import os

rule kraken_map:
    input:
        # directory(f"{kraken_db_path}"),
        fq = WORK_DIR + "/{sample}/single.fastq",
        p1 = WORK_DIR + "/{sample}/paired.1.fastq", 
        p2 = WORK_DIR + "/{sample}/paired.2.fastq"
    output:
        single_report = WORK_DIR + "/{sample}/contaminate_info/single_report.txt",
        single_output = temp(WORK_DIR + "/{sample}/contaminate_info/single_output.txt"),
        u_single = temp(WORK_DIR + "/{sample}/contaminate_info/useq_single.fq"),
        c_single = temp(WORK_DIR + "/{sample}/contaminate_info/cseq_single.fq"),
        paired_report = WORK_DIR + "/{sample}/contaminate_info/paired_report.txt",
        paired_output = temp(WORK_DIR + "/{sample}/contaminate_info/paired_output.txt"),
        c_paired_1 = temp(WORK_DIR + "/{sample}/contaminate_info/cseqs_1.fq"),
        c_paired_2 = temp(WORK_DIR + "/{sample}/contaminate_info/cseqs_2.fq"),
        u_paired_1 = temp(WORK_DIR + "/{sample}/contaminate_info/useqs_1.fq"),
        u_paired_2 = temp(WORK_DIR + "/{sample}/contaminate_info/useqs_2.fq")
    conda:
        os.path.join(WORKFLOW_PATH,"snakemodules/envs/kraken2.yml")
    resources:
        mem_mb = lambda wildcards, attempt: resources["kraken"]["mem"] * attempt,
        runtime = resources["kraken"]["time"]
    threads: 
        threads["multi"]["kraken"]
    container:
        containers["popins4snake"]
    log:
        out="logs/kraken/{sample}_kraken_map.out",
        err="logs/kraken/{sample}_kraken_map.err"
    benchmark:
        "benchmarks/kraken/{sample}_kraken_map.txt"
    shell:
        """
        mkdir -p {WORK_DIR}/{wildcards.sample}/contaminate_info;
        kraken2 \
        --db {kraken_db_path} \
        --threads {threads} \
        --report {output.single_report} \
        --output {output.single_output} \
        --classified-out {output.c_single} \
        --unclassified-out {output.u_single} \
        {input.fq} \
        > {log.out} 2> {log.err};
        kraken2 \
        --db {kraken_db_path} \
        --threads {threads} \
        --report {output.paired_report} \
        --output {output.paired_output} \
        --classified-out {WORK_DIR}/{wildcards.sample}/contaminate_info/cseqs#.fq \
        --unclassified-out {WORK_DIR}/{wildcards.sample}/contaminate_info/useqs#.fq \
        --paired {input.p1} {input.p2} \
        >> {log.out} 2>> {log.err}
        """

rule obtain_classified_human_reads:
    input:
        classified_single = rules.kraken_map.output.single_output,
        classified_paired = rules.kraken_map.output.paired_output,
        single_fq = rules.kraken_map.output.c_single,
        paired_fq1 = rules.kraken_map.output.c_paired_1,
        paired_fq2 = rules.kraken_map.output.c_paired_2,
        script= os.path.join(WORKFLOW_PATH,config["python_script"])
    output:
        hc_single = temp(WORK_DIR + "/{sample}/contaminate_info/human_classified_single.fastq"),
        hc_paired_1 = temp(WORK_DIR + "/{sample}/contaminate_info/human_classified_paired_1.fastq"),
        hc_paired_2 = temp(WORK_DIR + "/{sample}/contaminate_info/human_classified_paired_2.fastq")
    conda:
        os.path.join(WORKFLOW_PATH,"snakemodules/envs/biopython.yml")
    resources:
        mem_mb = resources["standard_2G"]["mem"],
        runtime = resources["standard_2G"]["time"]
    threads: 
        threads["single"]
    container:
        containers["popins4snake"]
    log:
        out="logs/kraken/{sample}_obtain_classified_human_reads.out",
        err="logs/kraken/{sample}_obtain_classified_human_reads.err"
    benchmark:
        "benchmarks/kraken/{sample}_obtain_classified_human_reads.txt"
    shell:
        """
        python3 {input.script} \
            --classified_single {input.classified_single} \
            --classified_paired {input.classified_paired} \
            --single_fq {input.single_fq} \
            --paired_fq1 {input.paired_fq1} \
            --paired_fq2 {input.paired_fq2} \
            --output_prefix {WORK_DIR}/{wildcards.sample}/contaminate_info \
            > {log.out} 2> {log.err}
        """

rule index_ref:
    input:
        ref = config["ALTREF"]
    output:
        f"{ALTREF}.amb",
        f"{ALTREF}.ann",
        f"{ALTREF}.bwt",
        f"{ALTREF}.pac",
        f"{ALTREF}.sa"
    conda:
        os.path.join(WORKFLOW_PATH,"snakemodules/envs/bwa.yml")
    resources:
        mem_mb = resources["standard_8G"]["mem"],
        runtime = resources["standard_8G"]["time"]
    threads: 
        threads["single"]
    container:
        containers["popins4snake"]
    log:
        err="logs/kraken/index_ref.err"
    benchmark:
        "benchmarks/kraken/index_ref.txt"
    shell:
        "bwa index {input.ref} 2> {log.err}"

rule bwa_remap_classified_human:
    input:
        index = rules.index_ref.output,
        ref = config["ALTREF"],
        single = rules.obtain_classified_human_reads.output.hc_single,
        fq1 = rules.obtain_classified_human_reads.output.hc_paired_1,
        fq2 = rules.obtain_classified_human_reads.output.hc_paired_1
    output:
        remapped_sam = temp(WORK_DIR + "/{sample}/contaminate_info/remapped_human.sam"),
    conda:
        os.path.join(WORKFLOW_PATH,"snakemodules/envs/bwa.yml")
    resources:
        mem_mb = resources["bwa"]["mem"],
        runtime = resources["bwa"]["time"]
    threads: 
        threads["multi"]["bwa"]
    container:
        containers["popins4snake"]
    log:
        err="logs/kraken/{sample}_bwa_remap_classified_human.err"
    benchmark:
        "benchmarks/kraken/{sample}_bwa_remap_classified_human.txt"
    shell:
        """
        bwa mem -t {threads} -M -R "@RG\\tID:{wildcards.sample}\\tPL:ILLUMINA\\tSM:{wildcards.sample}" \
        {input.ref} {input.single} > {output.remapped_sam} 2>> {log.err};
        bwa mem -t {threads} -M {input.ref} {input.fq1} {input.fq2} 2>> {log.err} | \
        awk '$1 !~ /@/' >> {output.remapped_sam}
        """

rule samtools_remap_classified_human:
    input:
        sam = rules.bwa_remap_classified_human.output.remapped_sam
    output:
        remapped_unsorted = temp(WORK_DIR + "/{sample}/contaminate_info/remapped_human_unsorted.bam"),
        remapped_bam = temp(WORK_DIR + "/{sample}/contaminate_info/remapped.bam")
    conda:
        os.path.join(WORKFLOW_PATH,"snakemodules/envs/samtools.yml")
    resources:
        mem_per_thread = resources["samtools_multithread"]["mem_per_thread"],
        mem_mb = lambda wildcards, input, threads, attempt: resources["samtools_multithread"]["mem_per_thread"] * threads,
        runtime = resources["samtools_multithread"]["time"]
    threads: 
        threads["multi"]["samtools"]
    container:
        containers["popins4snake"]
    log:
        err="logs/kraken/{sample}_samtools_remap_classified_human.err"
    benchmark:
        "benchmarks/kraken/{sample}_samtools_remap_classified_human.txt"
    shell:
        """
        samtools view -Sb {input.sam} > {output.remapped_unsorted} 2> {log.err};
        samtools sort -@ {threads} -m {resources.mem_per_thread}M -o {output.remapped_bam} {output.remapped_unsorted} 2>> {log.err}
        """
       
rule index_reads:
    input:
        bam = rules.samtools_remap_classified_human.output.remapped_bam
    output:
        bai = temp(WORK_DIR + "/{sample}/contaminate_info/remapped.bam.bai")
    conda:
        os.path.join(WORKFLOW_PATH,"snakemodules/envs/samtools.yml")
    resources:
        mem_mb = resources["standard_2G"]["mem"],
        runtime = resources["standard_2G"]["time"]
    threads: 
        threads["single"]
    container:
        containers["popins4snake"]
    log:
        err="logs/kraken/{sample}_index_reads.err"
    benchmark:
        "benchmarks/kraken/{sample}_index_reads.txt"
    shell:
        """
        samtools index {input.bam} 2> {log.err}
        """

rule crop_remapped:
    input:
        bam = rules.samtools_remap_classified_human.output.remapped_bam,
        bai = rules.index_reads.output.bai
    output:
        single = temp(WORK_DIR + "/{sample}/remapped_single.fastq"), 
        p1 = temp(WORK_DIR + "/{sample}/remapped_paired.1.fastq"), 
        p2 = temp(WORK_DIR + "/{sample}/remapped_paired.2.fastq"), 
        unsorted = temp(WORK_DIR + "/{sample}/remapped_mates.unsorted.bam")
    resources:
        mem_mb = resources["standard_2G"]["mem"],
        runtime = resources["standard_2G"]["time"]
    threads: 
        threads["single"]
    container:
        containers["popins4snake"]
    log:
        out="logs/kraken/{sample}_crop_remapped.out",
        err="logs/kraken/{sample}_crop_remapped.err"
    benchmark:
        "benchmarks/kraken/{sample}_crop_remapped.txt"    
    shell:
        "{POPINS4SNAKE} crop-unmapped {input.bam}"
        "   -m remapped_mates.unsorted.bam"
        "   -pe1 remapped_paired.1.fastq"
        "   -pe2 remapped_paired.2.fastq"
        "   -se remapped_single.fastq"
        "   --prefix {WORK_DIR}"
        "   --noSampleInfo"
        "   > {log.out} 2> {log.err}"


rule remapping_samsort_mates:
    input:
        WORK_DIR + "/{sample}/remapped_mates.unsorted.bam"
    output:
        temp(WORK_DIR + "/{sample}/remapped_mates.bam")
    conda:
        os.path.join(WORKFLOW_PATH,"snakemodules/envs/samtools.yml")
    resources:
        mem_per_thread = resources["samtools_multithread"]["mem_per_thread"],
        mem_mb = lambda wildcards, input, threads, attempt: resources["samtools_multithread"]["mem_per_thread"] * threads,
        runtime = resources["samtools_multithread"]["time"]
    threads: 
        threads["multi"]["samtools"]
    container:
        containers["popins4snake"]
    log:
        err="logs/kraken/{sample}_remapping_samsort_mates.err"
    benchmark:
        "benchmarks/kraken/{sample}_remapping_samsort_mates.txt"    
    shell:
        "samtools sort -n -@ {threads} -m {resources.mem_per_thread}M -o {output} {input} 2> {log.err}"


rule merge_set_mate:
    input:
        non_ref = WORK_DIR + "/{sample}/non_ref.bam",
        rm_ref = WORK_DIR + "/{sample}/remapped_mates.bam"
    output:
        temp(WORK_DIR + "/{sample}/remapped_non_ref.bam")
    resources:
        mem_mb = resources["standard_2G"]["mem"],
        runtime = resources["standard_2G"]["time"]
    threads: 
        threads["single"]
    container:
        containers["popins4snake"]
    log:
        out="logs/kraken/{sample}_merge_set_mate.out",
        err="logs/kraken/{sample}_merge_set_mate.err"
    benchmark:
        "benchmarks/kraken/{sample}_merge_set_mate.txt" 
    shell:
        "{POPINS4SNAKE} merge-bams non_ref.bam remapped_mates.bam"
        "   --prefix {WORK_DIR} "
        "   --sample {wildcards.sample} "
        "   -o remapped_non_ref.bam"
        "   > {log.out} 2> {log.err}"


rule contaminate_removed:
    input:
        single = rules.crop_remapped.output.single,
        p1 = rules.crop_remapped.output.p1,
        p2 = rules.crop_remapped.output.p2,
        u_single = rules.kraken_map.output.u_single,
        u_p1 = rules.kraken_map.output.u_paired_1,
        u_p2 = rules.kraken_map.output.u_paired_2
    output:
        single_clean = temp(WORK_DIR + "/{sample}/single_clean.fastq"),
        p1_clean = temp(WORK_DIR + "/{sample}/paired_1_clean.fastq"),
        p2_clean = temp(WORK_DIR + "/{sample}/paired_2_clean.fastq")
    resources:
        mem_mb = resources["standard_2G"]["mem"],
        runtime = resources["standard_2G"]["time"]
    threads: 
        threads["single"]
    log:
        "logs/kraken/{sample}_contaminate_removed.log"
    benchmark:
        "benchmarks/kraken/{sample}_contaminate_removed.txt" 
    shell:
        """
        cat {input.single} {input.u_single} > {output.single_clean};
        cat {input.p1} {input.u_p1} > {output.p1_clean};
        cat {input.p2} {input.u_p2} > {output.p2_clean};
        echo "Contaminates removed for {wildcards.sample}" >> {log}
        """

