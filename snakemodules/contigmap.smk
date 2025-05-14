import os
import tempfile


rule index_supercontigs:
    input:
        RESULTS_DIR + "/supercontigs.fa"
    output:
        temp(RESULTS_DIR + "/supercontigs.fa.amb"),
        temp(RESULTS_DIR + "/supercontigs.fa.ann"),
        temp(RESULTS_DIR + "/supercontigs.fa.bwt"),
        temp(RESULTS_DIR + "/supercontigs.fa.pac"),
        temp(RESULTS_DIR + "/supercontigs.fa.sa")
    conda:
        os.path.join(WORKFLOW_PATH,"snakemodules/envs/bwa.yml")
    resources:
        mem_mb = resources["standard_2G"]["mem"],
        runtime = resources["standard_2G"]["time"]
    threads: 
        threads["single"]
    container:
        containers["popins4snake"]
    log:
        err="logs/contigmap/index_supercontigs.err"
    benchmark:
        "benchmarks/contigmap/index_supercontigs.txt"
    shell:
        "bwa index {RESULTS_DIR}/supercontigs.fa 2> {log.err}"


if config["remove_contamination"] == 'yes':

    rule map_supercontigs:
        input:
            expand("{p}/{f}", p=RESULTS_DIR,f=["supercontigs.fa.amb","supercontigs.fa.ann","supercontigs.fa.bwt","supercontigs.fa.pac","supercontigs.fa.sa"]),
            RESULTS_DIR + "/supercontigs.fa",
            se  = WORK_DIR + "/{sample}/single_clean.fastq",
            pe1 = WORK_DIR + "/{sample}/paired_1_clean.fastq",
            pe2 = WORK_DIR + "/{sample}/paired_2_clean.fastq"
        output:
            temp(WORK_DIR + "/{sample}/contig_mapped_unsorted.sam")       
        conda:
            os.path.join(WORKFLOW_PATH,"snakemodules/envs/bwa.yml")
        resources:
            # Dynamically iterating mem
            mem_mb = lambda wildcards, attempt: resources["dynamic_schedule"]["mem_init"]* (2** (attempt - 1)),
            # mem_mb = resources["bwa"]["mem"],
            runtime = resources["bwa"]["time"]
        threads: 
            threads["multi"]["bwa"]
        container:
            containers["popins4snake"]
        log:
            err="logs/contigmap/{sample}_map_supercontigs.err"
        benchmark:
            "benchmarks/contigmap/{sample}_map_supercontigs_clean.txt"
        shell:
            "bwa mem -a -Y -t {threads} {RESULTS_DIR}/supercontigs.fa {input.pe1} {input.pe2} "
            "> {WORK_DIR}/{wildcards.sample}/contig_mapped_unsorted.sam 2>> {log.err};"
            "bwa mem -a -Y -t {threads} {RESULTS_DIR}/supercontigs.fa {input.se} 2>> {log.err} | "
            "awk '$1 !~ /@/' >> {WORK_DIR}/{wildcards.sample}/contig_mapped_unsorted.sam"


elif config["remove_contamination"] == 'no':  
    rule map_supercontigs:
        input:
            expand("{p}/{f}", p=RESULTS_DIR,f=["supercontigs.fa.amb","supercontigs.fa.ann","supercontigs.fa.bwt","supercontigs.fa.pac","supercontigs.fa.sa"]),
            RESULTS_DIR + "/supercontigs.fa",
            se  =WORK_DIR + "/{sample}/single.fastq",
            pe1 =WORK_DIR + "/{sample}/paired.1.fastq",
            pe2 =WORK_DIR + "/{sample}/paired.2.fastq"   
        output:
            temp(WORK_DIR + "/{sample}/contig_mapped_unsorted.sam")       
        conda:
            os.path.join(WORKFLOW_PATH,"snakemodules/envs/bwa.yml")
        resources:
            # Dynamically iterating mem
            mem_mb = lambda wildcards, attempt: resources["dynamic_schedule"]["mem_init"]* (2** (attempt - 1)),
            # mem_mb = resources["bwa"]["mem"],
            runtime = resources["bwa"]["time"]
        threads: 
            threads["multi"]["bwa"]
        container:
            containers["popins4snake"]
        log:
            err="logs/contigmap/{sample}_map_supercontigs.err"
        benchmark:
            "benchmarks/contigmap/{sample}_map_supercontigs.txt"
        shell:
            "bwa mem -a -Y -t {threads} {RESULTS_DIR}/supercontigs.fa {input.pe1} {input.pe2} "
            "> {WORK_DIR}/{wildcards.sample}/contig_mapped_unsorted.sam 2>> {log.err};"
            "bwa mem -a -Y -t {threads} {RESULTS_DIR}/supercontigs.fa {input.se} 2>> {log.err} | "
            "awk '$1 !~ /@/' >> {WORK_DIR}/{wildcards.sample}/contig_mapped_unsorted.sam"
        
################################################################################

rule name_sort_unsorted:
    input:
        WORK_DIR + "/{sample}/contig_mapped_unsorted.sam"  
    output:
        temp(WORK_DIR + "/{sample}/contig_mapped.bam")
    conda:
        os.path.join(WORKFLOW_PATH,"snakemodules/envs/samtools.yml")
    resources:
        mem_per_thread = resources["samtools_multithread"]["mem_per_thread"],
        mem_mb = lambda wildcards, input, threads, attempt: resources["samtools_multithread"]["mem_per_thread"] * (threads + 2),
        runtime = resources["samtools_multithread"]["time"]
    threads: 
        threads["multi"]["samtools"]
    container:
        containers["popins4snake"]
    log:
        err="logs/contigmap/{sample}_name_sort_unsorted.err"
    benchmark:
        "benchmarks/contigmap/{sample}_name_sort_unsorted.txt"
    shell:
        """
        samtools sort -n -@ {threads} -m {resources.mem_per_thread}M -o {output} {input} 2> {log.err}
        """

################################################################################

if config["remove_contamination"] == 'yes':

    rule remap_merge_set_mate:
        input:
            remapped = WORK_DIR + "/{sample}/contig_mapped.bam",
            non_ref = WORK_DIR + "/{sample}/remapped_non_ref.bam"
        output:
            temp(WORK_DIR + "/{sample}/merged.bam")
        resources:
            mem_mb = resources["standard_2G"]["mem"],
            runtime = resources["standard_2G"]["time"]
        threads: 
            threads["single"]
        container:
            containers["popins4snake"]
        log:
            out="logs/contigmap/{sample}_contaminate_removed_remap_merge_set_mate.out",
            err="logs/contigmap/{sample}_contaminate_removed_remap_merge_set_mate.err"
        benchmark:
            "benchmarks/contigmap/{sample}_contaminate_removed_remap_merge_set_mate.txt"
        shell:
            "{POPINS4SNAKE} merge-bams remapped_non_ref.bam contig_mapped.bam"
            "   --prefix {WORK_DIR} "
            "   --sample {wildcards.sample} "
            "   -o merged.bam "
            "   > {log.out} 2> {log.err} "


elif config["remove_contamination"] == 'no':

    rule popins2_merge_set_mate:
        input:
            remapped = WORK_DIR + "/{sample}/contig_mapped.bam",
            non_ref = WORK_DIR + "/{sample}/non_ref.bam"
        output:
            temp(WORK_DIR + "/{sample}/merged.bam")
        resources:
            mem_mb = resources["standard_2G"]["mem"],
            runtime = resources["standard_2G"]["time"]
        threads: 
            threads["single"]
        container:
            containers["popins4snake"]
        log:
            out="logs/contigmap/{sample}_remap_merge_set_mate.out",
            err="logs/contigmap/{sample}_remap_merge_set_mate.err"
        benchmark:
            "benchmarks/contigmap/{sample}_remap_merge_set_mate.txt"
        shell:
            "{POPINS4SNAKE} merge-bams non_ref.bam contig_mapped.bam"
            "   --prefix {WORK_DIR} "
            "   --sample {wildcards.sample} "
            "   -o merged.bam "
            "   > {log.out} 2> {log.err} "   

################################################################################

rule coordinate_sort_unsorted:
    input:
        WORK_DIR + "/{sample}/merged.bam" 
    output:
        temp(WORK_DIR + "/{sample}/non_ref_new.bam")
    conda:
        os.path.join(WORKFLOW_PATH,"snakemodules/envs/samtools.yml")
    resources:
        mem_per_thread = resources["samtools_multithread"]["mem_per_thread"],
        mem_mb = lambda wildcards, input, threads, attempt: resources["samtools_multithread"]["mem_per_thread"] * (threads + attempt),
        runtime = resources["samtools_multithread"]["time"]
    threads: 
        threads["multi"]["samtools"]
    container:
        containers["popins4snake"]
    log:
        err="logs/contigmap/{sample}_coordinate_sort_unsorted.err"
    benchmark:
        "benchmarks/contigmap/{sample}_coordinate_sort_unsorted.txt"
    shell:
        """
        samtools sort -@ {threads} -m {resources.mem_per_thread}M -o {output} -T {resources.tmpdir} {input} 2> {log.err}
        """


rule index_sorted:
    input:
        WORK_DIR + "/{sample}/non_ref_new.bam"
    output:
        temp(WORK_DIR + "/{sample}/non_ref_new.bam.bai")
    conda:
        os.path.join(WORKFLOW_PATH,"snakemodules/envs/samtools.yml")
    resources:
        mem_mb = resources["samtools"]["mem"],
        runtime = resources["samtools"]["time"]
    threads: 
        threads["single"]
    container:
        containers["popins4snake"]
    log:
        err="logs/contigmap/{sample}_index_sorted.err"
    benchmark:
        "benchmarks/contigmap/{sample}_index_sorted.txt"
    shell:
        "samtools index {input} 2> {log.err}"
