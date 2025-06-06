import os
def find_bam_file(wildcards):
    # Construct the direct path first:
    direct_path = os.path.join(INPUT_DIR, f"{wildcards.sample}.bam")
    if os.path.exists(direct_path):
        return direct_path
    else:
        # If not found directly in INPUT_DIR, search recursively in subdirectories
        search_pattern = os.path.join(INPUT_DIR, '**', f"{wildcards.sample}.bam")
        bam_matches = glob.glob(search_pattern, recursive=True)
        if bam_matches:
            # If multiple matches, choose the first or define another selection criterion
            return bam_matches[0]
        else:
            # Raise an error if no bam file found for this sample
            raise FileNotFoundError(f"No BAM file found for sample {wildcards.sample}")

if config["SICKLE"]=="no":
    rule popins2_crop_unmapped:
        input:
            find_bam_file
            # INPUT_DIR + "/{sample}/{sample}.bam"       
        output:
            temp(WORK_DIR + "/{sample}/single.fastq"), 
            temp(WORK_DIR + "/{sample}/paired.1.fastq"), 
            temp(WORK_DIR + "/{sample}/paired.2.fastq"), 
            temp(WORK_DIR + "/{sample}/mates.bam"), 
            WORK_DIR + "/{sample}/POPINS_SAMPLE_INFO"
        params:
            q = config["min-qual"],
            l = config["min-read-len"]
        resources:
            mem_mb = lambda wildcards, attempt: resources["dynamic_schedule"]["mem_init"]* (2** (attempt - 1)),
            # mem_mb = resources["standard_4G"]["mem"],
            runtime = resources["standard_4G"]["time"]
        threads: 
            threads["single"]
        log:
            out="logs/crop_unmapped/{sample}_crop_unmapped.out",
            err="logs/crop_unmapped/{sample}_crop_unmapped.err"
        benchmark:
            "benchmarks/crop_unmapped/{sample}_crop_unmapped.txt"
        shell:
            "mkdir -p {OUTPUT_PATH}/{WORK_DIR}/{wildcards.sample}; "
            "{POPINS2_BIN} crop-unmapped {input}"
            "   --prefix {WORK_DIR} "
            "   --sample {wildcards.sample}"
            "   --min-qual {params.q}"
            "   --min-read-len {params.l}"
            "   > {log.out} 2> {log.err}"
            
    rule popins2_sort:
        input:
            os.path.join(RESULTS_DIR, "read_numbers.pdf") if config["ANALYSIS"]=="yes" else [],
            mates = WORK_DIR + "/{sample}/mates.bam"
        output:
            temp(WORK_DIR + "/{sample}/non_ref.bam")
        conda:
            os.path.join(WORKFLOW_PATH,"snakemodules/envs/samtools.yml")
        resources:
            mem_per_thread = resources["samtools_multithread"]["mem_per_thread"],
            mem_mb = lambda wildcards, input, threads, attempt: resources["samtools_multithread"]["mem_per_thread"] * threads,
            runtime = resources["samtools_multithread"]["time"]
        threads: 
            threads["multi"]["samtools"]
        log:
            err="logs/crop_unmapped/{sample}_sort.err"
        benchmark:
            "benchmarks/crop_unmapped/{sample}_sort.txt"
        shell:
            "samtools sort -n -@ {threads} -m {resources.mem_per_thread}M -o {output} -T {resources.tmpdir} {input.mates} 2> {log.err}" 


elif config["SICKLE"]=="yes":
    rule popins2_crop_unmapped:
        input:
            INPUT_DIR + "/{sample}/{sample}.bam"       
        output:
            temp(WORK_DIR + "/{sample}/single.fastq"), 
            temp(WORK_DIR + "/{sample}/paired.1.fastq"), 
            temp(WORK_DIR + "/{sample}/paired.2.fastq"), 
            temp(WORK_DIR + "/{sample}/mates.bam"), 
            WORK_DIR + "/{sample}/POPINS_SAMPLE_INFO"
        params:
            q = config["min-qual"],
            l = config["min-read-len"]
        resources:
            mem_mb = lambda wildcards, attempt: resources["dynamic_schedule"]["mem_init"]* (2** (attempt - 1)),
            # mem_mb = resources["standard_4G"]["mem"],
            runtime = resources["standard_4G"]["time"]
        threads: 
            threads["single"]
        log:
            out="logs/crop_unmapped/{sample}_crop_unmapped_sickle.out",
            err="logs/crop_unmapped/{sample}_crop_unmapped_sickle.err"
        benchmark:
            "benchmarks/crop_unmapped/{sample}_crop_unmapped_sickle.txt"
        shell:
            "mkdir -p {WORK_DIR}/{wildcards.sample}; "
            "{POPINS2_BIN} crop-unmapped {input}"
            "   --prefix {WORK_DIR} "
            "   --sample {wildcards.sample}"
            "   > {log.out} 2> {log.err}"
            
    rule popins2_sort:
        input:
            # os.path.join(RESULTS_DIR, "read_numbers.pdf") if config["ANALYSIS"]=="yes" else [],
            mates = WORK_DIR + "/{sample}/mates.bam"
        output:
            temp(WORK_DIR + "/{sample}/non_ref.bam")
        resources:
            mem_per_thread = resources["samtools_multithread"]["mem_per_thread"],
            mem_mb = lambda wildcards, input, threads, attempt: resources["samtools_multithread"]["mem_per_thread"] * threads,
            runtime = resources["samtools_multithread"]["time"]
        threads: 
            threads["multi"]["samtools"]
        conda:
            os.path.join(WORKFLOW_PATH,"snakemodules/envs/samtools.yml")
        log:
            err="logs/crop_unmapped/{sample}_sort_sickle.err"
        benchmark:
            "benchmarks/crop_unmapped/{sample}_sort_sickle.txt"
        shell:
            "samtools sort -n -@ {threads} -m {resources.mem_per_thread}M -o {output} -T {resources.tmpdir} {input.mates} 2> {log.err}"     
     
    rule popins2_sickle:
        input:
            # os.path.join(RESULTS_DIR, "read_numbers.pdf") if config["ANALYSIS"]=="yes" else [],
            single = WORK_DIR + "/{sample}/single.fastq", 
            pair_1 = WORK_DIR + "/{sample}/paired.1.fastq", 
            pair_2 = WORK_DIR + "/{sample}/paired.2.fastq"
        output:
            single = temp(WORK_DIR + "/{sample}/sickle.single.fastq"), 
            pair_1 = temp(WORK_DIR + "/{sample}/sickle.paired.1.fastq"), 
            pair_2 = temp(WORK_DIR + "/{sample}/sickle.paired.2.fastq"),
            temp_single = temp(WORK_DIR + "/{sample}/sickle.single2.fastq")
        params:
            q = config["min-qual"],
            l = config["min-read-len"]
        resources:
            mem_mb = resources["standard_2G"]["mem"],
            runtime = resources["standard_2G"]["time"]
        threads: 
            threads["single"]
        log:
            out="logs/crop_unmapped/{sample}_sickle.out",
            err="logs/crop_unmapped/{sample}_sickle.err"
        benchmark:
            "benchmarks/crop_unmapped/{sample}_sickle.txt"
        shell:
            "{SICKLE}"
            "   pe -q {params.q} -l {params.l} -x -n -t sanger"
            "   -f {input.pair_1} -r {input.pair_2}"
            "   -o {output.pair_1} -p {output.pair_2} -s {output.single}"
            "   > {log.out} 2> {log.err}; "
            "{SICKLE}"
            "   se -q {params.q} -l {params.l} -x -n -t sanger"
            "   -f {input.single} -o {output.temp_single}"
            "   >> {log.out} 2>> {log.err}; "
            "cat {output.temp_single} >> {output.single}"
            
