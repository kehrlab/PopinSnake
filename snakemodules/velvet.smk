import os
if config["SICKLE"]=="yes":
    if config["remove_contamination"] == 'yes':
        rule popins2_velvet:
            input:
                single = WORK_DIR + "/{sample}/sickle.single_clean.fastq", 
                pair_1 = WORK_DIR + "/{sample}/sickle.paired_1_clean.fastq", 
                pair_2 = WORK_DIR + "/{sample}/sickle.paired_2_clean.fastq"
            output:
                temp(WORK_DIR + "/{sample}/{assemblr}.contigs.fa")
            params:
                kmerlength = config["kmerlength"]
            conda:
                os.path.join(WORKFLOW_PATH,"snakemodules/envs/velvet.yml")
            resources:
                mem_mb = resources["standard_4G"]["mem"],
                runtime = resources["standard_4G"]["time"]
            threads: 
                threads["single"]
            log:
                out="logs/{assemblr}/{sample}_sickle_clean.out",
                err="logs/{assemblr}/{sample}_sickle_clean.err"
            benchmark:
                "benchmarks/{assemblr}/{sample}_sickle_clean.txt"
            shell:
                "mkdir {WORK_DIR}/{wildcards.sample}/assembly/;"
                "{VELVET}/velveth {WORK_DIR}/{wildcards.sample}/assembly/ {params.kmerlength} -short -fastq {input.single} -shortPaired -fastq -separate {input.pair_1} {input.pair_2} > {log.out} 2> {log.err};"
                "{VELVET}/velvetg {WORK_DIR}/{wildcards.sample}/assembly/ -exp_cov auto -cov_cutoff 2 -max_coverage 100 -scaffolding no >> {log.out} 2>> {log.err};"    
                "mv {WORK_DIR}/{wildcards.sample}/assembly/contigs.fa {WORK_DIR}/{wildcards.sample}/{wildcards.assemblr}.contigs.fa;"
                "rm -r {WORK_DIR}/{wildcards.sample}/assembly/"

    elif config["remove_contamination"] == 'no':

        rule popins2_velvet:
            input:
                single = WORK_DIR + "/{sample}/sickle.single.fastq", 
                pair_1 = WORK_DIR + "/{sample}/sickle.paired.1.fastq", 
                pair_2 = WORK_DIR + "/{sample}/sickle.paired.2.fastq"
            output:
                temp(WORK_DIR + "/{sample}/{assemblr}.contigs.fa")
            params:
                kmerlength = config["kmerlength"]
            conda:
                os.path.join(WORKFLOW_PATH,"snakemodules/envs/velvet.yml")
            resources:
                mem_mb = resources["standard_4G"]["mem"],
                runtime = resources["standard_4G"]["time"]
            threads: 
                threads["single"]
            log:
                out="logs/{assemblr}/{sample}_sickle.out",
                err="logs/{assemblr}/{sample}_sickle.err"
            benchmark:
                "benchmarks/{assemblr}/{sample}_sickle.txt"
            shell:
                "mkdir {WORK_DIR}/{wildcards.sample}/assembly/;"
                "{VELVET}/velveth {WORK_DIR}/{wildcards.sample}/assembly/ {params.kmerlength} -short -fastq {input.single} -shortPaired -fastq -separate {input.pair_1} {input.pair_2} > {log.out} 2> {log.err};"
                "{VELVET}/velvetg {WORK_DIR}/{wildcards.sample}/assembly/ -exp_cov auto -cov_cutoff 2 -max_coverage 100 -scaffolding no >> {log.out} 2>> {log.err};"    
                "mv {WORK_DIR}/{wildcards.sample}/assembly/contigs.fa {WORK_DIR}/{wildcards.sample}/{wildcards.assemblr}.contigs.fa;"
                "rm -r {WORK_DIR}/{wildcards.sample}/assembly/"

elif config["SICKLE"]=="no":
    if config["remove_contamination"] == 'yes':
        rule popins2_velvet:
            input:
                single = WORK_DIR + "/{sample}/single_clean.fastq", 
                pair_1 = WORK_DIR + "/{sample}/paired_1_clean.fastq", 
                pair_2 = WORK_DIR + "/{sample}/paired_2_clean.fastq"
            output:
                temp(WORK_DIR + "/{sample}/{assemblr}.contigs.fa")
            params:
                kmerlength = config["kmerlength"]
            conda:
                os.path.join(WORKFLOW_PATH,"snakemodules/envs/velvet.yml")
            resources:
                mem_mb = resources["standard_4G"]["mem"],
                runtime = resources["standard_4G"]["time"]
            threads: 
                threads["single"]
            log:
                out="logs/{assemblr}/{sample}_clean.out",
                err="logs/{assemblr}/{sample}_clean.err"
            benchmark:
                "benchmarks/{assemblr}/{sample}_clean.txt"
            shell:
                "mkdir {WORK_DIR}/{wildcards.sample}/assembly/;"
                "{VELVET}/velveth {WORK_DIR}/{wildcards.sample}/assembly/ {params.kmerlength} -short -fastq {input.single} -shortPaired -fastq -separate {input.pair_1} {input.pair_2} > {log.out} 2> {log.err};"
                "{VELVET}/velvetg {WORK_DIR}/{wildcards.sample}/assembly/ -exp_cov auto -cov_cutoff 2 -max_coverage 100 -scaffolding no >> {log.out} 2>> {log.err};"    
                "mv {WORK_DIR}/{wildcards.sample}/assembly/contigs.fa {WORK_DIR}/{wildcards.sample}/{wildcards.assemblr}.contigs.fa;"
                "rm -r {WORK_DIR}/{wildcards.sample}/assembly/"

    elif config["remove_contamination"] == 'no':

        rule popins2_velvet:
            input:
                single = WORK_DIR + "/{sample}/single.fastq", 
                pair_1 = WORK_DIR + "/{sample}/paired.1.fastq", 
                pair_2 = WORK_DIR + "/{sample}/paired.2.fastq"
            output:
                temp(WORK_DIR + "/{sample}/{assemblr}.contigs.fa")
            params:
                kmerlength = config["kmerlength"]
            conda:
                os.path.join(WORKFLOW_PATH,"snakemodules/envs/velvet.yml")
            resources:
                mem_mb = resources["standard_4G"]["mem"],
                runtime = resources["standard_4G"]["time"]
            threads: 
                threads["single"]
            log:
                out="logs/{assemblr}/{sample}.out",
                err="logs/{assemblr}/{sample}.err"
            benchmark:
                "benchmarks/{assemblr}/{sample}.txt"
            shell:
                "mkdir {WORK_DIR}/{wildcards.sample}/assembly/;"
                "{VELVET}/velveth {WORK_DIR}/{wildcards.sample}/assembly/ {params.kmerlength} -short -fastq {input.single} -shortPaired -fastq -separate {input.pair_1} {input.pair_2} > {log.out} 2> {log.err};"
                "{VELVET}/velvetg {WORK_DIR}/{wildcards.sample}/assembly/ -exp_cov auto -cov_cutoff 2 -max_coverage 100 -scaffolding no >> {log.out} 2>> {log.err};"    
                "mv {WORK_DIR}/{wildcards.sample}/assembly/contigs.fa {WORK_DIR}/{wildcards.sample}/{wildcards.assemblr}.contigs.fa;"
                "rm -r {WORK_DIR}/{wildcards.sample}/assembly/"