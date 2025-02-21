import os

if config["SICKLE"]=="yes":
    if config["remove_contamination"] == 'yes':
        rule popins2_minia:
            input:
                single = WORK_DIR + "/{sample}/sickle.single_clean.fastq", 
                pair_1 = WORK_DIR + "/{sample}/sickle.paired_1_clean.fastq", 
                pair_2 = WORK_DIR + "/{sample}/sickle.paired_2_clean.fastq"
            output:
                temp(WORK_DIR + "/{sample}/{assemblr}.contigs.fa")
            conda:
                os.path.join(WORKFLOW_PATH,"snakemodules/envs/py27.yml")
            resources:
                mem_mb = resources["standard_8G"]["mem"],
                runtime = resources["standard_8G"]["time"]
            threads: 
                threads["multi"]["minia"]
            log:
                out="logs/{assemblr}/{sample}_sickle_clean.out",
                err="logs/{assemblr}/{sample}_sickle_clean.err"
            benchmark:
                "benchmarks/{assemblr}/{sample}_sickle_clean.txt"
            shell:
                "mkdir -p {WORK_DIR}/{wildcards.sample}/assembly/;"  
                "{GATB} --nb-cores {threads} -1 {input.pair_1} -2 {input.pair_2} -s {input.single}"
                " -o {WORK_DIR}/{wildcards.sample}/assembly/assembly --no-scaffolding > {log.out} 2> {log.err};"
                "cp {WORK_DIR}/{wildcards.sample}/assembly/assembly_final.contigs.fa {WORK_DIR}/{wildcards.sample}/{wildcards.assemblr}.contigs.fa;"     
                "rm -r {WORK_DIR}/{wildcards.sample}/assembly/"

    elif config["remove_contamination"] == 'no':

        rule popins2_minia:
            input:
                single = WORK_DIR + "/{sample}/sickle.single.fastq", 
                pair_1 = WORK_DIR + "/{sample}/sickle.paired.1.fastq", 
                pair_2 = WORK_DIR + "/{sample}/sickle.paired.2.fastq"
            output:
                temp(WORK_DIR + "/{sample}/{assemblr}.contigs.fa")
            conda:
                os.path.join(WORKFLOW_PATH,"snakemodules/envs/py27.yml")
            resources:
                mem_mb = resources["standard_8G"]["mem"],
                runtime = resources["standard_8G"]["time"]
            threads: 
                threads["multi"]["minia"]
            log:
                out="logs/{assemblr}/{sample}_sickle.out",
                err="logs/{assemblr}/{sample}_sickle.err"
            benchmark:
                "benchmarks/{assemblr}/{sample}_sickle.txt"
            shell:
                "mkdir -p {WORK_DIR}/{wildcards.sample}/assembly/;"  
                "{GATB} --nb-cores {threads} -1 {input.pair_1} -2 {input.pair_2} -s {input.single}"
                " -o {WORK_DIR}/{wildcards.sample}/assembly/assembly --no-scaffolding > {log.out} 2> {log.err};"
                "cp {WORK_DIR}/{wildcards.sample}/assembly/assembly_final.contigs.fa {WORK_DIR}/{wildcards.sample}/{wildcards.assemblr}.contigs.fa;"     
                "rm -r {WORK_DIR}/{wildcards.sample}/assembly/"

elif config["SICKLE"]=="no":
    if config["remove_contamination"] == 'yes':
        rule popins2_minia:
            input:
                single = WORK_DIR + "/{sample}/single_clean.fastq", 
                pair_1 = WORK_DIR + "/{sample}/paired_1_clean.fastq", 
                pair_2 = WORK_DIR + "/{sample}/paired_2_clean.fastq"
            output:
                temp(WORK_DIR + "/{sample}/{assemblr}.contigs.fa")
            conda:
                os.path.join(WORKFLOW_PATH,"snakemodules/envs/py27.yml")
            resources:
                mem_mb = resources["standard_8G"]["mem"],
                runtime = resources["standard_8G"]["time"]
            threads: 
                threads["multi"]["minia"]
            log:
                out="logs/{assemblr}/{sample}_clean.out",
                err="logs/{assemblr}/{sample}_clean.err"
            benchmark:
                "benchmarks/{assemblr}/{sample}_clean.txt"
            shell:
                "mkdir -p {WORK_DIR}/{wildcards.sample}/assembly/;"  
                "{GATB} --nb-cores {threads} -1 {input.pair_1} -2 {input.pair_2} -s {input.single}"
                " -o {WORK_DIR}/{wildcards.sample}/assembly/assembly --no-scaffolding > {log.out} 2> {log.err};"
                "cp {WORK_DIR}/{wildcards.sample}/assembly/assembly_final.contigs.fa {WORK_DIR}/{wildcards.sample}/{wildcards.assemblr}.contigs.fa;"     
                "rm -r {WORK_DIR}/{wildcards.sample}/assembly/"

    elif config["remove_contamination"] == 'no':

        rule popins2_minia:
            input:
                single = WORK_DIR + "/{sample}/single.fastq", 
                pair_1 = WORK_DIR + "/{sample}/paired.1.fastq", 
                pair_2 = WORK_DIR + "/{sample}/paired.2.fastq"
            output:
                temp(WORK_DIR + "/{sample}/{assemblr}.contigs.fa")
            conda:
                os.path.join(WORKFLOW_PATH,"snakemodules/envs/py27.yml")
            resources:
                mem_mb = resources["standard_8G"]["mem"],
                runtime = resources["standard_8G"]["time"]
            threads: 
                threads["multi"]["minia"]
            log:
                out="logs/{assemblr}/{sample}.out",
                err="logs/{assemblr}/{sample}.err"
            benchmark:
                "benchmarks/{assemblr}/{sample}.txt"
            shell:
                "mkdir -p {WORK_DIR}/{wildcards.sample}/assembly/;"  
                "{GATB} --nb-cores {threads} -1 {input.pair_1} -2 {input.pair_2} -s {input.single}"
                " -o {WORK_DIR}/{wildcards.sample}/assembly/assembly --no-scaffolding > {log.out} 2> {log.err};"
                "cp {WORK_DIR}/{wildcards.sample}/assembly/assembly_final.contigs.fa {WORK_DIR}/{wildcards.sample}/{wildcards.assemblr}.contigs.fa;"     
                "rm -r {WORK_DIR}/{wildcards.sample}/assembly/"
