import os
rule sort_vcf:
    input:
        RESULTS_DIR + "/place-finish.done",
        ins=RESULTS_DIR + "/insertions_unsorted.vcf"
    output:
        temp(RESULTS_DIR + "/insertions.vcf")
    conda:
        os.path.join(WORKFLOW_PATH,"snakemodules/envs/bcftools.yml")
    resources:
        mem_mb = resources["standard_2G"]["mem"],
        runtime = resources["standard_2G"]["time"]
    threads: 
        threads["single"]
    log:
        out="logs/genotype/sort_vcf.out",
        err="logs/genotype/sort_vcf.err"
    benchmark:
        "benchmarks/genotype/sort_vcf.txt"
    shell:
        "bcftools sort {input.ins} "
        "   --output-file {output}"
        "   --output-type v"
        "   > {log.out} 2> {log.err}"


rule popins2_genotype:
    input:
        ins = rules.sort_vcf.output,
        bam = rules.coordinate_sort_unsorted.output,
        bai = rules.index_sorted.output
    output:
        temp(WORK_DIR + "/{sample}/insertions.vcf")
    resources:
        mem_mb = resources["standard_2G"]["mem"],
        runtime = resources["standard_2G"]["time"]
    threads: 
        threads["single"]
    log:
        out="logs/genotype/{sample}_genotype.out",
        err="logs/genotype/{sample}_genotype.err"
    benchmark:
        "benchmarks/genotype/{sample}_genotype.txt"
    shell:
        "{POPINS2_BIN} genotype "
        "   --prefix {WORK_DIR} "
        "   --insertions {input.ins} "
        "   --contigs {rules.popins2_merge_contigs.output.supercontigs} "
        "   --reference {REFERENCE} "
        "   {wildcards.sample} "
        "   > {log.out} 2> {log.err}"


rule merge_vcfs:
    input:
        expand(WORK_DIR + "/{s}/insertions.vcf", s=SAMPLES)
    output:
        RESULTS_DIR + "/insertions_genotypes.vcf.gz"
    conda:
        os.path.join(WORKFLOW_PATH,"snakemodules/envs/bcftools.yml")
    resources:
        mem_mb = resources["standard_2G"]["mem"],
        runtime = resources["standard_2G"]["time"]
    threads: 
        threads["single"]
    log:
        out="logs/genotype/merge_vcfs.out",
        err="logs/genotype/merge_vcfs.err"

    benchmark:
        "benchmarks/genotype/merge_vcfs.txt"
    shell:
        "bcftools merge {input}"
        "  --output {output} "
        "  --output-type z"
        "  --no-index "
        "  > {log.out} 2> {log.err}"