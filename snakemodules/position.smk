

if config["remove_contamination"] == 'yes':

    rule remap_find_locations:
        input:
            expand(WORK_DIR + "/{s}/non_ref_new.bam", s=SAMPLES),
            os.path.join(RESULTS_DIR, "heatmap_coverage.png") if config["ANALYSIS"]=="yes" else [],
            WORK_DIR + "/{sample}/remapped_non_ref.bam"
        output:
            temp(WORK_DIR + "/{sample}/locations.txt")
        resources:
            mem_mb = resources["standard_4G"]["mem"],
            runtime = resources["standard_4G"]["time"]
        threads: 
            threads["single"]
        log:
            out="logs/position/{sample}_remap_find_locations.out",
            err="logs/position/{sample}_remap_find_locations.err"
        benchmark:
            "benchmarks/position/{sample}_remap_find_locations.txt"
        shell:
            "{POPINS2_BIN} find-locations {wildcards.sample}"
            "   -n remapped_non_ref.bam "
            "   --prefix {WORK_DIR} "
            "   --reference {REFERENCE} "
            "   > {log.out} 2> {log.err}"


elif config["remove_contamination"] == 'no':

    rule popins2_find_locations:
        input:
            os.path.join(RESULTS_DIR, "heatmap_coverage.png") if config["ANALYSIS"]=="yes" else [],
            expand(WORK_DIR + "/{s}/{f}", s=SAMPLES, f=["non_ref_new.bam","non_ref.bam"])
        output:
            temp(WORK_DIR + "/{sample}/locations.txt")
        resources:
            mem_mb = resources["standard_4G"]["mem"],
            runtime = resources["standard_4G"]["time"]
        threads: 
            threads["single"]
        log:
            out="logs/position/{sample}_find_locations.out",
            err="logs/position/{sample}_find_locations.err"
        benchmark:
            "benchmarks/position/{sample}_find_locations.txt"
        shell:
            "{POPINS2_BIN} find-locations {wildcards.sample}"
            "   --prefix {WORK_DIR} "
            "   --reference {REFERENCE} "
            "   > {log.out} 2> {log.err}"   

################################################################################

rule popins2_merge_locations:
    input:
        os.path.join(RESULTS_DIR, "heatmap_coverage.png") if config["ANALYSIS"]=="yes" else [],
        expand(WORK_DIR + "/{s}/locations.txt", s=SAMPLES)
    output:
        loc=RESULTS_DIR + "/locations.txt"
    resources:
        mem_mb = resources["standard_2G"]["mem"],
        runtime = resources["standard_2G"]["time"]
    threads: 
        threads["single"]
    log:
        out="logs/position/merge_locations.out",
        err="logs/position/merge_locations.err"
    benchmark:
        "benchmarks/position/merge_locations.txt"
    shell:
        "{POPINS2_BIN} merge-locations "
        "   --prefix {WORK_DIR} "
        "   --locations {RESULTS_DIR}/locations.txt "
        "   > {log.out} 2> {log.err}"


rule popins2_place_refalign:
    input:
        expand(WORK_DIR + "/{s}/non_ref_new.bam", s=SAMPLES),
        expand(WORK_DIR + "/{s}/non_ref_new.bam.bai", s=SAMPLES),
        loc=RESULTS_DIR + "/locations.txt"
    output:
        vcf=temp(RESULTS_DIR + "/insertions_unsorted.vcf"),
        grp=temp(RESULTS_DIR + "/groups.txt"),
        unp=temp(expand(WORK_DIR + "/{s}/locations_unplaced.txt", s=SAMPLES))
    params:
        readlen = config["readlen"]
    resources:
        mem_mb = resources["standard_2G"]["mem"],
        runtime = resources["standard_2G"]["time"]
    threads: 
        threads["single"]
    log:
        out="logs/position/place_refalign.out",
        err="logs/position/place_refalign.err"
    benchmark:
        "benchmarks/position/place_refalign.txt"
    shell:
        "{POPINS2_BIN} place-refalign "
        "   --prefix {WORK_DIR} "
        "   --locations {RESULTS_DIR}/locations.txt "
        "   --insertions {RESULTS_DIR}/insertions_unsorted.vcf "
        "   --groups {RESULTS_DIR}/groups.txt "
        "   --contigs {rules.popins2_merge_contigs.output.supercontigs} "
        "   --reference {REFERENCE} "
        "   --readLength {params.readlen} "
        "   > {log.out} 2> {log.err}"


rule popins2_place_splitalign:
    input:
        WORK_DIR + "/{sample}/locations_unplaced.txt"
    output:
        temp(WORK_DIR + "/{sample}/locations_placed.txt")
    params:
        readlen = config["readlen"]
    resources:
        mem_mb = resources["standard_2G"]["mem"],
        runtime = resources["standard_2G"]["time"]
    threads: 
        threads["single"]
    log:
        out="logs/position/{sample}_place_splitalign.out",
        err="logs/position/{sample}_place_splitalign.err"
    benchmark:
        "benchmarks/position/{sample}_place_splitalign.txt"
    shell:
        "{POPINS2_BIN} place-splitalign "
        "   --prefix {WORK_DIR} "
        "   --contigs {rules.popins2_merge_contigs.output.supercontigs} "
        "   --reference {REFERENCE} "
        "   --readLength {params.readlen} "
        "   {wildcards.sample} "
        "   > {log.out} 2> {log.err}"


rule popins2_place_finish:
    input:
        expand(WORK_DIR + "/{s}/locations_placed.txt", s=SAMPLES)
    output:
        temp(touch(RESULTS_DIR + "/place-finish.done"))
    resources:
        mem_mb = resources["standard_2G"]["mem"],
        runtime = resources["standard_2G"]["time"]
    threads: 
        threads["single"]
    log:
        out="logs/position/place_finish.out",
        err="logs/position/place_finish.err"
    benchmark:
        "benchmarks/position/place_finish.txt"
    shell:
        "{POPINS2_BIN} place-finish "
        "   --prefix {WORK_DIR} "
        "   --insertions {rules.popins2_place_refalign.output.vcf}"
        "   --reference {REFERENCE} "
        "   > {log.out} 2> {log.err}"