
rule popins2_merge_contigs:
    input:
        expand("{p}/{s}/{assemblr}.contigs.fa", p=WORK_DIR, s=SAMPLES, assemblr=ASSEMBLER),
        os.path.join(RESULTS_DIR, "analysis_table.pdf") if config["ANALYSIS"]=="yes" else [],
        os.path.join(RESULTS_DIR, "contig_length_dist.pdf") if config["ANALYSIS"]=="yes" else [],
        os.path.join(RESULTS_DIR, "gc_content_dist.pdf") if config["ANALYSIS"]=="yes" else []
    output:
        supercontigs = RESULTS_DIR + "/supercontigs.fa",
        gfa          = RESULTS_DIR + "/supercontigs.gfa",
        colors       = RESULTS_DIR + "/supercontigs.bfg_colors"
    resources:
        mem_mb = resources["standard_2G"]["mem"],
        runtime = resources["standard_2G"]["time"]
    threads: 
        threads["multi"]["merge"]
    params:
        k = config["k_merge"]
    log:
        out="logs/merge/merge_contigs.out",
        err="logs/merge/merge_contigs.err"
    benchmark:
        "benchmarks/merge/merge_contigs.txt"
    shell:
        "{POPINS2_BIN} merge-contigs"
        "   -k {params.k}"
        "   -f {ASSEMBLER}.contigs.fa"
        "   -t {threads}"
        "   -r {WORK_DIR}"
        "   -p {RESULTS_DIR}/supercontigs"
        "   -di"
        "   > {log.out} 2> {log.err}"
