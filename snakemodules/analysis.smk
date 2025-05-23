import os

if config["SICKLE"]=="no":
    rule email_unmapped_analysis:
        input:
            single=expand("{p}/{s}/single.fastq", p=WORK_DIR, s=SAMPLES),
            paired1=expand("{p}/{s}/paired.1.fastq", p=WORK_DIR, s=SAMPLES),
            paired2=expand("{p}/{s}/paired.2.fastq", p=WORK_DIR, s=SAMPLES)
        output:
            notify_file=temp("notify/unmapped_analysis.txt")
        params:
            email=config["email"],  # Replace with your email or retrieve from config
            subject="Snakemake Notification",
            message="Rule unmapped_analysis is ready to execute.",
            receive_email=config["receive_email_notifications"]  # "yes" or "no"
        resources:
            mem_mb = resources["email"]["mem"],
            runtime = resources["email"]["time"]
        threads: threads["single"]
        group: "unmapped_analysis"
        shell:
            """
            if [ "{params.receive_email}" = "yes" ]; then
                echo "{params.message}" | mail -s "{params.subject}" {params.email}
            fi
            touch {output.notify_file}
            """
    rule unmapped_analysis:
        input:
            notify="notify/unmapped_analysis.txt",
            single=expand("{p}/{s}/single.fastq", p=WORK_DIR, s=SAMPLES),
            paired1=expand("{p}/{s}/paired.1.fastq", p=WORK_DIR, s=SAMPLES),
            paired2=expand("{p}/{s}/paired.2.fastq", p=WORK_DIR, s=SAMPLES)
        output:
            RESULTS_DIR + "/read_numbers.pdf",
            RESULTS_DIR + "/read_numbers.txt"
        conda:
            os.path.join(WORKFLOW_PATH,"snakemodules/envs/notebooks.yml")
        benchmark:
            "benchmarks/analysis/unmapped_analysis.txt"
        resources:
            mem_mb = resources["analysis"]["mem"],
            runtime = resources["analysis"]["time"]
        threads: threads["single"]
        group: "unmapped_analysis"
        container:
            containers["popins4snake"]
        notebook:
            os.path.join(WORKFLOW_PATH,"snakemodules/notebooks/read_numbers_in_bar_graph.py.ipynb")



        

    if config["remove_contamination"]=="yes":
        rule email_unmapped_analysis_no_filtler_after_clean:
            input:
                single=expand("{p}/{s}/single_clean.fastq", p=WORK_DIR, s=SAMPLES),
                paired1=expand("{p}/{s}/paired_1_clean.fastq", p=WORK_DIR, s=SAMPLES),
                paired2=expand("{p}/{s}/paired_2_clean.fastq", p=WORK_DIR, s=SAMPLES)
            output:
                notify_file=temp("notify/unmapped_analysis_no_filtler_after_clean.txt")
            params:
                email=config["email"],  # Replace with your email or retrieve from config
                subject="Snakemake Notification",
                message="Rule unmapped_analysis_no_filtler_after_clean is ready to execute.",
                receive_email=config["receive_email_notifications"]  # "yes" or "no"
            resources:
                mem_mb = resources["email"]["mem"],
                runtime = resources["email"]["time"]
            threads: threads["single"]
            group: "unmapped_analysis_no_filtler_after_clean"
            shell:
                """
                if [ "{params.receive_email}" = "yes" ]; then
                    echo "{params.message}" | mail -s "{params.subject}" {params.email}
                fi
                touch {output.notify_file}
                """
            
        rule unmapped_analysis_no_filtler_after_clean:
            input:
                notify="notify/unmapped_analysis_no_filtler_after_clean.txt",
                single=expand("{p}/{s}/single_clean.fastq", p=WORK_DIR, s=SAMPLES),
                paired1=expand("{p}/{s}/paired_1_clean.fastq", p=WORK_DIR, s=SAMPLES),
                paired2=expand("{p}/{s}/paired_2_clean.fastq", p=WORK_DIR, s=SAMPLES)
            output:
                RESULTS_DIR + "/read_numbers_no_sickle_after_clean.pdf",
                RESULTS_DIR + "/read_numbers_no_sickle_after_clean.txt"
            conda:
                os.path.join(WORKFLOW_PATH,"snakemodules/envs/notebooks.yml")
            benchmark:
                "benchmarks/analysis/unmapped_analysis_after_cleaning.txt"
            resources:
                mem_mb = resources["analysis"]["mem"],
                runtime = resources["analysis"]["time"]
            threads: threads["single"]
            group: "unmapped_analysis_no_filtler_after_clean"
            container:
                containers["popins4snake"]
            notebook:
                os.path.join(WORKFLOW_PATH,"snakemodules/notebooks/read_numbers_in_bar_graph.py.ipynb")

            

elif config["SICKLE"]=="yes":
    if config["remove_contamination"]=="no":
        rule email_unmapped_analysis_before_filter:
            input:
                single=expand("{p}/{s}/single.fastq", p=WORK_DIR, s=SAMPLES),
                paired1=expand("{p}/{s}/paired.1.fastq", p=WORK_DIR, s=SAMPLES),
                paired2=expand("{p}/{s}/paired.2.fastq", p=WORK_DIR, s=SAMPLES)
            output:
                notify_file=temp("notify/unmapped_analysis_before_filter.txt")
            params:
                email=config["email"],  # Replace with your email or retrieve from config
                subject="Snakemake Notification",
                message="Rule unmapped_analysis_before_filter is ready to execute.",
                receive_email=config["receive_email_notifications"]  # "yes" or "no"
            resources:
                mem_mb = resources["email"]["mem"],
                runtime = resources["email"]["time"]
            threads: threads["single"]
            group: "unmapped_analysis_before_filter"
            shell:
                """
                if [ "{params.receive_email}" = "yes" ]; then
                    echo "{params.message}" | mail -s "{params.subject}" {params.email}
                fi
                touch {output.notify_file}
                """
            
        rule unmapped_analysis_before_filter:
            input:
                notify="notify/unmapped_analysis_before_filter.txt",
                single=expand("{p}/{s}/single.fastq", p=WORK_DIR, s=SAMPLES),
                paired1=expand("{p}/{s}/paired.1.fastq", p=WORK_DIR, s=SAMPLES),
                paired2=expand("{p}/{s}/paired.2.fastq", p=WORK_DIR, s=SAMPLES)
            output:
                RESULTS_DIR + "/read_numbers_before_filter.pdf",
                RESULTS_DIR + "/read_numbers_before_filter.txt",
            conda:
                os.path.join(WORKFLOW_PATH,"snakemodules/envs/notebooks.yml")
            benchmark:
                "benchmarks/analysis/unmapped_analysis_before_filter.txt"
            resources:
                mem_mb = resources["analysis"]["mem"],
                runtime = resources["analysis"]["time"]
            threads: threads["single"]
            group: "unmapped_analysis_before_filter"
            container:
                containers["popins4snake"]
            notebook:
                os.path.join(WORKFLOW_PATH,"snakemodules/notebooks/read_numbers_in_bar_graph.py.ipynb")


            

        rule email_unmapped_analysis_after_filter:
            input:
                single = expand("{p}/{s}/sickle.single.fastq", p=WORK_DIR, s=SAMPLES),
                paired1 = expand("{p}/{s}/sickle.paired.1.fastq", p=WORK_DIR, s=SAMPLES),
                paired2 = expand("{p}/{s}/sickle.paired.2.fastq", p=WORK_DIR, s=SAMPLES)
            output:
                notify_file=temp("notify/unmapped_analysis_after_filter.txt")
            params:
                email=config["email"],  # Replace with your email or retrieve from config
                subject="Snakemake Notification",
                message="Rule unmapped_analysis_after_filter is ready to execute.",
                receive_email=config["receive_email_notifications"]  # "yes" or "no"
            resources:
                mem_mb = resources["email"]["mem"],
                runtime = resources["email"]["time"]
            threads: threads["single"]
            group: "unmapped_analysis_after_filter"
            shell:
                """
                if [ "{params.receive_email}" = "yes" ]; then
                    echo "{params.message}" | mail -s "{params.subject}" {params.email}
                fi
                touch {output.notify_file}
                """
            
        rule unmapped_analysis_after_filter:
            input:
                notify="notify/unmapped_analysis_after_filter.txt",
                single = expand("{p}/{s}/sickle.single.fastq", p=WORK_DIR, s=SAMPLES),
                paired1 = expand("{p}/{s}/sickle.paired.1.fastq", p=WORK_DIR, s=SAMPLES),
                paired2 = expand("{p}/{s}/sickle.paired.2.fastq", p=WORK_DIR, s=SAMPLES)
            output:
                RESULTS_DIR + "/read_numbers_after_filter.pdf",
                RESULTS_DIR + "/read_numbers_after_filter.txt"
            conda:
                os.path.join(WORKFLOW_PATH,"snakemodules/envs/notebooks.yml")
            benchmark:
                "benchmarks/analysis/unmapped_analysis_after_filter.txt"
            resources:
                mem_mb = resources["analysis"]["mem"],
                runtime = resources["analysis"]["time"]
            threads: threads["single"]
            group: "unmapped_analysis_after_filter"
            container:
                containers["popins4snake"]
            notebook:
                os.path.join(WORKFLOW_PATH,"snakemodules/notebooks/read_numbers_in_bar_graph.py.ipynb")


            

    elif config["remove_contamination"]=="yes":
        rule email_unmapped_analysis_before_clean_filter:
            input:
                single = expand("{p}/{s}/single.fastq", p=WORK_DIR, s=SAMPLES),
                paired1 = expand("{p}/{s}/paired.1.fastq", p=WORK_DIR, s=SAMPLES),
                paired2 = expand("{p}/{s}/paired.2.fastq", p=WORK_DIR, s=SAMPLES)
            output:
                notify_file=temp("notify/unmapped_analysis_before_clean_filter.txt")
            params:
                email=config["email"],  # Replace with your email or retrieve from config
                subject="Snakemake Notification",
                message="Rule unmapped_analysis_before_clean_filter is ready to execute.",
                receive_email=config["receive_email_notifications"]  # "yes" or "no"
            resources:
                mem_mb = resources["email"]["mem"],
                runtime = resources["email"]["time"]
            threads: threads["single"]
            group: "unmapped_analysis_before_clean_filter"
            shell:
                """
                if [ "{params.receive_email}" = "yes" ]; then
                    echo "{params.message}" | mail -s "{params.subject}" {params.email}
                fi
                touch {output.notify_file}
                """
            
        rule unmapped_analysis_before_clean_filter:
            input:
                notify="notify/unmapped_analysis_before_clean_filter.txt",
                single = expand("{p}/{s}/single.fastq", p=WORK_DIR, s=SAMPLES),
                paired1 = expand("{p}/{s}/paired.1.fastq", p=WORK_DIR, s=SAMPLES),
                paired2 = expand("{p}/{s}/paired.2.fastq", p=WORK_DIR, s=SAMPLES)
            output:
                RESULTS_DIR + "/read_numbers_before_clean_filter.pdf",
                RESULTS_DIR + "/read_numbers_before_clean_filter.txt"
            conda:
                os.path.join(WORKFLOW_PATH,"snakemodules/envs/notebooks.yml")
            benchmark:
                "benchmarks/analysis/unmapped_analysis_before_clean_filte.txt"
            resources:
                mem_mb = resources["analysis"]["mem"],
                runtime = resources["analysis"]["time"]
            threads: threads["single"]
            group: "unmapped_analysis_before_clean_filter"
            container:
                containers["popins4snake"]
            notebook:
                os.path.join(WORKFLOW_PATH,"snakemodules/notebooks/read_numbers_in_bar_graph.py.ipynb")

            

        rule email_unmapped_analysis_after_clean:
            input:
                single = expand("{p}/{s}/single_clean.fastq", p=WORK_DIR, s=SAMPLES),
                paired1 = expand("{p}/{s}/paired_1_clean.fastq", p=WORK_DIR, s=SAMPLES),
                paired2 = expand("{p}/{s}/paired_2_clean.fastq", p=WORK_DIR, s=SAMPLES)
            output:
                notify_file=temp("notify/unmapped_analysis_after_clean.txt")
            params:
                email=config["email"],  # Replace with your email or retrieve from config
                subject="Snakemake Notification",
                message="Rule unmapped_analysis_after_clean is ready to execute.",
                receive_email=config["receive_email_notifications"]  # "yes" or "no"
            resources:
                mem_mb = resources["email"]["mem"],
                runtime = resources["email"]["time"]
            threads: threads["single"]
            group: "unmapped_analysis_after_clean"
            shell:
                """
                if [ "{params.receive_email}" = "yes" ]; then
                    echo "{params.message}" | mail -s "{params.subject}" {params.email}
                fi
                touch {output.notify_file}
                """
            
        rule unmapped_analysis_after_clean:
            input:
                notify="notify/unmapped_analysis_after_clean.txt",
                single = expand("{p}/{s}/single_clean.fastq", p=WORK_DIR, s=SAMPLES),
                paired1 = expand("{p}/{s}/paired_1_clean.fastq", p=WORK_DIR, s=SAMPLES),
                paired2 = expand("{p}/{s}/paired_2_clean.fastq", p=WORK_DIR, s=SAMPLES)
            output:
                RESULTS_DIR + "/read_numbers_after_clean.pdf",
                RESULTS_DIR + "/read_numbers_after_clean.txt"
            conda:
                os.path.join(WORKFLOW_PATH,"snakemodules/envs/notebooks.yml")
            benchmark:
                "benchmarks/analysis/unmapped_analysis_after_cleaning.txt"
            resources:
                mem_mb = resources["analysis"]["mem"],
                runtime = resources["analysis"]["time"]
            threads: threads["single"]
            group: "unmapped_analysis_after_clean"
            container:
                containers["popins4snake"]
            notebook:
                os.path.join(WORKFLOW_PATH,"snakemodules/notebooks/read_numbers_in_bar_graph.py.ipynb")
            

        rule email_unmapped_analysis_after_clean_filter:
            input:
                single = expand("{p}/{s}/sickle.single_clean.fastq", p=WORK_DIR, s=SAMPLES),
                paired1 = expand("{p}/{s}/sickle.paired_1_clean.fastq", p=WORK_DIR, s=SAMPLES),
                paired2 = expand("{p}/{s}/sickle.paired_2_clean.fastq", p=WORK_DIR, s=SAMPLES)
            output:
                notify_file=temp("notify/unmapped_analysis_after_clean_filter.txt")
            params:
                email=config["email"],  # Replace with your email or retrieve from config
                subject="Snakemake Notification",
                message="Rule unmapped_analysis_after_clean_filter is ready to execute.",
                receive_email=config["receive_email_notifications"]  # "yes" or "no"
            resources:
                mem_mb = resources["email"]["mem"],
                runtime = resources["email"]["time"]
            threads: threads["single"]
            group: "unmapped_analysis_after_clean_filter"
            shell:
                """
                if [ "{params.receive_email}" = "yes" ]; then
                    echo "{params.message}" | mail -s "{params.subject}" {params.email}
                fi
                touch {output.notify_file}
                """
            
        rule unmapped_analysis_after_clean_filter:
            input:
                notify="notify/unmapped_analysis_after_clean_filter.txt",
                single = expand("{p}/{s}/sickle.single_clean.fastq", p=WORK_DIR, s=SAMPLES),
                paired1 = expand("{p}/{s}/sickle.paired_1_clean.fastq", p=WORK_DIR, s=SAMPLES),
                paired2 = expand("{p}/{s}/sickle.paired_2_clean.fastq", p=WORK_DIR, s=SAMPLES)
            output:
                RESULTS_DIR + "/read_numbers_after_clean_filter.pdf",
                RESULTS_DIR + "/read_numbers_after_clean_filter.txt"
            conda:
                os.path.join(WORKFLOW_PATH,"snakemodules/envs/notebooks.yml")
            benchmark:
                "benchmarks/analysis/unmapped_analysis_after_clean_filter.txt"
            resources:
                mem_mb = resources["analysis"]["mem"],
                runtime = resources["analysis"]["time"]
            threads: threads["single"]
            group: "unmapped_analysis_after_clean_filter"
            container:
                containers["popins4snake"]
            notebook:
                os.path.join(WORKFLOW_PATH,"snakemodules/notebooks/read_numbers_in_bar_graph.py.ipynb")

            


rule email_assembly_analysis_table:
    input:
        contigs = expand("{p}/{s}/{assemblr}.contigs.fa", p=WORK_DIR, s=SAMPLES, assemblr=ASSEMBLER)
    output:
        notify_file=temp("notify/assembly_analysis_table.txt")
    params:
        email=config["email"],  # Replace with your email or retrieve from config
        subject="Snakemake Notification",
        message="Rule assembly_analysis_table is ready to execute.",
        receive_email=config["receive_email_notifications"]  # "yes" or "no"
    resources:
        mem_mb = resources["email"]["mem"],
        runtime = resources["email"]["time"]
    threads: threads["single"]
    group: "assembly_analysis_table"
    shell:
        """
        if [ "{params.receive_email}" = "yes" ]; then
            echo "{params.message}" | mail -s "{params.subject}" {params.email}
        fi
        touch {output.notify_file}
        """
    
rule assembly_analysis_table:
    input:
        notify="notify/assembly_analysis_table.txt",
        contigs = expand("{p}/{s}/{assemblr}.contigs.fa", p=WORK_DIR, s=SAMPLES, assemblr=ASSEMBLER)
    output:
        RESULTS_DIR + "/analysis_table.pdf",
        RESULTS_DIR + "/analysis_table.txt"
    conda:
        os.path.join(WORKFLOW_PATH,"snakemodules/envs/notebooks.yml")
    benchmark:
        "benchmarks/analysis/assembly_analysis_table.txt"
    resources:
        mem_mb = resources["analysis"]["mem"],
        runtime = resources["analysis"]["time"]
    threads: threads["single"]
    group: "assembly_analysis_table"
    container:
        containers["popins4snake"]
    notebook:
        os.path.join(WORKFLOW_PATH,"snakemodules/notebooks/sample_analysis_table.py.ipynb")
    
rule email_assembly_analysis:
    input:
        contigs = expand("{p}/{s}/{assemblr}.contigs.fa", p=WORK_DIR, s=SAMPLES, assemblr=ASSEMBLER)
    output:
        notify_file=temp("notify/assembly_analysis.txt")
    params:
        email=config["email"],  # Replace with your email or retrieve from config
        subject="Snakemake Notification",
        message="Rule assembly_analysis is ready to execute.",
        receive_email=config["receive_email_notifications"]  # "yes" or "no"
    resources:
        mem_mb = resources["email"]["mem"],
        runtime = resources["email"]["time"]
    threads: threads["single"]
    group: "assembly_analysis"
    shell:
        """
        if [ "{params.receive_email}" = "yes" ]; then
            echo "{params.message}" | mail -s "{params.subject}" {params.email}
        fi
        touch {output.notify_file}
        """
    
rule assembly_analysis:
    input:
        notify="notify/assembly_analysis.txt",
        contigs = expand("{p}/{s}/{assemblr}.contigs.fa", p=WORK_DIR, s=SAMPLES, assemblr=ASSEMBLER)
    output:
        len = RESULTS_DIR + "/contig_length_dist.pdf",
        gc = RESULTS_DIR + "/gc_content_dist.pdf",
        len_txt = RESULTS_DIR + "/distribution_length.txt",
        gc_txt = RESULTS_DIR + "/distribution_gc.txt",
        distribution = RESULTS_DIR + "/sample_distribution_info.txt"
    conda:
        os.path.join(WORKFLOW_PATH,"snakemodules/envs/notebooks.yml")
    benchmark:
        "benchmarks/analysis/assembly_analysis.txt"
    resources:
        mem_mb = resources["analysis"]["mem"],
        runtime = resources["analysis"]["time"]
    threads: threads["single"]
    group: "assembly_analysis"
    container:
        containers["popins4snake"]
    notebook:
        os.path.join(WORKFLOW_PATH,"snakemodules/notebooks/content_distribution.py.ipynb")
    

    
rule calculate_coverage:
    input:
        bams = WORK_DIR + "/{sample}/non_ref_new.bam",
        bais = WORK_DIR + "/{sample}/non_ref_new.bam.bai"
    output:
        cov = temp(WORK_DIR + "/{sample}/coverage.txt")
    conda:
        os.path.join(WORKFLOW_PATH,"snakemodules/envs/samtools21.yml")
    benchmark:
        "benchmarks/analysis/{sample}_calculate_coverage.txt"
    resources:
        mem_mb = resources["analysis"]["mem"],
        runtime = resources["analysis"]["time"]
    threads: threads["single"]
    container:
        containers["popins4snake"]
    shell:
        "samtools coverage {input.bams} > {output.cov}"


rule email_coverage_analysis:
    input:
        cov = expand(WORK_DIR + "/{s}/coverage.txt", s=SAMPLES),
        fasta = RESULTS_DIR + "/supercontigs.fa"
    output:
        notify_file=temp("notify/coverage_analysis.txt")
    params:
        email=config["email"],  # Replace with your email or retrieve from config
        subject="Snakemake Notification",
        message="Rule coverage_analysis is ready to execute.",
        receive_email=config["receive_email_notifications"]  # "yes" or "no"
    resources:
        mem_mb = resources["email"]["mem"],
        runtime = resources["email"]["time"]
    threads: threads["single"]
    group: "coverage_analysis"
    shell:
        """
        if [ "{params.receive_email}" = "yes" ]; then
            echo "{params.message}" | mail -s "{params.subject}" {params.email}
        fi
        touch {output.notify_file}
        """
    
rule coverage_analysis:
    input:
        notify="notify/coverage_analysis.txt",
        cov = expand(WORK_DIR + "/{s}/coverage.txt", s=SAMPLES),
        fasta = RESULTS_DIR + "/supercontigs.fa"
    output:
        RESULTS_DIR + "/heatmap_coverage.pdf",
        RESULTS_DIR + "/heatmap_coverage.txt",
        RESULTS_DIR + "/contigs_sorted_by_length.txt"
    conda:
        os.path.join(WORKFLOW_PATH,"snakemodules/envs/notebooks.yml")
    benchmark:
        "benchmarks/analysis/coverage_analysis.txt"
    resources:
        mem_mb = resources["analysis"]["mem"],
        runtime = resources["analysis"]["time"]
    threads: threads["single"]
    group: "coverage_analysis"
    container:
        containers["popins4snake"]
    notebook:
        os.path.join(WORKFLOW_PATH,"snakemodules/notebooks/contig_coverage_heatmap.py.ipynb")
