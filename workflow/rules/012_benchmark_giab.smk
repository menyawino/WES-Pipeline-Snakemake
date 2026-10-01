# Genome in a Bottle (GIAB) Benchmarking Rules

rule prepare_giab_truth:
    message:
        "Preparing GIAB {wildcards.sample} truth resources (GRCh38)"
    output:
        truth_vcf="resources/giab/{sample}/{sample}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        truth_tbi="resources/giab/{sample}/{sample}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi",
        truth_bed="resources/giab/{sample}/{sample}_GRCh38_1_22_v4.2.1_benchmark.bed"
    log:
        "resources/giab/{sample}/download.log"
    shell:
        """
        python3 workflow/scripts/download_giab_truth.py \
            --sample "{wildcards.sample}" \
            --outdir "resources/giab/{wildcards.sample}" > "{log}" 2>&1
        """

rule benchmark_sample_giab:
    message:
        "Benchmarking sample {wildcards.sample} variants against GIAB truth set"
    input:
        truth_vcf=lambda wildcards: config.get("giab_benchmark", {}).get("truth_vcf", f"resources/giab/{config.get('giab_benchmark', {}).get('sample', 'HG001')}/{config.get('giab_benchmark', {}).get('sample', 'HG001')}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"),
        truth_bed=lambda wildcards: config.get("giab_benchmark", {}).get("truth_bed", f"resources/giab/{config.get('giab_benchmark', {}).get('sample', 'HG001')}/{config.get('giab_benchmark', {}).get('sample', 'HG001')}_GRCh38_1_22_v4.2.1_benchmark.bed"),
        target_bed=config.get("icc_panel"),
        snp_vcfs=expand(config["outdir"] + "/analysis/006_variant_filtering/{{sample}}.{caller}.filtered.snp.vcf", caller=CALLERS),
        indel_vcfs=expand(config["outdir"] + "/analysis/006_variant_filtering/{{sample}}.{caller}.filtered.indel.vcf", caller=CALLERS)
    output:
        report_md=config["outdir"] + "/analysis/012_benchmark/giab/{sample}_giab_benchmark_report.md",
        summary_tsv=config["outdir"] + "/analysis/012_benchmark/giab/{sample}_giab_benchmark_summary.tsv",
        summary_json=config["outdir"] + "/analysis/012_benchmark/giab/{sample}_giab_benchmark_summary.json",
        dashboard_html=config["outdir"] + "/analysis/012_benchmark/giab/{sample}_giab_benchmark_dashboard.html"
    conda:
        "../envs/005_gatk_genomics.yml"
    threads:
        config.get("threads_low", 2)
    resources:
        mem_mb=config.get("mem_low", 4096),
        tmpdir=config.get("tmpdir", "/tmp")
    log:
        config["outdir"] + "/logs/012_benchmark/giab/{sample}_benchmark.log"
    benchmark:
        config["outdir"] + "/benchmarks/012_benchmark/giab/{sample}_benchmark.txt"
    shell:
        """
        python3 workflow/scripts/benchmark_giab.py \
            --outdir "{config[outdir]}" \
            --sample "{wildcards.sample}" \
            --callers {" ".join(CALLERS)} \
            --truth-vcf "{input.truth_vcf}" \
            --truth-bed "{input.truth_bed}" \
            --target-bed "{input.target_bed}" \
            --report-md "{output.report_md}" \
            --summary-tsv "{output.summary_tsv}" \
            --summary-json "{output.summary_json}" \
            --dashboard-html "{output.dashboard_html}" > "{log}" 2>&1
        """
