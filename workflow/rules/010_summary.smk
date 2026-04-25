# ClawBio-inspired variant summaries built directly from filtered SNP/indel VCFs.

rule summarize_variants:
    message:
        "Generating variant summary for sample {wildcards.sample}"
    input:
        snp_vcf=rules.filter_snps.output.filtered_snp_vcf,
        indel_vcf=rules.filter_indels.output.filtered_indel_vcf
    output:
        report_md=config["outdir"] + "/analysis/010_summary/{sample}/variant_summary.md",
        summary_tsv=config["outdir"] + "/analysis/010_summary/{sample}/variant_summary.tsv",
        summary_json=config["outdir"] + "/analysis/010_summary/{sample}/variant_summary.json"
    conda:
        "envs/010_summary.yml"
    params:
        top_variants=config.get("variant_summary", {}).get("top_variants", 10),
        top_chromosomes=config.get("variant_summary", {}).get("top_chromosomes", 10)
    log:
        config["outdir"] + "/logs/010_summary/{sample}_variant_summary.log"
    benchmark:
        config["outdir"] + "/benchmarks/010_summary/{sample}_variant_summary.txt"
    shell:
        """
        python workflow/scripts/variant_summary.py sample \
        --sample-name "{wildcards.sample}" \
        --snp-vcf {input.snp_vcf} \
        --indel-vcf {input.indel_vcf} \
        --report-md {output.report_md} \
        --summary-tsv {output.summary_tsv} \
        --summary-json {output.summary_json} \
        --top-variants {params.top_variants} \
        --top-chromosomes {params.top_chromosomes} \
        > {log} 2>&1
        """


rule aggregate_variant_summaries:
    message:
        "Aggregating cohort variant summaries"
    input:
        reports=expand(config["outdir"] + "/analysis/010_summary/{sample}/variant_summary.json", sample=sample_filename)
    output:
        cohort_report=config["outdir"] + "/analysis/010_summary/cohort_variant_report.md",
        cohort_table=config["outdir"] + "/analysis/010_summary/cohort_variant_summary.tsv",
        cohort_json=config["outdir"] + "/analysis/010_summary/cohort_variant_summary.json"
    conda:
        "envs/010_summary.yml"
    log:
        config["outdir"] + "/logs/010_summary/cohort_variant_summary.log"
    benchmark:
        config["outdir"] + "/benchmarks/010_summary/cohort_variant_summary.txt"
    shell:
        """
        python workflow/scripts/variant_summary.py cohort \
        --inputs {input.reports} \
        --report-md {output.cohort_report} \
        --summary-tsv {output.cohort_table} \
        --summary-json {output.cohort_json} \
        > {log} 2>&1
        """
