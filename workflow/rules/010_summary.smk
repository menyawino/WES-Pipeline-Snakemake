# ClawBio-inspired variant summaries built directly from filtered SNP/indel VCFs.

rule summarize_variants:
    message:
        "Generating {wildcards.caller} variant summary for sample {wildcards.sample}"
    input:
        snp_vcf=rules.filter_snps.output.filtered_snp_vcf,
        indel_vcf=rules.filter_indels.output.filtered_indel_vcf,
        annotated_vcf=rules.vep_genebe_annotate_variants.output.vep_vcf
    output:
        report_md=config["outdir"] + "/analysis/010_summary/{sample}/{caller}_variant_summary.md",
        summary_tsv=config["outdir"] + "/analysis/010_summary/{sample}/{caller}_variant_summary.tsv",
        summary_json=config["outdir"] + "/analysis/010_summary/{sample}/{caller}_variant_summary.json"
    conda:
        "icc_gatk"
    params:
        top_variants=config.get("variant_summary", {}).get("top_variants", 10),
        top_chromosomes=config.get("variant_summary", {}).get("top_chromosomes", 10)
    log:
        config["outdir"] + "/logs/010_summary/{sample}_{caller}_variant_summary.log"
    benchmark:
        config["outdir"] + "/benchmarks/010_summary/{sample}_{caller}_variant_summary.txt"
    shell:
        """
        python workflow/scripts/variant_summary.py sample \
        --sample-name "{wildcards.sample}_{wildcards.caller}" \
        --snp-vcf "{input.snp_vcf}" \
        --indel-vcf "{input.indel_vcf}" \
        --annotated-vcf "{input.annotated_vcf}" \
        --report-md "{output.report_md}" \
        --summary-tsv "{output.summary_tsv}" \
        --summary-json "{output.summary_json}" \
        --top-variants {params.top_variants} \
        --top-chromosomes {params.top_chromosomes} \
        > "{log}" 2>&1
        """

rule aggregate_variant_summaries:
    message:
        "Aggregating cohort variant summaries for {wildcards.caller}"
    input:
        reports=lambda wildcards: expand(config["outdir"] + "/analysis/010_summary/{sample}/" + wildcards.caller + "_variant_summary.json", sample=sample_filename)
    output:
        cohort_report=config["outdir"] + "/analysis/010_summary/cohort_{caller}_variant_report.md",
        cohort_table=config["outdir"] + "/analysis/010_summary/cohort_{caller}_variant_summary.tsv",
        cohort_json=config["outdir"] + "/analysis/010_summary/cohort_{caller}_variant_summary.json"
    conda:
        "icc_gatk"
    log:
        config["outdir"] + "/logs/010_summary/cohort_{caller}_variant_summary.log"
    benchmark:
        config["outdir"] + "/benchmarks/010_summary/cohort_{caller}_variant_summary.txt"
    shell:
        """
        python workflow/scripts/variant_summary.py cohort \
        --inputs {input.reports} \
        --report-md "{output.cohort_report}" \
        --summary-tsv "{output.cohort_table}" \
        --summary-json "{output.cohort_json}" \
        > "{log}" 2>&1
        """

rule compare_callers:
    message:
        "Evaluating concordance and discordance between GATK and DeepVariant calls"
    input:
        gatk_snps=expand(config["outdir"] + "/analysis/006_variant_filtering/{sample}.gatk.filtered.snp.vcf", sample=sample_filename),
        dv_snps=expand(config["outdir"] + "/analysis/006_variant_filtering/{sample}.deepvariant.filtered.snp.vcf", sample=sample_filename)
    output:
        report=config["outdir"] + "/analysis/010_summary/caller_concordance_report.md",
        table=config["outdir"] + "/analysis/010_summary/caller_concordance_matrix.tsv"
    conda:
        "icc_gatk"
    params:
        outdir=config["outdir"],
        samples=sample_filename
    log:
        config["outdir"] + "/logs/010_summary/caller_concordance.log"
    benchmark:
        config["outdir"] + "/benchmarks/010_summary/caller_concordance.txt"
    shell:
        """
        python3 workflow/scripts/compare_callers.py \
        --outdir "{params.outdir}" \
        --samples {params.samples} \
        --report-md "{output.report}" \
        --summary-tsv "{output.table}" \
        > "{log}" 2>&1
        """
