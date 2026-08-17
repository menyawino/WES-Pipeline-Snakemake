# ACMG/AMP 2015 Variant Annotation & Clinical Classification Rules

rule acmg_annotate_variants:
    message:
        "Classifying and annotating variants with ACMG/AMP 2015 guidelines for sample {wildcards.sample}"
    input:
        snp_vcf=rules.filter_snps.output.filtered_snp_vcf,
        indel_vcf=rules.filter_indels.output.filtered_indel_vcf
    output:
        acmg_vcf=config["outdir"] + "/analysis/007_annotation/{sample}.acmg_annotated.vcf",
        acmg_tsv=config["outdir"] + "/analysis/007_annotation/{sample}.acmg_variants.tsv",
        acmg_html=config["outdir"] + "/analysis/007_annotation/{sample}.acmg_report.html"
    conda:
        "icc_gatk"
    threads:
        config.get("threads_mid", 8)
    params:
        bed=config.get("icc_panel", "")
    log:
        config["outdir"] + "/logs/007_annotation/{sample}_acmg_annotation.log"
    benchmark:
        config["outdir"] + "/benchmarks/007_annotation/{sample}_acmg_annotation.txt"
    shell:
        """
        python3 workflow/scripts/acmg_annotator.py \
        --snp-vcf "{input.snp_vcf}" \
        --indel-vcf "{input.indel_vcf}" \
        --bed "{params.bed}" \
        --sample-name "{wildcards.sample}" \
        --output-vcf "{output.acmg_vcf}" \
        --output-tsv "{output.acmg_tsv}" \
        --output-html "{output.acmg_html}" \
        > "{log}" 2>&1
        """

rule aggregate_acmg_annotations:
    message:
        "Aggregating cohort-wide ACMG/AMP clinical variant classifications"
    input:
        tsvs=expand(config["outdir"] + "/analysis/007_annotation/{sample}.acmg_variants.tsv", sample=sample_filename)
    output:
        cohort_report=config["outdir"] + "/analysis/007_annotation/cohort_acmg_report.md",
        cohort_table=config["outdir"] + "/analysis/007_annotation/cohort_acmg_summary.tsv",
        cohort_json=config["outdir"] + "/analysis/007_annotation/cohort_acmg_summary.json",
        cohort_dashboard=config["outdir"] + "/analysis/007_annotation/cohort_acmg_dashboard.html"
    conda:
        "icc_gatk"
    log:
        config["outdir"] + "/logs/007_annotation/cohort_acmg_summary.log"
    benchmark:
        config["outdir"] + "/benchmarks/007_annotation/cohort_acmg_summary.txt"
    shell:
        """
        python3 workflow/scripts/aggregate_acmg.py \
        --inputs {input.tsvs} \
        --report-md "{output.cohort_report}" \
        --summary-tsv "{output.cohort_table}" \
        --summary-json "{output.cohort_json}" \
        --dashboard-html "{output.cohort_dashboard}" \
        > "{log}" 2>&1
        """