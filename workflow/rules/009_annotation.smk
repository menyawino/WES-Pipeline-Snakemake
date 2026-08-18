# Ensembl VEP & GeneBe ACMG 2015 Unified Variant Annotation Rules

rule vep_genebe_annotate_variants:
    message:
        "Annotating variants and classifying ACMG guidelines via Ensembl VEP & GeneBe plugin for sample {wildcards.sample}"
    input:
        snp_vcf=rules.filter_snps.output.filtered_snp_vcf,
        indel_vcf=rules.filter_indels.output.filtered_indel_vcf
    output:
        vep_vcf=config["outdir"] + "/analysis/007_annotation/{sample}.vep_annotated.vcf",
        acmg_tsv=config["outdir"] + "/analysis/007_annotation/{sample}.acmg_variants.tsv"
    conda:
        "icc_gatk"
    log:
        config["outdir"] + "/logs/007_annotation/{sample}_vep_genebe_annotation.log"
    benchmark:
        config["outdir"] + "/benchmarks/007_annotation/{sample}_vep_genebe_annotation.txt"
    shell:
        """
        python3 workflow/scripts/vep_online_annotator.py \
        --snp-vcf "{input.snp_vcf}" \
        --indel-vcf "{input.indel_vcf}" \
        --sample-name "{wildcards.sample}" \
        --output-vcf "{output.vep_vcf}" \
        --output-tsv "{output.acmg_tsv}" \
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
        cohort_json=config["outdir"] + "/analysis/007_annotation/cohort_acmg_summary.json"
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
        > "{log}" 2>&1
        """