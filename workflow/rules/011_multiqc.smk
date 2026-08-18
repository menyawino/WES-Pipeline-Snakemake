# MultiQC rule to aggregate all QC, alignment, variant calling, and annotation metrics

rule multiqc_report:
    message:
        "Generating final MultiQC report across all pipeline outputs"
    input:
        expand(rules.raw_fastqc.output.html, sample_filename=sample_filename, lane=lane, R=read),
        expand(rules.posttrim_fastqc.output.html, sample=sample_filename, lane=lane, R=read),
        expand(rules.qc_report.output.qc_metrics, sample=sample_filename),
        expand(rules.vep_genebe_annotate_variants.output.vep_vcf, sample=sample_filename),
        rules.aggregate_variant_summaries.output.cohort_report,
        rules.aggregate_acmg_annotations.output.cohort_report
    output:
        html=config["outdir"] + "/results/multiqc_report.html"
    conda:
        "multiqc_env"
    log:
        config["outdir"] + "/logs/multiqc.log"
    benchmark:
        config["outdir"] + "/benchmarks/multiqc.txt"
    shell:
        """
        multiqc "{config[outdir]}" \
        -o "{config[outdir]}/results" \
        -n "multiqc_report.html" \
        --force \
        > "{log}" 2>&1
        """
