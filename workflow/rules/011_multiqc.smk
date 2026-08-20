# MultiQC rule to aggregate all QC, alignment, variant calling, and annotation metrics

def get_multiqc_inputs(wildcards):
    items = [
        expand(rules.trimming_fp.output.json, sample=sample_filename, lane=lane),
        expand(rules.qc_report.output.qc_metrics, sample=sample_filename),
        expand(rules.vep_genebe_annotate_variants.output.vep_vcf, sample=sample_filename, caller=CALLERS),
        expand(rules.aggregate_variant_summaries.output.cohort_report, caller=CALLERS),
        expand(rules.aggregate_acmg_annotations.output.cohort_report, caller=CALLERS),
    ]
    if enable_comparison:
        items.append([rules.compare_callers.output.report])
    
    flat = []
    for x in items:
        if isinstance(x, list):
            flat.extend(x)
        else:
            flat.append(x)
    return flat

rule multiqc_report:
    message:
        "Generating final MultiQC report across all pipeline outputs"
    input:
        get_multiqc_inputs
    output:
        html=config["outdir"] + "/results/multiqc_report.html"
    container:
        "docker://quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0"
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
