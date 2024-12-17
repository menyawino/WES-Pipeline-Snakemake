# A rule to annotate variants using VEP

rule annotate_variants:
    message:
        "Annotating variants for sample {wildcards.sample}_{lane}"
    input:
        vcf=lambda wildcards: config["outdir"] + "/analysis/006_variant_filtering/{sample}_{lane}/{sample}_{lane}.filtered.vcf".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        annotated_vcf=lambda wildcards: config["outdir"] + "/analysis/007_annotation/{sample}_{lane}/{sample}_{lane}.annotated.vcf".format(sample=wildcards.sample, lane=wildcards.lane)
    conda:
        "icc_07_annotation"
    threads:
        config["threads"]
    params:
        cache_dir=config["vep"]["cache_dir"],
        fasta=config["vep"]["fasta"]
    log:
        lambda wildcards: config["outdir"] + "/logs/007_annotation/{sample}_{lane}_annotation.log".format(sample=wildcards.sample, lane=wildcards.lane)
    benchmark:
        lambda wildcards: config["outdir"] + "/benchmarks/007_annotation/{sample}_{lane}_annotation.txt".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        vep \
        --dir_cache {params.cache_dir} \
        --fasta {params.fasta} \
        --input_file {input.vcf} \
        --output_file {output.annotated_vcf} \
        --vcf \
        --fork {threads} \
        --everything \
        > {log} 2>&1
        """

rule extract_on_target_reads:
    message:
        "Extracting on-target reads for sample {wildcards.sample}_{lane}"
    input:
        recal_bam=lambda wildcards: config["outdir"] + "/analysis/006_variant_filtering/{sample}_{lane}/{sample}_{lane}.recal.bam".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        on_target_bam=lambda wildcards: config["outdir"] + "/analysis/007_annotation/{sample}_{lane}/{sample}_{lane}.on_target.bam".format(sample=wildcards.sample, lane=wildcards.lane)
    conda:
        "icc_07_annotation"
    params:
        target=config["target_file"]
    log:
        lambda wildcards: config["outdir"] + "/logs/007_annotation/{sample}_{lane}_extract_on_target_reads.log".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        bedtools intersect -abam {input.recal_bam} -b {params.target} > {output.on_target_bam} 2> {log}
        """

rule flagstat_report:
    message:
        "Generating flagstat report for sample {wildcards.sample}_{lane}"
    input:
        on_target_bam=lambda wildcards: config["outdir"] + "/analysis/007_annotation/{sample}_{lane}/{sample}_{lane}.on_target.bam".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        flagstat=lambda wildcards: config["outdir"] + "/analysis/007_annotation/{sample}_{lane}/{sample}_{lane}.flagstat.txt".format(sample=wildcards.sample, lane=wildcards.lane)
    conda:
        "icc_07_annotation"
    log:
        lambda wildcards: config["outdir"] + "/logs/007_annotation/{sample}_{lane}_flagstat.log".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        samtools flagstat {input.on_target_bam} > {output.flagstat} 2> {log}
        """

rule coverage_analysis:
    message:
        "Performing coverage analysis for sample {wildcards.sample}_{lane}"
    input:
        on_target_bam=lambda wildcards: config["outdir"] + "/analysis/007_annotation/{sample}_{lane}/{sample}_{lane}.on_target.bam".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        coverage_stats=lambda wildcards: config["outdir"] + "/analysis/007_annotation/{sample}_{lane}/{sample}_{lane}.coverage_stats.txt".format(sample=wildcards.sample, lane=wildcards.lane),
        coverage_per_base=lambda wildcards: config["outdir"] + "/analysis/007_annotation/{sample}_{lane}/{sample}_{lane}.coverage_per_base.txt".format(sample=wildcards.sample, lane=wildcards.lane),
        coverage_histogram=lambda wildcards: config["outdir"] + "/analysis/007_annotation/{sample}_{lane}/{sample}_{lane}.coverage_histogram.txt".format(sample=wildcards.sample, lane=wildcards.lane)
    conda:
        "icc_07_annotation"
    params:
        target=config["target_file"]
    log:
        lambda wildcards: config["outdir"] + "/logs/007_annotation/{sample}_{lane}_coverage_analysis.log".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        bedtools coverage -abam {input.on_target_bam} -b {params.target} > {output.coverage_stats} 2> {log}
        bedtools coverage -abam {input.on_target_bam} -b {params.target} -d > {output.coverage_per_base} 2>> {log}
        bedtools coverage -abam {input.on_target_bam} -b {params.target} -hist > {output.coverage_histogram} 2>> {log}
        """