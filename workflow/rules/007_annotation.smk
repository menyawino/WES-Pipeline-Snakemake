# A rule to annotate variants using VEP

rule annotate_variants:
    message:
        "Annotating variants for sample {wildcards.sample}_{lane}"
    input:
        vcf=config["outdir"] + "/analysis/006_variant_filtering/{sample}_{lane}/{sample}_{lane}.filtered.vcf"
    output:
        annotated_vcf=config["outdir"] + "/analysis/007_annotation/{sample}_{lane}/{sample}_{lane}.annotated.vcf"
    conda:
        "icc_07_annotation"
    threads:
        config["threads"]
    params:
        cache_dir=config["vep"]["cache_dir"],
        fasta=config["vep"]["fasta"]
    log:
        config["outdir"] + "/logs/007_annotation/{sample}_{lane}_annotation.log"
    benchmark:
        config["outdir"] + "/benchmarks/007_annotation/{sample}_{lane}_annotation.txt"
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
        recal_bam=config["outdir"] + "/analysis/006_variant_filtering/{sample}_{lane}/{sample}_{lane}.recal.bam"
    output:
        on_target_bam=config["outdir"] + "/analysis/007_annotation/{sample}_{lane}/{sample}_{lane}.on_target.bam"
    conda:
        "icc_07_annotation"
    params:
        target=config["target_file"]
    log:
        config["outdir"] + "/logs/007_annotation/{sample}_{lane}_extract_on_target_reads.log"
    shell:
        """
        bedtools intersect -abam {input.recal_bam} -b {params.target} > {output.on_target_bam} 2> {log}
        """

rule flagstat_report:
    message:
        "Generating flagstat report for sample {wildcards.sample}_{lane}"
    input:
        on_target_bam=config["outdir"] + "/analysis/007_annotation/{sample}_{lane}/{sample}_{lane}.on_target.bam"
    output:
        flagstat=config["outdir"] + "/analysis/007_annotation/{sample}_{lane}/{sample}_{lane}.flagstat.txt"
    conda:
        "icc_07_annotation"
    log:
        config["outdir"] + "/logs/007_annotation/{sample}_{lane}_flagstat.log"
    shell:
        """
        samtools flagstat {input.on_target_bam} > {output.flagstat} 2> {log}
        """

rule coverage_analysis:
    message:
        "Performing coverage analysis for sample {wildcards.sample}_{lane}"
    input:
        on_target_bam=config["outdir"] + "/analysis/007_annotation/{sample}_{lane}/{sample}_{lane}.on_target.bam"
    output:
        coverage_stats=config["outdir"] + "/analysis/007_annotation/{sample}_{lane}/{sample}_{lane}.coverage_stats.txt",
        coverage_per_base=config["outdir"] + "/analysis/007_annotation/{sample}_{lane}/{sample}_{lane}.coverage_per_base.txt",
        coverage_histogram=config["outdir"] + "/analysis/007_annotation/{sample}_{lane}/{sample}_{lane}.coverage_histogram.txt"
    conda:
        "icc_07_annotation"
    params:
        target=config["target_file"]
    log:
        config["outdir"] + "/logs/007_annotation/{sample}_{lane}_coverage_analysis.log"
    shell:
        """
        bedtools coverage -abam {input.on_target_bam} -b {params.target} > {output.coverage_stats} 2> {log}
        bedtools coverage -abam {input.on_target_bam} -b {params.target} -d > {output.coverage_per_base} 2>> {log}
        bedtools coverage -abam {input.on_target_bam} -b {params.target} -hist > {output.coverage_histogram} 2>> {log}
        """