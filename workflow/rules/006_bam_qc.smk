rule flagstat_original:
    message:
        "Flagstat for original BAM for sample {wildcards.sample}"
    input:
        bam=rules.apply_bqsr.output.bqsr_bam
    output:
        flagstat_original=config["outdir"] + "/analysis/004_bam_qc/{sample}.flagstat"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_flagstat_original.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_flagstat_original.txt"
    shell:
        """
        samtools flagstat \
        -@ {threads} \
        "{input.bam}" \
        > "{output.flagstat_original}" \
        2> "{log}"
        """

rule flagstat_target:
    message:
        "Flagstat for target BAM for sample {wildcards.sample}"
    input:
        bam_target=rules.filter_bam_target.output.bam_target
    output:
        flagstat_target=config["outdir"] + "/analysis/004_bam_qc/{sample}.target.flagstat"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_flagstat_target.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_flagstat_target.txt"
    shell:
        """
        samtools flagstat \
        -@ {threads} \
        "{input.bam_target}" \
        > "{output.flagstat_target}" \
        2> "{log}"
        """

rule flagstat_prot_coding:
    message:
        "Flagstat for protein coding BAM for sample {wildcards.sample}"
    input:
        bam_prot_coding=rules.filter_bam_prot_coding.output.bam_prot_coding
    output:
        flagstat_prot_coding=config["outdir"] + "/analysis/004_bam_qc/{sample}.prot_coding.flagstat"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_flagstat_prot_coding.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_flagstat_prot_coding.txt"
    shell:
        """
        samtools flagstat \
        -@ {threads} \
        "{input.bam_prot_coding}" \
        > "{output.flagstat_prot_coding}" \
        2> "{log}"
        """

rule flagstat_canon_tran:
    message:
        "Flagstat for canonical transcript BAM for sample {wildcards.sample}"
    input:
        bam_canon_tran=rules.filter_bam_canon_tran.output.bam_canon_tran
    output:
        flagstat_canon_tran=config["outdir"] + "/analysis/004_bam_qc/{sample}.canon_tran.flagstat"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_flagstat_canon_tran.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_flagstat_canon_tran.txt"
    shell:
        """
        samtools flagstat \
        -@ {threads} \
        "{input.bam_canon_tran}" \
        > "{output.flagstat_canon_tran}" \
        2> "{log}"
        """

rule coverage_stats:
    message:
        "Calculating coverage stats for protein coding target for sample {wildcards.sample}"
    input:
        bam_prot_coding=rules.filter_bam_prot_coding.output.bam_prot_coding
    output:
        coverage_stats=config["outdir"] + "/analysis/004_bam_qc/{sample}.prot_coding.coverage_stats.txt"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    params:
        cds_file=config["cds_panel"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_coverage_stats.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_coverage_stats.txt"
    shell:
        """
        bedtools coverage \
        -abam "{input.bam_prot_coding}" \
        -b "{params.cds_file}" \
        2> "{log}" \
        | sort -k 1,1 -k 2,2n \
        > "{output.coverage_stats}" \
        2>> "{log}"
        """

rule coverage_stats_target:
    message:
        "Calculating coverage stats for target BAM for sample {wildcards.sample}"
    input:
        bam_target=rules.filter_bam_target.output.bam_target
    output:
        coverage_stats_target=config["outdir"] + "/analysis/004_bam_qc/{sample}.target.coverage_stats.txt"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    params:
        cds_file=config["icc_panel"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_coverage_stats_target.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_coverage_stats_target.txt"
    shell:
        """
        bedtools coverage \
        -abam "{input.bam_target}" \
        -b "{params.cds_file}" \
        2> "{log}" \
        | sort -k 1,1 -k 2,2n \
        > "{output.coverage_stats_target}" \
        2>> "{log}"
        """

rule coverage_stats_per_base:
    message:
        "Calculating coverage stats per base for protein coding target for sample {wildcards.sample}"
    input:
        bam_prot_coding=rules.filter_bam_prot_coding.output.bam_prot_coding
    output:
        coverage_stats_per_base=config["outdir"] + "/analysis/004_bam_qc/{sample}.prot_coding.per-base.bed.gz"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    params:
        cds_file=config["cds_panel"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_coverage_stats_per_base.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_coverage_stats_per_base.txt"
    shell:
        """
        bedtools coverage \
        -abam "{input.bam_prot_coding}" \
        -b "{params.cds_file}" \
        -d 2> "{log}" \
        | gzip -c \
        > "{output.coverage_stats_per_base}" \
        2>> "{log}"
        """
        
rule coverage_stats_per_base_target:
    message:
        "Calculating coverage stats per base for target BAM for sample {wildcards.sample}"
    input:
        bam_target=rules.filter_bam_target.output.bam_target
    output:
        coverage_stats_per_base_target=config["outdir"] + "/analysis/004_bam_qc/{sample}.target.coverage_per_base.txt"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    params:
        cds_file=config["icc_panel"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_coverage_stats_per_base_target.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_coverage_stats_per_base_target.txt"
    shell:
        """
        bedtools coverage \
        -abam "{input.bam_target}" \
        -b "{params.cds_file}" \
        -d \
        > "{output.coverage_stats_per_base_target}" \
        2> "{log}"
        """

rule coverage_hist:
    message:
        "Calculating coverage histogram for protein coding target for sample {wildcards.sample}"
    input:
        bam_prot_coding=rules.filter_bam_prot_coding.output.bam_prot_coding
    output:
        coverage_hist=config["outdir"] + "/analysis/004_bam_qc/{sample}.prot_coding.coverage_hist.txt"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    params:
        cds_file=config["cds_panel"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_coverage_hist.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_coverage_hist.txt"
    shell:
        """
        bedtools coverage \
        -abam "{input.bam_prot_coding}" \
        -b "{params.cds_file}" \
        -hist \
        > "{output.coverage_hist}" \
        2> "{log}"
        """

rule coverage_hist_target:
    message:
        "Calculating coverage histogram for target BAM for sample {wildcards.sample}"
    input:
        bam_target=rules.filter_bam_target.output.bam_target
    output:
        coverage_hist_target=config["outdir"] + "/analysis/004_bam_qc/{sample}.target.coverage_hist.txt"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    params:
        cds_file=config["icc_panel"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_coverage_hist_target.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_coverage_hist_target.txt"
    shell:
        """
        bedtools coverage \
        -abam "{input.bam_target}" \
        -b "{params.cds_file}" \
        -hist \
        > "{output.coverage_hist_target}" \
        2> "{log}"
        """

rule fast_bam_qc_prot_coding:
    message:
        "Fast BAM QC (Depth and Picard) for protein coding target for sample {wildcards.sample}"
    input:
        bam=rules.filter_bam_prot_coding.output.bam_prot_coding,
        bai=rules.filter_bam_prot_coding.output.bai_prot_coding
    output:
        depth_of_coverage=config["outdir"] + "/analysis/004_bam_qc/{sample}.prot_coding.depth_of_coverage.sample_summary",
        alignment_summary_metrics=config["outdir"] + "/analysis/004_bam_qc/{sample}.prot_coding.align_sum_metrics.txt"
    conda:
        "005_omni_bam"
    threads:
        config["threads_mid"]
    params:
        cds_file=config["cds_panel"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_fast_qc_prot_coding.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_fast_qc_prot_coding.txt"
    shell:
        """
        python workflow/scripts/fast_bam_qc.py \
        --bam "{input.bam}" \
        --bed "{params.cds_file}" \
        --out-depth "{output.depth_of_coverage}" \
        --out-metrics "{output.alignment_summary_metrics}" \
        --threads {threads} \
        > "{log}" 2>&1
        """

rule fast_bam_qc_target:
    message:
        "Fast BAM QC (Depth and Picard) for target BAM for sample {wildcards.sample}"
    input:
        bam=rules.filter_bam_target.output.bam_target,
        bai=rules.filter_bam_target.output.bai_target
    output:
        depth_of_coverage_target=config["outdir"] + "/analysis/004_bam_qc/{sample}.target.depth_of_coverage.sample_summary",
        alignment_summary_metrics_target=config["outdir"] + "/analysis/004_bam_qc/{sample}.target.align_sum_metrics.txt"
    conda:
        "005_omni_bam"
    threads:
        config["threads_mid"]
    params:
        cds_file=config["icc_panel"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_fast_qc_target.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_fast_qc_target.txt"
    shell:
        """
        python workflow/scripts/fast_bam_qc.py \
        --bam "{input.bam}" \
        --bed "{params.cds_file}" \
        --out-depth "{output.depth_of_coverage_target}" \
        --out-metrics "{output.alignment_summary_metrics_target}" \
        --threads {threads} \
        > "{log}" 2>&1
        """

rule mean_coverage_per_exon:
    message:
        "Calculating mean coverage per exon for protein coding target for sample {wildcards.sample}"
    input:
        bam_prot_coding=rules.filter_bam_prot_coding.output.bam_prot_coding
    output:
        mean_coverage=config["outdir"] + "/analysis/004_bam_qc/{sample}.prot_coding.mean_coverage.bed"
    conda:
        "icc_gatk"
    params:
        cds_file=config["cds_panel"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_mean_coverage_per_exon.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_mean_coverage_per_exon.txt"
    shell:
        """
        bedtools coverage \
        -a "{params.cds_file}" \
        -b "{input.bam_prot_coding}" \
        -mean \
        > "{output.mean_coverage}" \
        2> "{log}"
        """

rule mean_coverage_per_exon_target:
    message:
        "Calculating mean coverage per exon for target BAM for sample {wildcards.sample}"
    input:
        bam_target=rules.filter_bam_target.output.bam_target
    output:
        mean_coverage_target=config["outdir"] + "/analysis/004_bam_qc/{sample}.target.mean_coverage.bed"
    conda:
        "icc_gatk"
    params:
        cds_file=config["icc_panel"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_mean_coverage_per_exon_target.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_mean_coverage_per_exon_target.txt"
    shell:
        """
        bedtools coverage \
        -a "{params.cds_file}" \
        -b "{input.bam_target}" \
        -mean \
        > "{output.mean_coverage_target}" \
        2> "{log}"
        """



rule qc_report:
    input:
        flagstat_original = rules.flagstat_original.output.flagstat_original,
        flagstat_target = rules.flagstat_target.output.flagstat_target,
        coverage_stats = rules.coverage_stats.output.coverage_stats,
        coverage_stats_target = rules.coverage_stats_target.output.coverage_stats_target,
        coverage_hist = rules.coverage_hist.output.coverage_hist,
        coverage_hist_target = rules.coverage_hist_target.output.coverage_hist_target,
        depth_of_coverage = rules.fast_bam_qc_prot_coding.output.depth_of_coverage,
        depth_of_coverage_target = rules.fast_bam_qc_target.output.depth_of_coverage_target,
        mean_coverage = rules.mean_coverage_per_exon.output.mean_coverage,
        mean_coverage_target = rules.mean_coverage_per_exon_target.output.mean_coverage_target,
        alignment_summary_metrics = rules.fast_bam_qc_prot_coding.output.alignment_summary_metrics,
        alignment_summary_metrics_target = rules.fast_bam_qc_target.output.alignment_summary_metrics_target
    output:
        qc_metrics = config["outdir"] + "/analysis/004_bam_qc/{sample}.qc_metrics.tsv"
    conda:
        "icc_gatk"
    script:
        "../scripts/qc_report.py"
