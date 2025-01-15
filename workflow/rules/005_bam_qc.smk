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
        sambamba view \
        -t {threads} \
        -F "not secondary_alignment" \
        -f bam \
        -l 0 \
        {input.bam} \
        2> {log} \
        | sambamba flagstat \
        -t {threads} \
        /dev/stdin \
        > {output.flagstat_original} \
        2> {log}
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
        sambamba view \
        -t {threads} \
        -F "not secondary_alignment" \
        -f bam \
        -l 0 \
        {input.bam_target} \
        2> {log} \
        | sambamba flagstat \
        -t {threads} \
        /dev/stdin \
        > {output.flagstat_target} \
        2> {log}
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
        sambamba view \
        -t {threads} \
        -F "not secondary_alignment" \
        -f bam \
        -l 0 \
        {input.bam_prot_coding} \
        2> {log} \
        | sambamba flagstat \
        -t {threads} \
        /dev/stdin \
        > {output.flagstat_prot_coding} \
        2> {log}
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
        sambamba view \
        -t {threads} \
        -F "not secondary_alignment" \
        -f bam \
        -l 0 \
        {input.bam_canon_tran} \
        2> {log} \
        | sambamba flagstat \
        -t {threads} \
        /dev/stdin \
        > {output.flagstat_canon_tran} \
        2> {log}
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
        sambamba view \
        -t {threads} \
        -F "not duplicate" \
        -f bam \
        {input.bam_prot_coding} \
        2> {log} \
        | bedtools coverage \
        -abam /dev/stdin \
        -b {params.cds_file} \
        2> {log} \
        | sort -k 1,1 -k 2,2n \
        > {output.coverage_stats} \
        2> {log}
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
        cds_file=config["cds_panel"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_coverage_stats_target.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_coverage_stats_target.txt"
    shell:
        """
        sambamba view \
        -t {threads} \
        -F "not duplicate" \
        -f bam \
        {input.bam_target} \
        2> {log} \
        | bedtools coverage \
        -abam /dev/stdin \
        -b {params.cds_file} \
        2> {log} \
        | sort -k 1,1 -k 2,2n \
        > {output.coverage_stats_target} \
        2> {log}
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
        sambamba view \
        -t {threads} \
        -F "not duplicate" \
        -f bam \
        {input.bam_prot_coding} \
        2> {log} \
        | bedtools coverage \
        -abam stdin \
        -b {params.cds_file} \
        -d \
        > {output.coverage_stats_per_base} \
        2> {log}
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
        cds_file=config["cds_panel"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_coverage_stats_per_base_target.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_coverage_stats_per_base_target.txt"
    shell:
        """
        sambamba view \
        -t {threads} \
        -F "not duplicate" \
        -f bam \
        {input.bam_target} \
        2> {log} \
        | bedtools coverage \
        -abam stdin \
        -b {params.cds_file} \
        -d \
        > {output.coverage_stats_per_base_target} \
        2> {log}
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
        sambamba view \
        -t {threads} \
        -F "not duplicate" \
        -f bam \
        {input.bam_prot_coding} \
        2> {log} \
        | bedtools coverage \
        -abam stdin \
        -b {params.cds_file} \
        -hist \
        > {output.coverage_hist} \
        2> {log}
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
        cds_file=config["cds_panel"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_coverage_hist_target.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_coverage_hist_target.txt"
    shell:
        """
        sambamba view \
        -t {threads} \
        -F "not duplicate" \
        -f bam \
        {input.bam_target} \
        2> {log} \
        | bedtools coverage \
        -abam stdin \
        -b {params.cds_file} \
        -hist \
        > {output.coverage_hist_target} \
        2> {log}
        """

rule depth_of_coverage:
    message:
        "Calculating depth of coverage for sample {wildcards.sample}"
    input:
        bam_prot_coding=rules.filter_bam_prot_coding.output.bam_prot_coding
    output:
        depth_of_coverage=config["outdir"] + "/analysis/004_bam_qc/{sample}.prot_coding.depth_of_coverage"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    params:
        ref=config["reference_genome"],
        cds_file=config["cds_panel"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_depth_of_coverage.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_depth_of_coverage.txt"
    shell:
        """
        gatk DepthOfCoverage \
        -R {params.ref} \
        -I {input.bam_prot_coding} \
        -O {output.depth_of_coverage} \
        -L {params.cds_file} \
        &> {log}
        """

rule depth_of_coverage_target:
    message:
        "Calculating depth of coverage for target BAM for sample {wildcards.sample}"
    input:
        bam_target=rules.filter_bam_target.output.bam_target
    output:
        depth_of_coverage_target=config["outdir"] + "/analysis/004_bam_qc/{sample}.target.depth_of_coverage"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    params:
        ref=config["reference_genome"],
        cds_file=config["cds_panel"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_depth_of_coverage_target.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_depth_of_coverage_target.txt"
    shell:
        """
        gatk DepthOfCoverage \
        -R {params.ref} \
        -I {input.bam_target} \
        -O {output.depth_of_coverage_target} \
        -L {params.cds_file} \
        &> {log}
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
    threads:
        config["threads_mid"]
    params:
        cds_file=config["cds_panel"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_mean_coverage_per_exon.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_mean_coverage_per_exon.txt"
    shell:
        """
        bedtools coverage \
        -a {params.cds_file} \
        -b {input.bam_prot_coding} \
        -mean \
        > {output.mean_coverage}
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
    threads:
        config["threads_mid"]
    params:
        cds_file=config["cds_panel"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_mean_coverage_per_exon_target.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_mean_coverage_per_exon_target.txt"
    shell:
        """
        bedtools coverage \
        -a {params.cds_file} \
        -b {input.bam_target} \
        -mean \
        > {output.mean_coverage_target}
        """

rule collect_alignment_summary_metrics:
    message:
        "Collecting alignment summary metrics for protein coding target for sample {wildcards.sample}"
    input:
        bam_prot_coding=rules.filter_bam_prot_coding.output.bam_prot_coding
    output:
        alignment_summary_metrics=config["outdir"] + "/analysis/004_bam_qc/{sample}.prot_coding.align_sum_metrics.txt"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    params:
        ref=config["reference_genome"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_collect_alignment_summary_metrics.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_collect_alignment_summary_metrics.txt"
    shell:
        """
        gatk CollectAlignmentSummaryMetrics \
        INPUT={input.bam_prot_coding} \
        OUTPUT={output.alignment_summary_metrics} \
        REFERENCE_SEQUENCE={params.ref} \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT \
        &> {log}
        """

rule collect_alignment_summary_metrics_target:
    message:
        "Collecting alignment summary metrics for target BAM for sample {wildcards.sample}"
    input:
        bam_target=rules.filter_bam_target.output.bam_target
    output:
        alignment_summary_metrics_target=config["outdir"] + "/analysis/004_bam_qc/{sample}.target.align_sum_metrics.txt"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    params:
        ref=config["reference_genome"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_collect_alignment_summary_metrics_target.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_collect_alignment_summary_metrics_target.txt"
    shell:
        """
        gatk CollectAlignmentSummaryMetrics \
        INPUT={input.bam_target} \
        OUTPUT={output.alignment_summary_metrics_target} \
        REFERENCE_SEQUENCE={params.ref} \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT \
        &> {log}
        """