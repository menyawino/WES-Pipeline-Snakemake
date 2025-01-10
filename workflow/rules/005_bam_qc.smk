rule flagstat_original:
    message:
        "Flagstat for original BAM for sample {wildcards.sample}"
    input:
        bam=rules.apply_bqsr.output.bqsr_bam
    output:
        flagstat_original=config["outdir"] + "/analysis/004_bam_qc/{sample}.bam.flagstat"
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
        {input.bam} \
        2> {log} \
        | sambamba flagstat \
        -t {threads} \
        > {output.flagstat_original} \
        2> {log}
        """

rule flagstat_target:
    message:
        "Flagstat for target BAM for sample {wildcards.sample}"
    input:
        bam_target=rules.filter_bam_target.output.bam_target
    output:
        flagstat_target=config["outdir"] + "/analysis/004_bam_qc/{sample}.target.bam.flagstat"
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
        {input.bam_target} \
        2> {log} \
        | sambamba flagstat \
        -t {threads} \
        > {output.flagstat_target} \
        2> {log}
        """

rule flagstat_prot_coding:
    message:
        "Flagstat for protein coding BAM for sample {wildcards.sample}"
    input:
        bam_prot_coding=rules.filter_bam_prot_coding.output.bam_prot_coding
    output:
        flagstat_prot_coding=config["outdir"] + "/analysis/004_bam_qc/{sample}.prot_coding.bam.flagstat"
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
        {input.bam_prot_coding} \
        2> {log} \
        | sambamba flagstat \
        -t {threads} \
        > {output.flagstat_prot_coding} \
        2> {log}
        """

rule flagstat_canon_tran:
    message:
        "Flagstat for canonical transcript BAM for sample {wildcards.sample}"
    input:
        bam_canon_tran=rules.filter_bam_canon_tran.output.bam_canon_tran
    output:
        flagstat_canon_tran=config["outdir"] + "/analysis/004_bam_qc/{sample}.canon_tran.bam.flagstat"
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
        {input.bam_canon_tran} \
        2> {log} \
        | sambamba flagstat \
        -t {threads} \
        > {output.flagstat_canon_tran} \
        2> {log}
        """

rule coverage_stats:
    message:
        "Calculating coverage stats for protein coding target for sample {wildcards.sample}"
    input:
        bam_prot_coding=rules.filter_bam_prot_coding.output.bam_prot_coding
    output:
        coverage_stats=config["outdir"] + "/analysis/004_bam_qc/{sample}.prot_coding.bam.coverage_stats.txt"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    params:
        CDSFile=config["cds_panel"]
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
        -abam stdin \
        -b {params.CDSFile} \
        2> {log} \
        | sort -k 1,1 -k 2,2n \
        > {output.coverage_stats} \
        2> {log}
        """

rule coverage_stats_per_base:
    message:
        "Calculating coverage stats per base for protein coding target for sample {wildcards.sample}"
    input:
        bam_prot_coding=rules.filter_bam_prot_coding.output.bam_prot_coding
    output:
        coverage_stats_per_base=config["outdir"] + "/analysis/004_bam_qc/{sample}.prot_coding.bam.coverage_per_base.txt"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    params:
        CDSFile=config["cds_panel"]
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
        -b {params.CDSFile} \
        -d \
        | sort -k 1,1 -k 2,2n -k 5,5n \
        > {output.coverage_stats_per_base} \
        2> {log}
        """

rule coverage_hist:
    message:
        "Calculating coverage histogram for protein coding target for sample {wildcards.sample}"
    input:
        bam_prot_coding=rules.filter_bam_prot_coding.output.bam_prot_coding
    output:
        coverage_hist=config["outdir"] + "/analysis/004_bam_qc/{sample}.prot_coding.bam.coverage_hist.txt"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    params:
        CDSFile=config["cds_panel"]
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
        -b {params.CDSFile} \
        -hist \
        > {output.coverage_hist} \
        2> {log}
        """

rule depth_of_coverage:
    message:
        "Calculating depth of coverage for sample {wildcards.sample}"
    input:
        bam_prot_coding=rules.filter_bam_prot_coding.output.bam_prot_coding
    output:
        depth_of_coverage=config["outdir"] + "/analysis/004_bam_qc/{sample}.prot_coding.bam.depth_of_coverage.txt"
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
        --omitDepthOutputAtEachBase false \
        --includeDeletions true \
        --minMappingQuality 10 \
        --minBaseQuality 20 \
        --stop 2000 \
        --includeDeletions \
        --summaryCoverageThreshold 10\
        --outputFormat table \
        --omitDepthOutputAtEachBase \
        &> {log}
        """

rule mean_coverage_per_exon:
    message:
        "Calculating mean coverage per exon for protein coding target for sample {wildcards.sample}"
    input:
        coverage_stats=rules.coverage_stats.output.coverage_stats,
        depth_of_coverage=rules.depth_of_coverage.output.depth_of_coverage
    output:
        mean_coverage=config["outdir"] + "/analysis/004_bam_qc/{sample}.prot_coding.bam.mean_coverage.bed"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    log:
        config["outdir"] + "/logs/004_bam_qc/{sample}_mean_coverage_per_exon.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_bam_qc/{sample}_mean_coverage_per_exon.txt"
    shell:
        """
        awk 'BEGIN {{FS=":"}} {{OFS="\t"}} {{print $1,$2}}' {input.depth_of_coverage} | \
        awk 'BEGIN {{FS="-"}} {{OFS="\t"}} {{print $1,$2}}' | \
        awk 'BEGIN {{FS=" "}} {{OFS="\t"}} {{print $1,$2-1,$3,$4,$5,$9}}' | \
        sed '1d' | \
        sort -k 1,1 -k 2,2n > {output}.temp.depth1.bed

        cut -f4 {input.coverage_stats} | \
        awk 'BEGIN {{FS="|"}} {{print $1}}' > {output}.temp.depth2.bed

        paste <(cut -f1-3 {output}.temp.depth1.bed) \
        <(cut -f4 {output}.temp.depth2.bed) \
        <(cut -f4- {output}.temp.depth1.bed) | \
        cut -f1-4,6 > {output}

        rm -rf {output}.temp*
        """

rule collect_alignment_summary_metrics:
    message:
        "Collecting alignment summary metrics for protein coding target for sample {wildcards.sample}"
    input:
        bam_prot_coding=rules.filter_bam_prot_coding.output.bam_prot_coding
    output:
        alignment_summary_metrics=config["outdir"] + "/analysis/004_bam_qc/{sample}.prot_coding.bam.align_sum_metrics.txt"
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
