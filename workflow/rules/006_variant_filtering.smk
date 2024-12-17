# A rule to filter variants using GATK VariantFiltration

rule base_recalibrator:
    message:
        "Base recalibration for sample {wildcards.sample}_{lane}"
    input:
        realigned_bam=lambda wildcards: config["outdir"] + "/analysis/005_variant_calling/{sample}_{lane}/{sample}_{lane}.realigned.bam".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        recal_table=lambda wildcards: config["outdir"] + "/analysis/006_variant_filtering/{sample}_{lane}/{sample}_{lane}.recal_data.table".format(sample=wildcards.sample, lane=wildcards.lane)
    conda:
        "gatk"
    params:
        known_sites=config["known_sites"],
        target=config["target_file"]
    log:
        lambda wildcards: config["outdir"] + "/logs/006_variant_filtering/{sample}_{lane}_base_recalibrator.log".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        gatk BaseRecalibrator \
        -R {params.target} \
        -I {input.realigned_bam} \
        --known-sites {params.known_sites} \
        -O {output.recal_table} \
        > {log} 2>&1
        """

rule print_reads:
    message:
        "Printing reads for sample {wildcards.sample}_{lane}"
    input:
        realigned_bam=lambda wildcards: config["outdir"] + "/analysis/005_variant_calling/{sample}_{lane}/{sample}_{lane}.realigned.bam".format(sample=wildcards.sample, lane=wildcards.lane),
        recal_table=lambda wildcards: config["outdir"] + "/analysis/006_variant_filtering/{sample}_{lane}/{sample}_{lane}.recal_data.table".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        recal_bam=lambda wildcards: config["outdir"] + "/analysis/006_variant_filtering/{sample}_{lane}/{sample}_{lane}.recal.bam".format(sample=wildcards.sample, lane=wildcards.lane)
    conda:
        "gatk"
    log:
        lambda wildcards: config["outdir"] + "/logs/006_variant_filtering/{sample}_{lane}_print_reads.log".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        gatk PrintReads \
        -R {params.target} \
        -I {input.realigned_bam} \
        -BQSR {input.recal_table} \
        -O {output.recal_bam} \
        > {log} 2>&1
        """

rule filter_variants:
    message:
        "Filtering variants for sample {wildcards.sample}_{lane}"
    input:
        vcf=lambda wildcards: config["outdir"] + "/analysis/005_variant_calling/{sample}_{lane}/{sample}_{lane}.haplotypecaller.vcf".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        filtered_vcf=lambda wildcards: config["outdir"] + "/analysis/006_variant_filtering/{sample}_{lane}/{sample}_{lane}.filtered.vcf".format(sample=wildcards.sample, lane=wildcards.lane)
    conda:
        "gatk"
    threads:
        config["threads"]
    params:
        ref=config["reference"]
    log:
        lambda wildcards: config["outdir"] + "/logs/006_variant_filtering/{sample}_{lane}_filtering.log".format(sample=wildcards.sample, lane=wildcards.lane)
    benchmark:
        lambda wildcards: config["outdir"] + "/benchmarks/006_variant_filtering/{sample}_{lane}_filtering.txt".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        gatk VariantFiltration \
        -R {params.ref} \
        -V {input.vcf} \
        -O {output.filtered_vcf} \
        --filter-expression "QD < 2.0" \
        --filter-expression "FS > 60.0" \
        --filter-expression "MQ < 40.0" \
        --filter-expression "MQRankSum < -12.5" \
        --filter-expression "ReadPosRankSum < -8.0" \
        --filter-name "QDFilter" \
        --filter-name "FSFilter" \
        --filter-name "MQFilter" \
        --filter-name "MQRankSumFilter" \
        --filter-name "ReadPosFilter" \
        > {log} 2>&1
        """