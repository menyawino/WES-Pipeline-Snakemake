# A rule to filter variants using GATK VariantFiltration

rule base_recalibrator:
    message:
        "Base recalibration for sample {wildcards.sample}"
    input:
        realigned_bam=config["outdir"] + "/analysis/005_variant_calling/{sample}.realigned.bam"
    output:
        recal_table=config["outdir"] + "/analysis/006_variant_filtering/{sample}.recal_data.table"
    conda:
        "gatk"
    params:
        known_sites=config["gatk"]["BaseRecalibrator"]["known_sites"],
        target=config["target_file"]
    log:
        config["outdir"] + "/logs/006_variant_filtering/{sample}_base_recalibrator.log"
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
        "Printing reads for sample {wildcards.sample}"
    input:
        realigned_bam=config["outdir"] + "/analysis/005_variant_calling/{sample}.realigned.bam",
        recal_table=config["outdir"] + "/analysis/006_variant_filtering/{sample}.recal_data.table"
    output:
        recal_bam=config["outdir"] + "/analysis/006_variant_filtering/{sample}.recal.bam"
    conda:
        "gatk"
    log:
        config["outdir"] + "/logs/006_variant_filtering/{sample}_print_reads.log"
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
        "Filtering variants for sample {wildcards.sample}"
    input:
        vcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.haplotypecaller.vcf"
    output:
        filtered_vcf=config["outdir"] + "/analysis/006_variant_filtering/{sample}.filtered.vcf"
    conda:
        "gatk"
    threads:
        config["threads"]
    params:
        ref=config["reference"]
    log:
        config["outdir"] + "/logs/006_variant_filtering/{sample}_filtering.log"
    benchmark:
        config["outdir"] + "/benchmarks/006_variant_filtering/{sample}_filtering.txt"
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