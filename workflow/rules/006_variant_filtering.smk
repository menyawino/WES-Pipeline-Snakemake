# Variant Filtering

rule filter_snps:
    message:
        "Filtering SNPs for sample {wildcards.sample}"
    input:
        snp_vcf=rules.split_vcfs.output.snp_vcf
    output:
        filtered_snp_vcf=config["outdir"] + "/analysis/006_variant_filtering/{sample}.filtered.snp.vcf"
    conda:
        "icc_gatk"
    params:
        ref=config["reference_genome"],
        target=config["icc_panel"]
    log:
        config["outdir"] + "/logs/006_variant_filtering/{sample}_filter_snps.log"
    benchmark:
        config["outdir"] + "/benchmarks/006_variant_filtering/{sample}_filter_snps.txt"
    shell:
        """
        gatk VariantFiltration \
        -R {params.ref} \
        -V {input.snp_vcf} \
        -O {output.filtered_snp_vcf} \
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
        --intervals {params.target} \
        &> {log}
        """

rule filter_indels:
    message:
        "Filtering Indels for sample {wildcards.sample}"
    input:
        indel_vcf=rules.split_vcfs.output.indel_vcf
    output:
        filtered_indel_vcf=config["outdir"] + "/analysis/006_variant_filtering/{sample}.filtered.indel.vcf"
    conda:
        "icc_gatk"
    params:
        ref=config["reference_genome"],
        target=config["icc_panel"]
    log:
        config["outdir"] + "/logs/006_variant_filtering/{sample}_filter_indels.log"
    benchmark:
        config["outdir"] + "/benchmarks/006_variant_filtering/{sample}_filter_indels.txt"
    shell:
        """
        gatk VariantFiltration \
        -R {params.ref} \
        -V {input.indel_vcf} \
        -O {output.filtered_indel_vcf} \
        --filter-expression "QD < 2.0" \
        --filter-expression "ReadPosRankSum < -20.0" \
        --filter-expression "FS > 200.0" \
        --filter-name "QDFilter" \
        --filter-name "ReadPosFilter" \
        --filter-name "FSFilter" \
        --intervals {params.target} \
        &> {log}
        """