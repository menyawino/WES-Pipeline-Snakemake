# Variant Filtering (Unified for GATK and DeepVariant)

rule filter_snps:
    message:
        "Filtering {wildcards.caller} SNPs for sample {wildcards.sample}"
    input:
        snp_vcf=rules.split_vcfs.output.snp_vcf,
        snp_tbi=rules.split_vcfs.output.snp_tbi
    output:
        filtered_snp_vcf=config["outdir"] + "/analysis/006_variant_filtering/{sample}.{caller}.filtered.snp.vcf"
    container:
        "docker://broadinstitute/gatk:4.4.0.0"
    threads:
        config.get("threads_mid", 8)
    resources:
        mem_mb=config.get("mem_mid", 16384),
        tmpdir=config.get("tmpdir", "/tmp")
    params:
        ref="/dev/shm/wes_ref_grch38/GRCh38.primary_assembly.genome.fa",
        target=config["icc_panel"]
    log:
        config["outdir"] + "/logs/006_variant_filtering/{sample}_{caller}_filter_snps.log"
    benchmark:
        config["outdir"] + "/benchmarks/006_variant_filtering/{sample}_{caller}_filter_snps.txt"
    shell:
        """
        mkdir -p "{resources.tmpdir}"
        if [ "{wildcards.caller}" = "gatk" ]; then
            gatk --java-options "-XX:+UseParallelGC -Xmx{resources.mem_mb}m" VariantFiltration \
            -R "{params.ref}" \
            -V "{input.snp_vcf}" \
            -O "{output.filtered_snp_vcf}" \
            --filter-name "QDFilter" --filter-expression "QD < 2.0" \
            --filter-name "FSFilter" --filter-expression "FS > 60.0" \
            --filter-name "MQFilter" --filter-expression "MQ < 40.0" \
            --filter-name "MQRankSumFilter" --filter-expression "MQRankSum < -12.5" \
            --filter-name "ReadPosFilter" --filter-expression "ReadPosRankSum < -8.0" \
            --intervals "{params.target}" \
            --create-output-variant-index false \
            --tmp-dir "{resources.tmpdir}" \
            &> "{log}"
        else
            # DeepVariant standard quality filtering (tag non-PASS variants)
            bcftools filter -e 'FILTER != "PASS" && QUAL < 15.0' -s "LowQual" -m + "{input.snp_vcf}" -O v -o "{output.filtered_snp_vcf}" > "{log}" 2>&1
        fi
        """

rule filter_indels:
    message:
        "Filtering {wildcards.caller} Indels for sample {wildcards.sample}"
    input:
        indel_vcf=rules.split_vcfs.output.indel_vcf,
        indel_tbi=rules.split_vcfs.output.indel_tbi
    output:
        filtered_indel_vcf=config["outdir"] + "/analysis/006_variant_filtering/{sample}.{caller}.filtered.indel.vcf"
    container:
        "docker://broadinstitute/gatk:4.4.0.0"
    threads:
        config.get("threads_mid", 8)
    resources:
        mem_mb=config.get("mem_mid", 16384),
        tmpdir=config.get("tmpdir", "/tmp")
    params:
        ref="/dev/shm/wes_ref_grch38/GRCh38.primary_assembly.genome.fa",
        target=config["icc_panel"]
    log:
        config["outdir"] + "/logs/006_variant_filtering/{sample}_{caller}_filter_indels.log"
    benchmark:
        config["outdir"] + "/benchmarks/006_variant_filtering/{sample}_{caller}_filter_indels.txt"
    shell:
        """
        mkdir -p "{resources.tmpdir}"
        if [ "{wildcards.caller}" = "gatk" ]; then
            gatk --java-options "-XX:+UseParallelGC -Xmx{resources.mem_mb}m" VariantFiltration \
            -R "{params.ref}" \
            -V "{input.indel_vcf}" \
            -O "{output.filtered_indel_vcf}" \
            --filter-name "QDFilter" --filter-expression "QD < 2.0" \
            --filter-name "ReadPosFilter" --filter-expression "ReadPosRankSum < -20.0" \
            --filter-name "FSFilter" --filter-expression "FS > 200.0" \
            --intervals "{params.target}" \
            --create-output-variant-index false \
            --tmp-dir "{resources.tmpdir}" \
            &> "{log}"
        else
            # DeepVariant standard quality filtering
            bcftools filter -e 'FILTER != "PASS" && QUAL < 15.0' -s "LowQual" -m + "{input.indel_vcf}" -O v -o "{output.filtered_indel_vcf}" > "{log}" 2>&1
        fi
        """