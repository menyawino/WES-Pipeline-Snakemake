# A rule to call variants using GATK HaplotypeCaller and GenotypeGVCFs

rule haplotypecaller:
    message:
        "Calling variants with HaplotypeCaller for sample {wildcards.sample}"
    input:
        bam=rules.filter_bam_target.output.bam_target
    output:
        gvcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.haplotypecaller.g.vcf"
    conda:
        "icc_gatk"
    threads:
        config["threads_high"]
    params:
        ref=config["reference_genome"],
        target=config["icc_panel"]
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_haplotypecaller.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_haplotypecaller.txt"
    shell:
        """
        gatk HaplotypeCallerSpark \
        -R {params.ref} \
        -I {input.bam} \
        -O {output.gvcf} \
        -ERC GVCF \
        -OVI true \
        --intervals {params.target} \
        --spark-master local[{threads}] \
        &> {log}
        """
    
rule index_gvcf:
    message:
        "Indexing GVCF for sample {wildcards.sample}"
    input:
        gvcf=rules.haplotypecaller.output.gvcf
    output:
        gvcf_index=config["outdir"] + "/analysis/005_variant_calling/{sample}.haplotypecaller.g.vcf.idx"
    conda:
        "icc_gatk"
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_index_gvcf.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_index_gvcf.txt"
    shell:
        """
        gatk IndexFeatureFile \
        -I {input.gvcf} \
        &> {log}
        """

rule variant_annotation:
    message:
        "Annotating variants for sample {wildcards.sample}"
    input:
        vcf=rules.haplotypecaller.output.gvcf,
        idx=rules.index_gvcf.output.gvcf_index
    output:
        annotated_vcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.annotated.vcf"
    conda:
        "icc_gatk"
    params:
        ref=config["reference_genome"],
        target=config["icc_panel"],
        dbsnp=config["dbsnp"]
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_variant_annotation.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_variant_annotation.txt"
    shell:
        """
        gatk VariantAnnotator \
        -R {params.ref} \
        -V {input.vcf} \
        -O {output.annotated_vcf} \
        -A Coverage \
        -A DepthPerAlleleBySample \
        -A QualByDepth \
        -A InbreedingCoeff \
        -G StandardAnnotation \
        --dbsnp {params.dbsnp} \
        --intervals {params.target} \
        &> {log}
        """

rule genotype_gvcfs:
    message:
        "Genotyping GVCFs for sample {wildcards.sample}"
    input:
        gvcf=rules.variant_annotation.output.annotated_vcf
    output:
        vcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.genotyped.vcf"
    conda:
        "icc_gatk"
    params:
        ref=config["reference_genome"],
        target=config["icc_panel"],
        dbsnp=config["dbsnp"]
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_genotype_gvcfs.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_genotype_gvcfs.txt"
    shell:
        """
        gatk GenotypeGVCFs \
        -R {params.ref} \
        -V {input.gvcf} \
        -O {output.vcf} \
        -A DepthPerAlleleBySample \
        -A Coverage \
        -A InbreedingCoeff \
        --intervals {params.target} \
        --dbsnp {params.dbsnp} \
        &> {log}
        """

rule vqsr_snp_recalibration:
    message:
        "Recalibrating SNPs for sample {wildcards.sample}"
    input:
        vcf=rules.genotype_gvcfs.output.vcf
    output:
        recal_snp=config["outdir"] + "/analysis/005_variant_calling/{sample}.recalibrate_SNP.recal",
        tranches_snp=config["outdir"] + "/analysis/005_variant_calling/{sample}.recalibrate_SNP.tranches"
    conda:
        "icc_gatk"
    params:
        ref=config["reference_genome"],
        hapmap=config["hapmap"],
        omni=config["omni"],
        tenk=config["tenk"],
        dbsnp=config["dbsnp"]
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_vqsr_snp_recalibration.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_vqsr_snp_recalibration.txt"
    shell:
        """
        gatk VariantRecalibrator \
        -R {params.ref} \
        -V {input.vcf} \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} \
        --resource:omni,known=false,training=true,truth=true,prior=12.0 {params.omni} \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 {params.tenk} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \
        -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
        -mode SNP \
        -O {output.recal_snp} \
        --tranches-file {output.tranches_snp} \
        &> {log}
        """

rule apply_vqsr_snp:
    message:
        "Applying VQSR for SNPs for sample {wildcards.sample}"
    input:
        vcf=rules.genotype_gvcfs.output.vcf,
        recal_snp=rules.vqsr_snp_recalibration.output.recal_snp,
        tranches_snp=rules.vqsr_snp_recalibration.output.tranches_snp
    output:
        recal_vcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.recalibrated.vcf"
    conda:
        "icc_gatk"
    params:
        ref=config["reference_genome"]
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_apply_vqsr_snp.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_apply_vqsr_snp.txt"
    shell:
        """
        gatk ApplyVQSR \
        -R {params.ref} \
        -V {input.vcf} \
        -O {output.recal_vcf} \
        --recal-file {input.recal_snp} \
        --tranches-file {input.tranches_snp} \
        -mode SNP \
        --truth-sensitivity-filter-level 99.0 \
        &> {log}
        """

rule vqsr_indel_recalibration:
    message:
        "Recalibrating Indels for sample {wildcards.sample}"
    input:
        vcf=rules.apply_vqsr_snp.output.recal_vcf
    output:
        recal_indel=config["outdir"] + "/analysis/005_variant_calling/{sample}.recalibrate_INDEL.recal",
        tranches_indel=config["outdir"] + "/analysis/005_variant_calling/{sample}.recalibrate_INDEL.tranches"
    conda:
        "icc_gatk"
    params:
        ref=config["reference_genome"],
        mills=config["mills"],
        dbsnp=config["dbsnp"]
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_vqsr_indel_recalibration.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_vqsr_indel_recalibration.txt"
    shell:
        """
        gatk VariantRecalibrator \
        -R {params.ref} \
        -V {input.vcf} \
        --resource:mills,known=true,training=true,truth=true,prior=12.0 {params.mills} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \
        -an QD -an ReadPosRankSum -an FS -an SOR -an DP \
        -mode INDEL \
        -O {output.recal_indel} \
        --tranches-file {output.tranches_indel} \
        &> {log}
        """

rule apply_vqsr_indel:
    message:
        "Applying VQSR for Indels for sample {wildcards.sample}"
    input:
        vcf=rules.apply_vqsr_snp.output.recal_vcf,
        recal_indel=rules.vqsr_indel_recalibration.output.recal_indel,
        tranches_indel=rules.vqsr_indel_recalibration.output.tranches_indel
    output:
        recal_vcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.recalibrated.indel.vcf"  # Update output filename to avoid conflict
    conda:
        "icc_gatk"
    params:
        ref=config["reference_genome"]
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_apply_vqsr_indel.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_apply_vqsr_indel.txt"
    shell:
        """
        gatk ApplyVQSR \
        -R {params.ref} \
        -V {input.vcf} \
        -O {output.recal_vcf} \
        --recal-file {input.recal_indel} \
        --tranches-file {input.tranches_indel} \
        -mode INDEL \
        --truth-sensitivity-filter-level 99.0 \
        &> {log}
        """

rule split_vcfs:
    message:
        "Splitting VCFs into SNPs and Indels for sample {wildcards.sample}"
    input:
        vcf=rules.apply_vqsr_indel.output.recal_vcf
    output:
        snp_vcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.genotyped.snp.vcf",
        indel_vcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.genotyped.indel.vcf"
    conda:
        "icc_gatk"
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_split_vcfs.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_split_vcfs.txt"
    shell:
        """
        gatk SplitVcfs \
        I={input.vcf} \
        SNP_OUTPUT={output.snp_vcf} \
        INDEL_OUTPUT={output.indel_vcf} \
        STRICT=false \
        &> {log}
        """