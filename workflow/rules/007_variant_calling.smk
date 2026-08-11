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
    resources:
        mem_mb=config.get("mem_high", 32768)
    params:
        ref=config["reference_genome"],
        target=config["icc_panel"]
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_haplotypecaller.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_haplotypecaller.txt"
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_mb}m" HaplotypeCallerSpark \
        -R {params.ref} \
        -I {input.bam} \
        -O {output.gvcf} \
        -ERC GVCF \
        -OVI true \
        --intervals {params.target} \
        --spark-master local[{threads}] \
        &> {log}
        """
    
rule genotype_gvcfs:
    message:
        "Genotyping GVCFs for sample {wildcards.sample}"
    input:
        gvcf=rules.haplotypecaller.output.gvcf
    output:
        vcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.genotyped.vcf"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    resources:
        mem_mb=config.get("mem_mid", 16384)
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
        gatk --java-options "-Xmx{resources.mem_mb}m" GenotypeGVCFs \
        -R {params.ref} \
        -V {input.gvcf} \
        -O {output.vcf} \
        -G StandardAnnotation \
        -A DepthPerAlleleBySample \
        -A Coverage \
        -A InbreedingCoeff \
        -A QualByDepth \
        -A FS \
        -A SOR \
        -A ReadPosRankSum \
        -A MQRankSum \
        --intervals {params.target} \
        --dbsnp {params.dbsnp} \
        &> {log}
        """

rule split_vcfs:
    message:
        "Splitting VCFs into SNPs and Indels for sample {wildcards.sample}"
    input:
        vcf=rules.genotype_gvcfs.output.vcf
    output:
        snp_vcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.genotyped.snp.vcf",
        indel_vcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.genotyped.indel.vcf"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
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

