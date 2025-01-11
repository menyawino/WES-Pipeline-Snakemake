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
    threads:
        config["threads_low"]
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
    threads:
        config["threads_low"]
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
    threads:
        config["threads_low"]
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
        config["threads_low"]
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