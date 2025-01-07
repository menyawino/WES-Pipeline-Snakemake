# A rule to call variants using GATK HaplotypeCaller and GenotypeGVCFs

rule haplotypecaller:
    message:
        "Calling variants with HaplotypeCaller for sample {wildcards.sample}"
    input:
        bam=rules.apply_bqsr.output.bqsr_bam
    output:
        gvcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.haplotypecaller.g.vcf",
        bam=config["outdir"] + "/analysis/005_variant_calling/{sample}.haplotypecaller.bam"
    conda:
        "icc_gatk"
    threads:
        config["threads"]
    params:
        ref=config["reference_genome"],
        dbsnp=config["dbsnp"],
        dcovg=config["gatk"]["HaplotypeCaller"]["dcovg"],
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
        -A Coverage \
	    -A AlleleBalance \
	    -G Standard \
        -ERC GVCF \
        --downsample_to_coverage {params.dcovg} \
        --intervals {params.target} \
        --dbsnp {params.dbsnp} \
        --spark-master local[{threads}] \
        --bamOutput {output.bam} \
        > {log} 2>&1
        """

rule annotate_variants:
    message:
        "Annotating variants for sample {wildcards.sample}"
    input:
        vcf=rules.haplotypecaller.output.gvcf
    output:
        annotated_vcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.annotated.vcf"
    conda:
        "gatk"
    threads:
        config["threads"]
    params:
        ref=config["reference_genome"],
        dbsnp=config["dbsnp"],
        target=config["icc_panel"]
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_annotate_variants.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_annotate_variants.txt"
    shell:
        """
        gatk VariantAnnotator \
        -R {params.ref} \
        -V {input.vcf} \
        -O {output.annotated_vcf} \
        -A Coverage \
        -A AlleleBalance \
        -A HaplotypeScore \
        -A InbreedingCoeff \
        -A HomopolymerRun \
        -A HardyWeinberg \
        -A GCContent \
        --dbsnp {params.dbsnp} \
        --intervals {params.target} \
        > {log} 2>&1
        """

rule genotype_gvcfs:
    message:
        "Genotyping GVCFs for sample {wildcards.sample}"
    input:
        gvcf=rules.annotate_variants.output.annotated_vcf
    output:
        vcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.genotyped.vcf"
    conda:
        "icc_gatk"
    threads:
        config["threads"]
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
        -A AlleleBalance \
        -A DepthPerAlleleBySample \
        -A Coverage \
        --intervals {params.target} \
        --dbsnp {params.dbsnp} \
        > {log} 2>&1
        """