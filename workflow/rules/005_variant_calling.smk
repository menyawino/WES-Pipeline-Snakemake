# A rule to call variants using GATK HaplotypeCaller and UnifiedGenotyper

rule haplotypecaller:
    message:
        "Calling variants with HaplotypeCaller for sample {wildcards.sample}"
    input:
        bam=rules.apply_bqsr.output.bqsr_bam
    output:
        vcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.haplotypecaller.vcf"
    conda:
        "icc_gatk"
    threads:
        config["threads"]
    params:
        ref=config["reference_genome"]
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_haplotypecaller.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_haplotypecaller.txt"
    shell:
        """
        gatk HaplotypeCallerSpark \
        -R {params.ref} \
        -I {input.bam} \
        -O {output.vcf} \
        -ERC GVCF \
        --spark-master local[{threads}] \
        > {log} 2>&1
        """

rule unifiedgenotyper:
    message:
        "Calling variants with UnifiedGenotyper for sample {wildcards.sample}"
    input:
        bam=rules.apply_bqsr.output.bqsr_bam
    output:
        vcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.unifiedgenotyper.vcf"
    conda:
        "icc_gatk"
    threads:
        config["threads"]
    params:
        ref=config["reference_genome"]
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_unifiedgenotyper.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_unifiedgenotyper.txt"
    shell:
        """
        gatk UnifiedGenotyper \
        -R {params.ref} \
        -I {input.bam} \
        -O {output.vcf} \
        > {log} 2>&1
        """
