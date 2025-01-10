# A rule to call variants using GATK HaplotypeCaller and GenotypeGVCFs

rule haplotypecaller:
    message:
        "Calling variants with HaplotypeCaller for sample {wildcards.sample}"
    input:
        bam=rules.apply_bqsr.output.bqsr_bam
    output:
        gvcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.haplotypecaller.g.vcf"
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
	    -A AlleleBalance \
        -A Coverage \
        -A HaplotypeScore \
        -A InbreedingCoeff \
        -A HomopolymerRun \
        -A HardyWeinberg \
        -A GCContent \
	    -G Standard \
        -ERC GVCF \
        --downsample_to_coverage {params.dcovg} \
        --intervals {params.target} \
        --dbsnp {params.dbsnp} \
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
        -A HaplotypeScore \
        -A InbreedingCoeff \
        -A HomopolymerRun \
        -A HardyWeinberg \
        -A GCContent \
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
        config["threads"]
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