# A rule to call variants using GATK HaplotypeCaller and UnifiedGenotyper

rule add_read_groups:
    message:
        "Adding read groups to BAM for sample {wildcards.sample}_{lane}"
    input:
        sorted_bam=lambda wildcards: config["outdir"] + "/analysis/004_alignment/bwa/{sample}_{lane}/{sample}_{lane}_Aligned.sortedByCoord.out.bam".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        rg_bam=lambda wildcards: config["outdir"] + "/analysis/005_variant_calling/{sample}_{lane}/{sample}_{lane}.rg.bam".format(sample=wildcards.sample, lane=wildcards.lane)
    conda:
        "gatk"
    params:
        rgid=config["rgid"],
        rglb=config["rglb"],
        rgpl=config["rgpl"],
        rgpu=config["rgpu"],
        rgsm=config["rgsm"],
        rgcn=config["rgcn"]
    log:
        lambda wildcards: config["outdir"] + "/logs/005_variant_calling/{sample}_{lane}_add_rg.log".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        gatk AddOrReplaceReadGroups \
        I={input.sorted_bam} \
        O={output.rg_bam} \
        RGID={params.rgid} \
        RGLB={params.rglb} \
        RGPL={params.rgpl} \
        RGPU={params.rgpu} \
        RGSM={params.rgsm} \
        RGCN={params.rgcn} \
        > {log} 2>&1
        """

rule mark_duplicates:
    message:
        "Marking duplicates in BAM for sample {wildcards.sample}_{lane}"
    input:
        rg_bam=lambda wildcards: config["outdir"] + "/analysis/005_variant_calling/{sample}_{lane}/{sample}_{lane}.rg.bam".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        markdup_bam=lambda wildcards: config["outdir"] + "/analysis/005_variant_calling/{sample}_{lane}/{sample}_{lane}.markdup.bam".format(sample=wildcards.sample, lane=wildcards.lane),
        metrics=lambda wildcards: config["outdir"] + "/analysis/005_variant_calling/{sample}_{lane}/{sample}_{lane}.markdup.metrics.txt".format(sample=wildcards.sample, lane=wildcards.lane)
    conda:
        "gatk"
    log:
        lambda wildcards: config["outdir"] + "/logs/005_variant_calling/{sample}_{lane}_markdup.log".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        gatk MarkDuplicates \
        I={input.rg_bam} \
        O={output.markdup_bam} \
        M={output.metrics} \
        > {log} 2>&1
        """

rule index_markdup_bam:
    message:
        "Indexing markdup BAM for sample {wildcards.sample}_{lane}"
    input:
        markdup_bam=lambda wildcards: config["outdir"] + "/analysis/005_variant_calling/{sample}_{lane}/{sample}_{lane}.markdup.bam".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        indexed_markdup_bam=lambda wildcards: config["outdir"] + "/analysis/005_variant_calling/{sample}_{lane}/{sample}_{lane}.markdup.bam.bai".format(sample=wildcards.sample, lane=wildcards.lane)
    conda:
        "gatk"
    log:
        lambda wildcards: config["outdir"] + "/logs/005_variant_calling/{sample}_{lane}_index_markdup.log".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        samtools index {input.markdup_bam} > {log} 2>&1
        """

rule realigner_target_creator:
    message:
        "Creating realignment targets for sample {wildcards.sample}_{lane}"
    input:
        markdup_bam=lambda wildcards: config["outdir"] + "/analysis/005_variant_calling/{sample}_{lane}/{sample}_{lane}.markdup.bam".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        intervals=lambda wildcards: config["outdir"] + "/analysis/005_variant_calling/{sample}_{lane}/{sample}_{lane}.markdup.intervals".format(sample=wildcards.sample, lane=wildcards.lane)
    conda:
        "gatk"
    params:
        known=config["known_sites"],
        target=config["target_file"]
    log:
        lambda wildcards: config["outdir"] + "/logs/005_variant_calling/{sample}_{lane}_realigner_target_creator.log".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        gatk RealignerTargetCreator \
        -R {params.target} \
        -I {input.markdup_bam} \
        -known {params.known} \
        -o {output.intervals} \
        > {log} 2>&1
        """

rule indel_realigner:
    message:
        "Realigning indels for sample {wildcards.sample}_{lane}"
    input:
        markdup_bam=lambda wildcards: config["outdir"] + "/analysis/005_variant_calling/{sample}_{lane}/{sample}_{lane}.markdup.bam".format(sample=wildcards.sample, lane=wildcards.lane),
        intervals=lambda wildcards: config["outdir"] + "/analysis/005_variant_calling/{sample}_{lane}/{sample}_{lane}.markdup.intervals".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        realigned_bam=lambda wildcards: config["outdir"] + "/analysis/005_variant_calling/{sample}_{lane}/{sample}_{lane}.realigned.bam".format(sample=wildcards.sample, lane=wildcards.lane)
    conda:
        "gatk"
    log:
        lambda wildcards: config["outdir"] + "/logs/005_variant_calling/{sample}_{lane}_indel_realigner.log".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        gatk IndelRealigner \
        -R {params.target} \
        -I {input.markdup_bam} \
        -targetIntervals {input.intervals} \
        -o {output.realigned_bam} \
        > {log} 2>&1
        """

rule haplotypecaller:
    message:
        "Calling variants with HaplotypeCaller for sample {wildcards.sample}_{lane}"
    input:
        bam=lambda wildcards: config["outdir"] + "/analysis/005_variant_calling/{sample}_{lane}/{sample}_{lane}.realigned.bam".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        vcf=lambda wildcards: config["outdir"] + "/analysis/005_variant_calling/{sample}_{lane}/{sample}_{lane}.haplotypecaller.vcf".format(sample=wildcards.sample, lane=wildcards.lane)
    conda:
        "gatk"
    threads:
        config["threads"]
    params:
        ref=config["reference"]
    log:
        lambda wildcards: config["outdir"] + "/logs/005_variant_calling/{sample}_{lane}_haplotypecaller.log".format(sample=wildcards.sample, lane=wildcards.lane)
    benchmark:
        lambda wildcards: config["outdir"] + "/benchmarks/005_variant_calling/{sample}_{lane}_haplotypecaller.txt".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        gatk HaplotypeCaller \
        -R {params.ref} \
        -I {input.bam} \
        -O {output.vcf} \
        -ERC GVCF \
        > {log} 2>&1
        """

rule unifiedgenotyper:
    message:
        "Calling variants with UnifiedGenotyper for sample {wildcards.sample}_{lane}"
    input:
        bam=lambda wildcards: config["outdir"] + "/analysis/005_variant_calling/{sample}_{lane}/{sample}_{lane}.realigned.bam".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        vcf=lambda wildcards: config["outdir"] + "/analysis/005_variant_calling/{sample}_{lane}/{sample}_{lane}.unifiedgenotyper.vcf".format(sample=wildcards.sample, lane=wildcards.lane)
    conda:
        "gatk"
    threads:
        config["threads"]
    params:
        ref=config["reference"]
    log:
        lambda wildcards: config["outdir"] + "/logs/005_variant_calling/{sample}_{lane}_unifiedgenotyper.log".format(sample=wildcards.sample, lane=wildcards.lane)
    benchmark:
        lambda wildcards: config["outdir"] + "/benchmarks/005_variant_calling/{sample}_{lane}_unifiedgenotyper.txt".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        gatk UnifiedGenotyper \
        -R {params.ref} \
        -I {input.bam} \
        -O {output.vcf} \
        > {log} 2>&1
        """
