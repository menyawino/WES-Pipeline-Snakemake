# A rule to call variants using GATK HaplotypeCaller and UnifiedGenotyper

rule add_read_groups:
    message:
        "Adding read groups to BAM for sample {wildcards.sample}_{lane}"
    input:
        sorted_bam=rules.sort_bam.output.sorted_bam
    output:
        rg_bam=config["outdir"] + "/analysis/004_variant_calling/{sample}_{lane}.rg.bam"
    conda:
        "gatk"

    params:
        rgid=config["gatk"]["AddOrReplaceReadGroups"]["RGID"],
        rglb=config["gatk"]["AddOrReplaceReadGroups"]["RGLB"],
        rgpl=config["gatk"]["AddOrReplaceReadGroups"]["RGPL"],
        rgpu=config["gatk"]["AddOrReplaceReadGroups"]["RGPU"],
        rgsm=config["gatk"]["AddOrReplaceReadGroups"]["RGSM"],
        rgcn=config["gatk"]["AddOrReplaceReadGroups"]["RGCN"],
        rgds=config["gatk"]["AddOrReplaceReadGroups"]["RGDS"],
        validation_stringency=config["gatk"]["AddOrReplaceReadGroups"]["validation_stringency"]
    log:
        config["outdir"] + "/logs/004_variant_calling/{sample}_{lane}_add_rg.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_variant_calling/{sample}_{lane}_add_rg.txt"
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
        RGDS={params.rgds} \
        VALIDATION_STRINGENCY={params.validation_stringency} \
        > {log} 2>&1
        """


rule mark_duplicates:
    message:
        "Marking duplicates in BAM for sample {wildcards.sample}_{lane}"
    input:
        rg_bam=rules.add_read_groups.output.rg_bam
    output:
        markdup_bam=config["outdir"] + "/analysis/004_variant_calling/{sample}_{lane}.markdup.bam",
        metrics=config["outdir"] + "/analysis/004_variant_calling/{sample}_{lane}.markdup.metrics.txt"
    conda:
        "gatk"
    log:
        config["outdir"] + "/logs/004_variant_calling/{sample}_{lane}_markdup.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_variant_calling/{sample}_{lane}_markdup.txt"
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
        markdup_bam=rules.mark_duplicates.output.markdup_bam
    output:
        indexed_markdup_bam=config["outdir"] + "/analysis/004_variant_calling/{sample}_{lane}.markdup.bam.bai"
    conda:
        "gatk"
    log:
        config["outdir"] + "/logs/004_variant_calling/{sample}_{lane}_index_markdup.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_variant_calling/{sample}_{lane}_index_markdup.txt"
    shell:
        """
        samtools index \
        {input.markdup_bam} \
        > {output.indexed_markdup_bam}
        """

rule realigner_target_creator:
    message:
        "Creating realignment targets for sample {wildcards.sample}_{lane}"
    input:
        markdup_bam=rules.mark_duplicates.output.markdup_bam
    output:
        intervals=config["outdir"] + "/analysis/004_variant_calling/{sample}_{lane}.markdup.intervals"
    conda:
        "gatk"
    params:
        known=config["gatk"]["BaseRecalibrator"]["known_sites"],
        target=config["gatk"]["BaseRecalibrator"]["target"]
    log:
        config["outdir"] + "/logs/004_variant_calling/{sample}_{lane}_realigner_target_creator.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_variant_calling/{sample}_{lane}_realigner_target_creator.txt"
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
        markdup_bam=rules.mark_duplicates.output.markdup_bam
        intervals=config["outdir"] + "/analysis/004_variant_calling/{sample}_{lane}.markdup.intervals"
    output:
        realigned_bam=config["outdir"] + "/analysis/004_variant_calling/{sample}_{lane}.realigned.bam"
    conda:
        "gatk"
    params:
        known=config["gatk"]["BaseRecalibrator"]["known_sites"],
        target=config["gatk"]["BaseRecalibrator"]["target"]
    log:
        config["outdir"] + "/logs/004_variant_calling/{sample}_{lane}_indel_realigner.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_variant_calling/{sample}_{lane}_indel_realigner.txt"
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
        bam=rules.indel_realigner.output.realigned_bam
    output:
        vcf=config["outdir"] + "/analysis/004_variant_calling/{sample}_{lane}.haplotypecaller.vcf"
    conda:
        "gatk"
    threads:
        config["threads"]
    params:
        ref=config["reference_genome"]
    log:
        config["outdir"] + "/logs/004_variant_calling/{sample}_{lane}_haplotypecaller.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_variant_calling/{sample}_{lane}_haplotypecaller.txt"
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
        bam=rules.indel_realigner.output.realigned_bam
    output:
        vcf=config["outdir"] + "/analysis/004_variant_calling/{sample}_{lane}.unifiedgenotyper.vcf"
    conda:
        "gatk"
    threads:
        config["threads"]
    params:
        ref=config["reference_genome"]
    log:
        config["outdir"] + "/logs/004_variant_calling/{sample}_{lane}_unifiedgenotyper.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_variant_calling/{sample}_{lane}_unifiedgenotyper.txt"
    shell:
        """
        gatk UnifiedGenotyper \
        -R {params.ref} \
        -I {input.bam} \
        -O {output.vcf} \
        > {log} 2>&1
        """
