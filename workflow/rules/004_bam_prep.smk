rule add_read_groups:
    message:
        "Adding read groups to BAM for sample {wildcards.sample}_{lane}"
    input:
        sorted_bam=rules.sort_bam.output.sorted_bam
    output:
        rg_bam=config["outdir"] + "/analysis/003_alignment/03_read_grouped/{sample}_{lane}.rg.bam"
    conda:
        "icc_gatk"
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
        config["outdir"] + "/logs/003_alignment/03_read_grouped/{sample}_{lane}_add_rg.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/03_read_grouped/{sample}_{lane}_add_rg.txt"
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
        markdup_bam=config["outdir"] + "/analysis/003_alignment/04_markduped/{sample}_{lane}.markdup.bam",
        metrics=config["outdir"] + "/analysis/003_alignment/04_markduped/{sample}_{lane}.markdup.metrics.txt"
    conda:
        "icc_gatk"
    threads:
        config["threads"]
    log:
        config["outdir"] + "/logs/003_alignment/04_markduped/{sample}_{lane}_markdup.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/04_markduped/{sample}_{lane}_markdup.txt"
    shell:
        """
        gatk MarkDuplicatesSpark \
        -I {input.rg_bam} \
        -O {output.markdup_bam} \
        -M {output.metrics} \
        --spark-master local[{threads}] \
        > {log} 2>&1
        """

rule index_markdup_bam:
    message:
        "Indexing markdup BAM for sample {wildcards.sample}_{lane}"
    input:
        markdup_bam=rules.mark_duplicates.output.markdup_bam
    output:
        indexed_markdup_bam=config["outdir"] + "/analysis/003_alignment/04_markduped/{sample}_{lane}.markdup.bam.bai"
    conda:
        "icc_gatk"
    log:
        config["outdir"] + "/logs/003_alignment/04_markduped/{sample}_{lane}_index_markdup.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/04_markduped/{sample}_{lane}_index_markdup.txt"
    shell:
        """
        samtools index \
        {input.markdup_bam} \
        > {output.indexed_markdup_bam}
        """

rule base_recalibrator:
    message:
        "Running BaseRecalibrator for sample {wildcards.sample}_{lane}"
    input:
        bam=rules.mark_duplicates.output.markdup_bam,
    output:
        recal_table=config["outdir"] + "/analysis/003_alignment/05_bqsr/{sample}_{lane}.recal_data.table"
    conda:
        "icc_gatk"
    threads:
        config["threads"]
    params:
        ref=config["reference_genome"],
        known_sites=config["gatk"]["BaseRecalibrator"]["known_sites"]
    log:
        config["outdir"] + "/logs/003_alignment/05_bqsr/{sample}_{lane}_base_recalibrator.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/05_bqsr/{sample}_{lane}_base_recalibrator.txt"
    shell:
        """
        gatk BaseRecalibratorSpark \
        -I {input.bam} \
        -R {params.ref} \
        -O {output.recal_table} \
        --known-sites {params.known_sites} \
        --spark-master local[{threads}] \
        > {log} 2>&1
        """

rule apply_bqsr:
    message:
        "Applying BQSR for sample {wildcards.sample}_{lane}"
    input:
        bam=rules.mark_duplicates.output.markdup_bam,
        recal_table=rules.base_recalibrator.output.recal_table
    output:
        bqsr_bam=config["outdir"] + "/analysis/003_alignment/05_bqsr/{sample}_{lane}.bqsr.bam"
    conda:
        "icc_gatk"
    threads:
        config["threads"]
    params:
        ref=config["reference_genome"]
    log:
        config["outdir"] + "/logs/003_alignment/05_bqsr/{sample}_{lane}_apply_bqsr.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/05_bqsr/{sample}_{lane}_apply_bqsr.txt"
    shell:
        """
        gatk ApplyBQSRSpark  \
        -R {params.ref} \
        -I {input.bam} \
        --bqsr-recal-file {input.recal_table} \
        -O {output.bqsr_bam} \
        --spark-master local[{threads}] \
        > {log} 2>&1
        """
