rule add_read_groups:
    message:
        "Adding read groups to BAM for sample {wildcards.sample}"
    input:
        bam=rules.merge_bams.output.merged_bam
    output:
        rg_bam=config["outdir"] + "/analysis/003_alignment/03_read_grouped/{sample}.rg.bam"
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
        config["outdir"] + "/logs/003_alignment/03_read_grouped/{sample}_add_rg.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/03_read_grouped/{sample}_add_rg.txt"
    shell:
        """
        gatk AddOrReplaceReadGroups \
        I={input.bam} \
        O={output.rg_bam} \
        RGID={params.rgid} \
        RGLB={params.rglb} \
        RGPL={params.rgpl} \
        RGPU={params.rgpu} \
        RGSM={params.rgsm} \
        RGCN={params.rgcn} \
        RGDS={params.rgds} \
        VALIDATION_STRINGENCY={params.validation_stringency} \
        &> {log}
        """

rule mark_duplicates:
    message:
        "Marking duplicates in BAM for sample {wildcards.sample}"
    input:
        rg_bam=rules.add_read_groups.output.rg_bam
    output:
        markdup_bam=config["outdir"] + "/analysis/003_alignment/04_markduped/{sample}.markdup.bam",
        metrics=config["outdir"] + "/analysis/003_alignment/04_markduped/{sample}.markdup.metrics.txt"
    conda:
        "icc_gatk"
    params:
        validation_stringency=config["gatk"]["MarkDuplicates"]["validation_stringency"]
    threads:
        config["threads_high"]
    log:
        config["outdir"] + "/logs/003_alignment/04_markduped/{sample}_markdup.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/04_markduped/{sample}_markdup.txt"
    shell:
        """
        gatk MarkDuplicatesSpark \
        -I {input.rg_bam} \
        -O {output.markdup_bam} \
        -M {output.metrics} \
        -VS {params.validation_stringency} \
        --spark-master local[{threads}] \
        &> {log}
        """

rule base_recalibrator:
    message:
        "Running BaseRecalibrator for sample {wildcards.sample}"
    input:
        bam=rules.mark_duplicates.output.markdup_bam,
    output:
        recal_table=config["outdir"] + "/analysis/003_alignment/05_bqsr/{sample}.recal_data.table"
    conda:
        "icc_gatk"
    threads:
        config["threads_high"]
    params:
        ref=config["reference_genome"],
        dbsnp=config["dbsnp"],
        omni=config["omni"],
        tenk=config["tenk"],
        hapmap=config["hapmap"],
        mills=config["mills"],
        tenk_indel=config["tenk_indel"]
    log:
        config["outdir"] + "/logs/003_alignment/05_bqsr/{sample}_base_recalibrator.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/05_bqsr/{sample}_base_recalibrator.txt"
    shell:
        """
        gatk BaseRecalibratorSpark \
        -I {input.bam} \
        -R {params.ref} \
        -O {output.recal_table} \
        --known-sites {params.dbsnp} \
        --known-sites {params.mills} \
        --known-sites {params.tenk_indel} \
        --spark-master local[{threads}] \
        &> {log}
        """

rule apply_bqsr:
    message:
        "Applying BQSR for sample {wildcards.sample}"
    input:
        bam=rules.mark_duplicates.output.markdup_bam,
        recal_table=rules.base_recalibrator.output.recal_table
    output:
        bqsr_bam=config["outdir"] + "/analysis/003_alignment/05_bqsr/{sample}.bqsr.bam"
    conda:
        "icc_gatk"
    threads:
        config["threads_high"]
    params:
        ref=config["reference_genome"]
    log:
        config["outdir"] + "/logs/003_alignment/05_bqsr/{sample}_apply_bqsr.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/05_bqsr/{sample}_apply_bqsr.txt"
    shell:
        """
        gatk ApplyBQSRSpark  \
        -R {params.ref} \
        -I {input.bam} \
        --bqsr-recal-file {input.recal_table} \
        -O {output.bqsr_bam} \
        --spark-master local[{threads}] \
        &> {log}
        """

rule filter_bam_target:
    message:
        "Filtering BAM for target regions for sample {wildcards.sample}"
    input:
        bam=rules.apply_bqsr.output.bqsr_bam
    output:
        bam_target=config["outdir"] + "/analysis/003_alignment/06_filtering/{sample}.target.bam"
    conda:
        "icc_gatk"
    threads:
        config["threads_high"]
    params:
        TargetFile=config["icc_panel"]
    log:
        config["outdir"] + "/logs/003_alignment/06_filtering/{sample}_filter_bam_target.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/06_filtering/{sample}_filter_bam_target.txt"
    shell:
        """
        sambamba view \
        -t {threads} \
        -L {params.TargetFile} \
        -f bam -F "mapping_quality > 8" \
        {input.bam} \
        -o {output.bam_target} \
        2> {log}
        """

rule filter_bam_prot_coding:
    message:
        "Filtering BAM for protein coding regions for sample {wildcards.sample}"
    input:
        bam=rules.apply_bqsr.output.bqsr_bam
    output:
        bam_prot_coding=config["outdir"] + "/analysis/003_alignment/06_filtering/{sample}.prot_coding.bam"
    conda:
        "icc_gatk"
    threads:
        config["threads_high"]
    params:
        CDSFile=config["cds_panel"]
    log:
        config["outdir"] + "/logs/003_alignment/06_filtering/{sample}_filter_bam_prot_coding.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/06_filtering/{sample}_filter_bam_prot_coding.txt"
    shell:
        """
        sambamba view \
        -t {threads} \
        -L {params.CDSFile} \
        -f bam -F "mapping_quality > 8" \
        {input.bam} \
        -o {output.bam_prot_coding} \
        2> {log}
        """

rule filter_bam_canon_tran:
    message:
        "Filtering BAM for canonical transcript regions for sample {wildcards.sample}"
    input:
        bam=rules.apply_bqsr.output.bqsr_bam
    output:
        bam_canon_tran=config["outdir"] + "/analysis/003_alignment/06_filtering/{sample}.canon_tran.bam"
    conda:
        "icc_gatk"
    threads:
        config["threads_high"]
    params:
        CanonTranFile=config["canontran_panel"]
    log:
        config["outdir"] + "/logs/003_alignment/06_filtering/{sample}_filter_bam_canon_tran.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/06_filtering/{sample}_filter_bam_canon_tran.txt"
    shell:
        """
        sambamba view \
        -t {threads} \
        -L {params.CanonTranFile} \
        -f bam -F "mapping_quality > 8" \
        {input.bam} \
        -o {output.bam_canon_tran} \
        2> {log}
        """