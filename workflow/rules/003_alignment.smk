# A rule to align the trimmed reads to the reference genome using bwa and convert to BAM using sambamba view

rule bwa_mem:
    message:
        "Aligning and converting to BAM for sample {wildcards.sample}_{lane}"
    input:
        fq1=rules.trimming_fp.output.fq1,
        fq2=rules.trimming_fp.output.fq2
    output:
        bam=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_{lane}.bam"
    conda:
        "icc_04_alignment"
    threads:
        config["threads_high"]
    params: 
        ref=config["reference_genome"]
    log:
        bwa=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_{lane}_bwa.log",
        sambamba_view=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_{lane}_sambamba_view.log",
        sambamba_sort=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_{lane}_sambamba_sort.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_{lane}_alignment.txt"
    shell:
        """
        bwa-mem2 mem \
        -t {threads} \
        {params.ref} \
        {input.fq1} \
        {input.fq2} \
        2> {log.bwa} \
        | sambamba view \
        -S -f bam \
        -t {threads} \
        -l 0 \
        /dev/stdin \
        2> {log.sambamba_view} |
        sambamba sort \
        -t {threads} \
        -o {output.bam} \
        /dev/stdin \
        2> {log.sambamba_sort}
        """

rule merge_bams:
    message:
        "Merging BAM files for sample {wildcards.sample}"
    input:
        bams=expand(rules.bwa_mem.output.bam, sample="{sample}", lane=("L001", "L002", "L003", "L004")),
    output:
        merged_bam=config["outdir"] + "/analysis/003_alignment/02_merged/{sample}.merged.bam"
    conda:
        "icc_04_alignment"
    threads:
        config["threads_mid"]
    log:
        config["outdir"] + "/logs/003_alignment/02_merged/{sample}_merge.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/02_merged/{sample}_merge.txt"
    shell:
        """
        sambamba merge \
        -t {threads} \
        {output.merged_bam} \
        {input.bams} \
        > {log} 2>&1
        """
