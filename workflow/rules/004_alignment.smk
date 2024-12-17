# A rule to align the trimmed reads to the reference genome using bwa

rule bwa_alignment:
    message:
        "Aligning sample {wildcards.sample}_{lane} to the reference genome"
    input:
        fq1=lambda wildcards: config["outdir"] + "/analysis/002_trimming/{sample}/{sample}_{lane}_R1_trimmed.fastq.gz".format(sample=wildcards.sample, lane=wildcards.lane),
        fq2=lambda wildcards: config["outdir"] + "/analysis/002_trimming/{sample}/{sample}_{lane}_R2_trimmed.fastq.gz".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        bam=lambda wildcards: config["outdir"] + "/analysis/004_alignment/bwa/{sample}_{lane}/{sample}_{lane}.bam".format(sample=wildcards.sample, lane=wildcards.lane)
    conda:
        "icc_04_alignment"
    threads:
        config["threads"]
    params: 
        path=lambda wildcards: config["outdir"] + "/results/{}".format(wildcards.sample),
        idx=config["aligner"]["index_bwa"]
    log:
        lambda wildcards: config["outdir"] + "/logs/004_alignment/alignment/{sample}_{lane}_alignment.log".format(sample=wildcards.sample, lane=wildcards.lane)
    benchmark:
        lambda wildcards: config["outdir"] + "/benchmarks/004_alignment/alignment/{sample}_{lane}_alignment.txt".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        bwa mem \
        -t {threads} \
        {params.idx} \
        {input.fq1} \
        {input.fq2} \
        | samtools view -bS - \
        > {output} \
        2> {log}
        """

rule sort_bam:
    message: 
        "Sorting aligned sample {wildcards.sample}_{lane}"
    input:
        bam=lambda wildcards: config["outdir"] + "/analysis/004_alignment/bwa/{sample}_{lane}/{sample}_{lane}.bam".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        sorted_bam=lambda wildcards: config["outdir"] + "/analysis/004_alignment/bwa/{sample}_{lane}/{sample}_{lane}_Aligned.sortedByCoord.out.bam".format(sample=wildcards.sample, lane=wildcards.lane)
    conda: 
        "icc_04_alignment"
    threads: 
        config["threads"]
    params: 
        path=lambda wildcards: config["outdir"] + "/results/{}".format(wildcards.sample)
    log:
        lambda wildcards: config["outdir"] + "/logs/004_alignment/sort/{sample}_{lane}_sort.log".format(sample=wildcards.sample, lane=wildcards.lane)
    benchmark:
        lambda wildcards: config["outdir"] + "/benchmarks/004_alignment/sort/{sample}_{lane}_sort.txt".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        samtools sort \
        -@ {threads} \
        -o {output} \
        {input.bam} \
        > {log} 2>&1
        """

rule index_bam:
    message:
        "Indexing sorted BAM for sample {wildcards.sample}_{lane}"
    input:
        sorted_bam=lambda wildcards: config["outdir"] + "/analysis/004_alignment/bwa/{sample}_{lane}/{sample}_{lane}_Aligned.sortedByCoord.out.bam".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        indexed_bam=lambda wildcards: config["outdir"] + "/analysis/004_alignment/bwa/{sample}_{lane}/{sample}_{lane}_Aligned.sortedByCoord.out.bam.bai".format(sample=wildcards.sample, lane=wildcards.lane)
    conda:
        "icc_04_alignment"
    log:
        lambda wildcards: config["outdir"] + "/logs/004_alignment/index/{sample}_{lane}_index.log".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        samtools index {input.sorted_bam} > {log} 2>&1
        """