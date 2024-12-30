# A rule to align the trimmed reads to the reference genome using bwa

rule bwa_alignment:
    message:
        "Aligning sample {wildcards.sample}_{lane} to the reference genome"
    input:
        fq1=rules.trimming.output.fq1,
        fq2=rules.trimming.output.fq2
    output:
        bam=config["outdir"] + "/analysis/003_alignment/bwa/{sample}_{lane}.bam"
    conda:
        "icc_04_alignment"
    threads:
        config["threads"]
    params: 
        idx=config["reference_genome"]
    log:
        config["outdir"] + "/logs/003_alignment/alignment/{sample}_{lane}_alignment.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/alignment/{sample}_{lane}_alignment.txt"
    shell:
        """
        bwa mem \
        -t {threads} \
        {params.idx} \
        {input.fq1} \
        {input.fq2} \
        | samtools view -bS - \
        > {output} \
        &> {log}
        """

rule sort_bam:
    message: 
        "Sorting aligned sample {wildcards.sample}_{lane}"
    input:
        bam=config["outdir"] + "/analysis/003_alignment/bwa/{sample}_{lane}/{sample}_{lane}.bam"
    output:
        sorted_bam=config["outdir"] + "/analysis/003_alignment/bwa/{sample}_{lane}/{sample}_{lane}_Aligned.sortedByCoord.out.bam"
    conda: 
        "icc_04_alignment"
    threads: 
        config["threads"]
    log:
        config["outdir"] + "/logs/003_alignment/sort/{sample}_{lane}_sort.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/sort/{sample}_{lane}_sort.txt"
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
        sorted_bam=config["outdir"] + "/analysis/003_alignment/bwa/{sample}_{lane}/{sample}_{lane}_Aligned.sortedByCoord.out.bam"
    output:
        indexed_bam=config["outdir"] + "/analysis/003_alignment/bwa/{sample}_{lane}/{sample}_{lane}_Aligned.sortedByCoord.out.bam.bai"
    conda:
        "icc_04_alignment"
    log:
        config["outdir"] + "/logs/003_alignment/index/{sample}_{lane}_index.log"
    shell:
        """
        samtools index {input.sorted_bam} > {log} 2>&1
        """