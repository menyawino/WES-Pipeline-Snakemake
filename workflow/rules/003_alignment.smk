# A rule to align the trimmed reads to the reference genome using bwa

rule bwa_alignment:
    message:
        "Aligning sample {wildcards.sample}_{lane} to the reference genome"
    input:
        fq1=rules.trimming.output.fq1,
        fq2=rules.trimming.output.fq2
    output:
        sam=config["outdir"] + "/analysis/003_alignment/bwa/{sample}_{lane}.sam"
    conda:
        "icc_04_alignment"
    threads:
        config["threads"]
    params: 
        idx=config["reference_genome"]
    log:
        config["outdir"] + "/logs/003_alignment/bwa/{sample}_{lane}_alignment.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/bwa/{sample}_{lane}_alignment.txt"
    shell:
        """
        bwa mem \
        -t {threads} \
        {params.idx} \
        {input.fq1} \
        {input.fq2} \
        > {output.sam}
        """

rule sam_to_bam:
    message:
        "Converting SAM to BAM for sample {wildcards.sample}_{lane}"
    input:
        sam=rules.bwa_alignment.output.sam
    output:
        bam=config["outdir"] + "/analysis/003_alignment/bwa/{sample}_{lane}.bam"
    conda:
        "icc_04_alignment"
    threads:
        config["threads"]
    log:
        config["outdir"] + "/logs/003_alignment/bwa/{sample}_{lane}_sam_to_bam.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/bwa/{sample}_{lane}_sam_to_bam.txt"
    shell:
        """
        samtools view \
        -bS {input.sam} \
        > {output.bam}
        """

rule sort_bam:
    message: 
        "Sorting aligned sample {wildcards.sample}_{lane}"
    input:
        bam=rules.sam_to_bam.output.bam
    output:
        sorted_bam=config["outdir"] + "/analysis/003_alignment/sorted/{sample}_{lane}_Aligned.sortedByCoord.bam"
    conda: 
        "icc_04_alignment"
    threads: 
        config["threads"]
    log:
        config["outdir"] + "/logs/003_alignment/sorted/{sample}_{lane}_sort.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/sorted/{sample}_{lane}_sort.txt"
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
        sorted_bam=rules.sort_bam.output.sorted_bam
    output:
        indexed_bam=config["outdir"] + "/analysis/003_alignment/sorted/{sample}_{lane}_Aligned.sortedByCoord.bam.bai"
    conda:
        "icc_04_alignment"
    log:
        config["outdir"] + "/logs/003_alignment/sorted/{sample}_{lane}_index.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/sorted/{sample}_{lane}_index.txt"
    shell:
        """
        samtools index \
        {input.sorted_bam} \
        > {log} 2>&1
        """