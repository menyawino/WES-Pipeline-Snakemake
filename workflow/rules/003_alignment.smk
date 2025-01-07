# A rule to align the trimmed reads to the reference genome using bwa

rule bwa_alignment:
    message:
        "Aligning sample {wildcards.sample}_{lane} to the reference genome"
    input:
        fq1=rules.trimming.output.fq1,
        fq2=rules.trimming.output.fq2
    output:
        bam=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_{lane}.bam"
    conda:
        "icc_04_alignment"
    threads:
        config["threads_high"]
    params: 
        ref=config["reference_genome"]
    log:
        config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_{lane}_alignment.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_{lane}_alignment.txt"
    shell:
        """
        bwa mem \
        -t {threads} \
        {params.ref} \
        {input.fq1} \
        {input.fq2} \
        2> {log} \
        | sambamba view \
        -S -f bam \
        -t {threads} \
        -l 0 \
        /dev/stdout \
        /dev/stdin \
        2> {log} \
        | sambamba sort \
        -t {threads} \
        -o {output.bam} \
        /dev/stdin \
        2> {log}
        """

        # bwa mem \
        # -t {threads} \
        # {params.ref} \
        # {input.fq1} \
        # {input.fq2} \
        # > {output.sam}


# rule index_bam:
#     message:
#         "Indexing sorted BAM for sample {wildcards.sample}_{lane}"
#     input:
#         sorted_bam=rules.sort_bam.output.sorted_bam
#     output:
#         indexed_bam=config["outdir"] + "/analysis/003_alignment/02_sorted/{sample}_{lane}_Aligned.sortedByCoord.bam.bai"
#     conda:
#         "icc_04_alignment"
#     log:
#         config["outdir"] + "/logs/003_alignment/02_sorted/{sample}_{lane}_index.log"
#     benchmark:
#         config["outdir"] + "/benchmarks/003_alignment/02_sorted/{sample}_{lane}_index.txt"
#     shell:
#         """
#         samtools index \
#         {input.sorted_bam} \
#         > {log} 2>&1
#         """

# rule merge_bams:
#     message:
#         "Merging BAM files for sample {wildcards.sample}"
#     input:
#         bams=expand(config["outdir"] + "/analysis/003_alignment/02_sorted/{sample}_L00{lane}_Aligned.sortedByCoord.bam", sample="{sample}", lane=(1, 2, 3, 4))
#     output:
#         merged_bam=config["outdir"] + "/analysis/003_alignment/03_merged/{sample}.merged.bam"
#     conda:
#         "icc_04_alignment"
#     threads:
#         config["threads"]
#     log:
#         config["outdir"] + "/logs/003_alignment/03_merged/{sample}_merge.log"
#     benchmark:
#         config["outdir"] + "/benchmarks/003_alignment/03_merged/{sample}_merge.txt"
#     shell:
#         """
#         samtools merge \
#         -@ {threads} \
#         {output.merged_bam} \
#         {input.bams} \
#         > {log} 2>&1
#         """

# rule index_merged_bam:
#     message:
#         "Indexing merged BAM for sample {wildcards.sample}"
#     input:
#         merged_bam=rules.merge_bams.output.merged_bam
#     output:
#         indexed_merged_bam=config["outdir"] + "/analysis/003_alignment/03_merged/{sample}.merged.bam.bai"
#     conda:
#         "icc_04_alignment"
#     log:
#         config["outdir"] + "/logs/003_alignment/03_merged/{sample}_index.log"
#     benchmark:
#         config["outdir"] + "/benchmarks/003_alignment/03_merged/{sample}_index.txt"
#     shell:
#         """
#         samtools index \
#         {input.merged_bam} \
#         > {log} 2>&1
#         """