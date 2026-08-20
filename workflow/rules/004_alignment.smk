rule bwa_mem:
    message:
        "Aligning and converting to BAM for sample {wildcards.sample}_{lane}"
    input:
        fq1=rules.trimming_fp.output.fq1,
        fq2=rules.trimming_fp.output.fq2
    output:
        bam=temp(config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_{lane}.bam")
    conda:
        "../envs/004_alignment.yml"
    threads:
        config["threads_high"]
    resources:
        mem_mb=config.get("mem_high", 32768),
        tmpdir=config.get("tmpdir", "/tmp")
    params: 
        ref=config["reference_genome"]
    log:
        bwa=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_{lane}_bwa.log",
        sort=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_{lane}_sort.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_{lane}_alignment.txt"
    shell:
        """
        mkdir -p {resources.tmpdir}
        sample_name=$(basename {wildcards.sample})
        rg_header="@RG\\tID:${{sample_name}}_{wildcards.lane}\\tSM:${{sample_name}}\\tLB:lib1\\tPL:illumina\\tPU:unit1"
        
        bwa-mem2 mem \
        -t {threads} \
        -K 100000000 -Y \
        -R "$rg_header" \
        "{params.ref}" \
        <(pigz -dc "{input.fq1}") \
        <(pigz -dc "{input.fq2}") \
        2> "{log.bwa}" \
        | sambamba sort \
        -t {threads} \
        -m 2G \
        -l 1 \
        --tmpdir "/dev/shm/sort_${{sample_name}}_{wildcards.lane}" \
        -o "{output.bam}" \
        /dev/stdin \
        2> "{log.sort}"
        """

rule merge_bams:
    message:
        "Merging BAM files for sample {wildcards.sample}"
    input:
        bams=lambda wildcards: expand(rules.bwa_mem.output.bam, sample=wildcards.sample, lane=lane)
    output:
        merged_bam=temp(config["outdir"] + "/analysis/003_alignment/02_merged/{sample}.merged.bam")
    conda:
        "../envs/004_alignment.yml"
    threads:
        config["threads_mid"]
    resources:
        mem_mb=config.get("mem_mid", 16384)
    log:
        config["outdir"] + "/logs/003_alignment/02_merged/{sample}_merge.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/02_merged/{sample}_merge.txt"
    shell:
        """
        bam_count=$(echo {input.bams} | wc -w)
        if [ "$bam_count" -eq 1 ]; then
            ln -f {input.bams} {output.merged_bam} 2>/dev/null || cp {input.bams} {output.merged_bam}
        else
            sambamba merge \
            -t {threads} \
            {output.merged_bam} \
            {input.bams} \
            > {log} 2>&1
        fi
        """
