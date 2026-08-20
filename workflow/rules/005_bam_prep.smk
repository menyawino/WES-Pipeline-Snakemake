rule mark_duplicates:
    message:
        "Marking duplicates in BAM for sample {wildcards.sample}"
    input:
        bam=rules.merge_bams.output.merged_bam
    output:
        markdup_bam=temp(config["outdir"] + "/analysis/003_alignment/04_markduped/{sample}.markdup.bam"),
        metrics=config["outdir"] + "/analysis/003_alignment/04_markduped/{sample}.markdup.metrics.txt"
    conda:
        "icc_04_alignment"
    threads:
        config["threads_mid"]
    resources:
        mem_mb=config.get("mem_mid", 16384),
        tmpdir=config.get("tmpdir", "/tmp")
    log:
        config["outdir"] + "/logs/003_alignment/04_markduped/{sample}_markdup.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/04_markduped/{sample}_markdup.txt"
    shell:
        """
        mkdir -p "{resources.tmpdir}"
        sambamba markdup \
        -t {threads} \
        --hash-table-size=1048576 \
        --overflow-list-size=1048576 \
        --io-buffer-size=256 \
        --tmpdir "/dev/shm" \
        "{input.bam}" \
        "{output.markdup_bam}" \
        > "{log}" 2>&1
        
        touch "{output.metrics}"
        """

rule index_markdup_bam:
    message:
        "Indexing markdup BAM for sample {wildcards.sample}"
    input:
        markdup_bam=rules.mark_duplicates.output.markdup_bam
    output:
        indexed_markdup_bam=temp(config["outdir"] + "/analysis/003_alignment/04_markduped/{sample}.markdup.bam.bai")
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    log:
        config["outdir"] + "/logs/003_alignment/04_markduped/{sample}_index_markdup.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/04_markduped/{sample}_index_markdup.txt"
    shell:
        """
        samtools index -@ {threads} \
        "{input.markdup_bam}" \
        "{output.indexed_markdup_bam}" \
        > "{log}" 2>&1
        """

rule base_recalibrator:
    message:
        "Running BaseRecalibrator for sample {wildcards.sample}"
    input:
        bam=rules.mark_duplicates.output.markdup_bam,
        bai=rules.index_markdup_bam.output.indexed_markdup_bam
    output:
        recal_table=config["outdir"] + "/analysis/003_alignment/05_bqsr/{sample}.recal_data.table"
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    resources:
        mem_mb=config.get("mem_mid", 16384),
        tmpdir=config.get("tmpdir", "/tmp")
    params:
        ref=config["reference_genome"],
        known_sites=config["dbsnp"],
        target=config["icc_panel"]
    log:
        config["outdir"] + "/logs/003_alignment/05_bqsr/{sample}_base_recalibrator.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/05_bqsr/{sample}_base_recalibrator.txt"
    shell:
        """
        mkdir -p "{resources.tmpdir}"
        gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads={threads} -XX:ConcGCThreads={threads} -Xmx{resources.mem_mb}m" BaseRecalibrator \
        -I "{input.bam}" \
        -R "{params.ref}" \
        -O "{output.recal_table}" \
        --known-sites "{params.known_sites}" \
        -L "{params.target}" \
        --interval-padding 100 \
        --tmp-dir "/dev/shm" \
        > "{log}" 2>&1
        """

rule apply_bqsr:
    message:
        "Applying BQSR for sample {wildcards.sample}"
    input:
        bam=rules.mark_duplicates.output.markdup_bam,
        bai=rules.index_markdup_bam.output.indexed_markdup_bam,
        recal_table=rules.base_recalibrator.output.recal_table
    output:
        bqsr_bam=temp(config["outdir"] + "/analysis/003_alignment/05_bqsr/{sample}.bqsr.bam"),
        bqsr_bai=temp(config["outdir"] + "/analysis/003_alignment/05_bqsr/{sample}.bqsr.bai")
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    resources:
        mem_mb=config.get("mem_mid", 16384),
        tmpdir=config.get("tmpdir", "/tmp")
    params:
        ref=config["reference_genome"],
        target=config["icc_panel"]
    log:
        config["outdir"] + "/logs/003_alignment/05_bqsr/{sample}_apply_bqsr.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/05_bqsr/{sample}_apply_bqsr.txt"
    shell:
        """
        mkdir -p "{resources.tmpdir}"
        gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads={threads} -XX:ConcGCThreads={threads} -Xmx{resources.mem_mb}m" ApplyBQSR \
        -R "{params.ref}" \
        -I "{input.bam}" \
        --bqsr-recal-file "{input.recal_table}" \
        -L "{params.target}" \
        --interval-padding 100 \
        -O "{output.bqsr_bam}" \
        --create-output-bam-index true \
        --tmp-dir "/dev/shm" \
        > "{log}" 2>&1
        """

rule filter_bam_target:
    message:
        "Filtering BAM for target regions for sample {wildcards.sample}"
    input:
        bam=rules.apply_bqsr.output.bqsr_bam,
        bai=rules.apply_bqsr.output.bqsr_bai
    output:
        bam_target=config["outdir"] + "/analysis/003_alignment/06_filtering/{sample}.target.bam",
        bai_target=config["outdir"] + "/analysis/003_alignment/06_filtering/{sample}.target.bam.bai"
    conda:
        "icc_04_alignment"
    threads:
        config["threads_mid"]
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
        -L "{params.TargetFile}" \
        -F "mapping_quality >= 8" \
        -f bam \
        -o "{output.bam_target}" \
        "{input.bam}" \
        2> "{log}"
        
        sambamba index -t {threads} "{output.bam_target}" "{output.bai_target}" >> "{log}" 2>&1
        """

rule filter_bam_prot_coding:
    message:
        "Filtering BAM for protein coding regions for sample {wildcards.sample}"
    input:
        bam=rules.apply_bqsr.output.bqsr_bam,
        bai=rules.apply_bqsr.output.bqsr_bai
    output:
        bam_prot_coding=config["outdir"] + "/analysis/003_alignment/06_filtering/{sample}.prot_coding.bam",
        bai_prot_coding=config["outdir"] + "/analysis/003_alignment/06_filtering/{sample}.prot_coding.bam.bai"
    conda:
        "icc_04_alignment"
    threads:
        config["threads_mid"]
    params:
        CDSFile=config["cds_panel"],
        TargetFile=config["icc_panel"],
        target_bam=rules.filter_bam_target.output.bam_target,
        target_bai=rules.filter_bam_target.output.bai_target
    log:
        config["outdir"] + "/logs/003_alignment/06_filtering/{sample}_filter_bam_prot_coding.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/06_filtering/{sample}_filter_bam_prot_coding.txt"
    shell:
        """
        if [ "{params.CDSFile}" = "{params.TargetFile}" ] && [ -f "{params.target_bam}" ]; then
            ln -f "{params.target_bam}" "{output.bam_prot_coding}" 2>/dev/null || cp "{params.target_bam}" "{output.bam_prot_coding}"
            ln -f "{params.target_bai}" "{output.bai_prot_coding}" 2>/dev/null || cp "{params.target_bai}" "{output.bai_prot_coding}"
        else
            sambamba view \
            -t {threads} \
            -L "{params.CDSFile}" \
            -F "mapping_quality >= 8" \
            -f bam \
            -o "{output.bam_prot_coding}" \
            "{input.bam}" \
            2> "{log}"
            
            sambamba index -t {threads} "{output.bam_prot_coding}" "{output.bai_prot_coding}" >> "{log}" 2>&1
        fi
        """

rule filter_bam_canon_tran:
    message:
        "Filtering BAM for canonical transcript regions for sample {wildcards.sample}"
    input:
        bam=rules.apply_bqsr.output.bqsr_bam,
        bai=rules.apply_bqsr.output.bqsr_bai
    output:
        bam_canon_tran=config["outdir"] + "/analysis/003_alignment/06_filtering/{sample}.canon_tran.bam",
        bai_canon_tran=config["outdir"] + "/analysis/003_alignment/06_filtering/{sample}.canon_tran.bam.bai"
    conda:
        "icc_04_alignment"
    threads:
        config["threads_mid"]
    params:
        CanonTranFile=config["canontran_panel"],
        TargetFile=config["icc_panel"],
        target_bam=rules.filter_bam_target.output.bam_target,
        target_bai=rules.filter_bam_target.output.bai_target
    log:
        config["outdir"] + "/logs/003_alignment/06_filtering/{sample}_filter_bam_canon_tran.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/06_filtering/{sample}_filter_bam_canon_tran.txt"
    shell:
        """
        if [ "{params.CanonTranFile}" = "{params.TargetFile}" ] && [ -f "{params.target_bam}" ]; then
            ln -f "{params.target_bam}" "{output.bam_canon_tran}" 2>/dev/null || cp "{params.target_bam}" "{output.bam_canon_tran}"
            ln -f "{params.target_bai}" "{output.bai_canon_tran}" 2>/dev/null || cp "{params.target_bai}" "{output.bai_canon_tran}"
        else
            sambamba view \
            -t {threads} \
            -L "{params.CanonTranFile}" \
            -F "mapping_quality >= 8" \
            -f bam \
            -o "{output.bam_canon_tran}" \
            "{input.bam}" \
            2> "{log}"
            
            sambamba index -t {threads} "{output.bam_canon_tran}" "{output.bai_canon_tran}" >> "{log}" 2>&1
        fi
        """