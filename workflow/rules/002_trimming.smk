# Rule to trim and remove adapters from raw FASTQ data using fastp

rule trimming_fp:
    message: 
        "Trimming and removing adapters from sample {wildcards.sample}_{lane}"
    container: 
        "docker://quay.io/biocontainers/fastp:0.23.4--h5f740d0_0"
    input:
        fq1=get_trimming_r1,
        fq2=get_trimming_r2
    output:
        fq1=temp(config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_R1.fastq.gz"),
        fq2=temp(config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_R2.fastq.gz"),
        report=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_report.html",
        json=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_report.json"
    threads:
        config["threads_mid"]
    resources:
        mem_mb=config.get("mem_mid", 16384)
    params:
        path=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}",
        min_length=config["fastp"]["min_read_length"],
        window_size=config["fastp"]["window_size"]
    log:
        config["outdir"] + "/logs/002_trimming/{sample}/{sample}_{lane}.log"
    benchmark:
        config["outdir"] + "/benchmarks/002_trimming/{sample}/{sample}_{lane}.txt"
    shell:
        """
        fastp \
        -i "{input.fq1}" \
        -I "{input.fq2}" \
        -j "{output.json}" \
        -o "{output.fq1}" \
        -O "{output.fq2}" \
        -h "{output.report}" \
        -w {threads} \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --length_required {params.min_length} \
        --cut_window_size {params.window_size} \
        &> "{log}"
        """


        