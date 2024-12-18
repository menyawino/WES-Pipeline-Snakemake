# A rule to run FastQC on the raw data

rule raw_fastqc:
    message: 
        "Running FastQC on sample {wildcards.sample}_{lane}_{wildcards.R}"
    conda: 
        "fastqc_env"
    input:
        config["inputdir"] + "/samples/{sample}_{lane}_{R}_001.fastq.gz"
    output:
        html=config["outdir"] + "/analysis/001_QC/{sample}/{sample}_{lane}_{R}_fastqc.html",
        zip=config["outdir"] + "/analysis/001_QC/{sample}/{sample}_{lane}_{R}_fastqc.zip"
    threads: 
        config["np_threads"]
    params: 
        path=config["outdir"] + "/analysis/001_QC/{}"
    log:
        config["outdir"] + "/logs/001_QC/{sample}/{sample}_{lane}_{R}.log"
    benchmark:
        config["outdir"] + "/benchmarks/001_QC/{sample}/{sample}_{lane}_{R}.txt"
    shell:
        """
        mkdir -p {params.path}
        fastqc {input} \
        -t {threads} \
        -o {params.path} \
        > {log} 2>&1
        mv {params.path}/{wildcards.sample}_{wildcards.lane}_{wildcards.R}_001_fastqc.html {output.html}
        mv {params.path}/{wildcards.sample}_{wildcards.lane}_{wildcards.R}_001_fastqc.zip {output.zip}
        """