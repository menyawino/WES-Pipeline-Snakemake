# A rule to run FastQC on the raw data

rule raw_fastqc:
    message: 
        "Running FastQC on sample {wildcards.sample}_{lane}_{wildcards.R}"
    conda: 
        "fastqc_env"
    input:
        lambda wildcards: config["inputdir"] + "/samples/{sample}_{lane}_{R}_001.fastq.gz".format(sample=wildcards.sample, lane=wildcards.lane, R=wildcards.R)
    output:
        html=lambda wildcards: config["outdir"] + "/analysis/001_QC/{sample}/{sample}_{lane}_{R}_fastqc.html".format(sample=wildcards.sample, lane=wildcards.lane, R=wildcards.R),
        zip=lambda wildcards: config["outdir"] + "/analysis/001_QC/{sample}/{sample}_{lane}_{R}_fastqc.zip".format(sample=wildcards.sample, lane=wildcards.lane, R=wildcards.R)
    threads: 
        config["np_threads"]
    params: 
        path=lambda wildcards: config["outdir"] + "/analysis/001_QC/{}".format(wildcards.sample)
    log:
        lambda wildcards: config["outdir"] + "/logs/001_QC/{sample}/{sample}_{lane}_{R}.log".format(sample=wildcards.sample, lane=wildcards.lane, R=wildcards.R)
    benchmark:
        lambda wildcards: config["outdir"] + "/benchmarks/001_QC/{sample}/{sample}_{lane}_{R}.txt".format(sample=wildcards.sample, lane=wildcards.lane, R=wildcards.R)
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