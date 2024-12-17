# A rule to perform post-trimming quality control using FastQC

rule posttrim_fastqc:
    message:
        "Running post-trimming FastQC for sample {wildcards.sample}_{lane}"
    input:
        fq1=lambda wildcards: config["outdir"] + "/analysis/002_trimming/{sample}/{sample}_{lane}_R1_trimmed.fastq.gz".format(sample=wildcards.sample, lane=wildcards.lane),
        fq2=lambda wildcards: config["outdir"] + "/analysis/002_trimming/{sample}/{sample}_{lane}_R2_trimmed.fastq.gz".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        html1=lambda wildcards: config["outdir"] + "/qc/{sample}_{lane}_R1_trimmed_fastqc.html".format(sample=wildcards.sample, lane=wildcards.lane),
        html2=lambda wildcards: config["outdir"] + "/qc/{sample}_{lane}_R2_trimmed_fastqc.html".format(sample=wildcards.sample, lane=wildcards.lane)
    conda:
        "fastqc_env"
    threads:
        config["threads"]
    log:
        lambda wildcards: config["outdir"] + "/logs/003_posttrim_qc/{sample}_{lane}_posttrim_fastqc.log".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        fastqc -t {threads} {input.fq1} {input.fq2} -o {config[outdir]}/qc > {log} 2>&1
        """