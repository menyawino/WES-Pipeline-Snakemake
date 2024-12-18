# A rule to perform post-trimming quality control using FastQC

rule posttrim_fastqc:
    message:
        "Running post-trimming FastQC for sample {wildcards.sample}_{lane}"
    input:
        fq1=config["outdir"] + "/analysis/002_trimming/{sample}/{sample}_{lane}_R1_trimmed.fastq.gz",
        fq2=config["outdir"] + "/analysis/002_trimming/{sample}/{sample}_{lane}_R2_trimmed.fastq.gz"
    output:
        html1=config["outdir"] + "/qc/{sample}_{lane}_R1_trimmed_fastqc.html",
        html2=config["outdir"] + "/qc/{sample}_{lane}_R2_trimmed_fastqc.html"
    conda:
        "fastqc_env"
    threads:
        config["threads"]
    log:
        config["outdir"] + "/logs/003_posttrim_qc/{sample}_{lane}_posttrim_fastqc.log"
    shell:
        """
        fastqc -t {threads} {input.fq1} {input.fq2} -o {config[outdir]}/qc > {log} 2>&1
        """