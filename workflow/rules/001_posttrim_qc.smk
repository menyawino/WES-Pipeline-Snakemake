# A rule to perform post-trimming quality control using FastQC

rule posttrim_fastqc:
    message:
        "Running post-trimming FastQC for sample {wildcards.sample}_{lane}"
    conda:
        "fastqc_env"
    input:
        fq=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_{R}.fastq"
    output:
        zip=config["outdir"] + "/analysis/001_QC/posttrim/{sample}_{lane}_{R}_posttrim_fastqc.zip",
        html=config["outdir"] + "/analysis/001_QC/posttrim/{sample}_{lane}_{R}_posttrim_fastqc.html"
    threads:
        config["threads"]
    params:
        path=config["outdir"] + "/analysis/001_QC/posttrim/{sample}"
    log:
        config["outdir"] + "/logs/001_QC/posttrim/{sample}_{lane}_{R}_posttrim_fastqc.log"
    benchmark:
        config["outdir"] + "/benchmarks/001_QC/posttrim/{sample}_{lane}_{R}_posttrim_fastqc.txt"
    shell:
        """
        # Generate parent directory path
        parent_path=$(dirname {params.path})
        
        # Create the output directory
        mkdir -p "$parent_path"
        
        fastqc {input.fq} \
        -t {threads} \
        -o "$parent_path" \
        &> {log}

        mv {params.path}_{wildcards.lane}_{wildcards.R}_fastqc.zip {output.zip}
        mv {params.path}_{wildcards.lane}_{wildcards.R}_fastqc.html {output.html}
        """