# A rule to run FastQC on the raw data

rule raw_fastqc:
    message: 
        "Running FastQC on sample {lane}_{wildcards.R}"
    conda: 
        "fastqc_env"
    input:
        config["inputdir"] + "/{sample_filename}_{lane}_{R}_001.fastq.gz"
    output:
        html=config["outdir"] + "/analysis/001_QC/{sample_filename}_{lane}_{R}_fastqc.html",
        zip=config["outdir"] + "/analysis/001_QC/{sample_filename}_{lane}_{R}_fastqc.zip"
    threads: 
        config["np_threads"]
    params: 
        path=config["outdir"] + "/analysis/001_QC/{sample_filename}",
    log:
        config["outdir"] + "/logs/001_QC/{sample_filename}_{lane}_{R}.log"
    benchmark:
        config["outdir"] + "/benchmarks/001_QC/{sample_filename}_{lane}_{R}.txt"
    shell:
        """
        # Generate parent directory path
        parent_path=$(dirname {params.path})

        # Create the output directory
        mkdir -p "$parent_path"

        # Run FastQC
        fastqc {input} \
        -t {threads} \
        -o "$parent_path" \
        > {log} 2>&1

        # Move the FastQC outputs to the desired location
        mv {params.path}_{wildcards.lane}_{wildcards.R}_001_fastqc.html {output.html}
        mv {params.path}_{wildcards.lane}_{wildcards.R}_001_fastqc.zip {output.zip}
        """