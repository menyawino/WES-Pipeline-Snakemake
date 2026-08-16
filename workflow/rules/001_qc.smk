# A rule to run FastQC on the raw data

rule raw_fastqc:
    message: 
        "Running FastQC on sample {lane}_{wildcards.R}"
    conda: 
        "fastqc_env"
    input:
        get_raw_fastqc_input
    output:
        html=config["outdir"] + "/analysis/001_QC/pretrim/{sample_filename}_{lane}_{R}_pretrim_fastqc.html",
        zip=config["outdir"] + "/analysis/001_QC/pretrim/{sample_filename}_{lane}_{R}_pretrim_fastqc.zip"
    threads: 
        config["threads_mid"]
    params: 
        path=config["outdir"] + "/analysis/001_QC/pretrim/{sample_filename}",
    log:
        config["outdir"] + "/logs/001_QC/pretrim/{sample_filename}_{lane}_{R}.log"
    benchmark:
        config["outdir"] + "/benchmarks/001_QC/pretrim/{sample_filename}_{lane}_{R}.txt"
    shell:
        """
        # Generate parent directory path
        parent_path=$(dirname "{params.path}")

        # Create the output directory
        mkdir -p "$parent_path"

        fastqc "{input}" \
        -t {threads} \
        -o "$parent_path" \
        &> "{log}"

        # FastQC names output based on input filename
        in_base=$(basename "{input}" .fastq.gz)
        in_base=${{in_base%.fq.gz}}

        # Move the FastQC outputs to the desired location
        mv "$parent_path/${{in_base}}_fastqc.html" "{output.html}"
        mv "$parent_path/${{in_base}}_fastqc.zip" "{output.zip}"
        """

 
