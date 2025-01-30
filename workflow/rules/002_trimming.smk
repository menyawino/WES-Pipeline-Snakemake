# A rule to trim and remove adapters from the raw data using prinseq

rule trimming:
    message: 
        "Trimming and removing adapters from sample {wildcards.sample}_{lane}"
    conda: 
        "icc_02_trimming"
    input:
        fq1=config["inputdir"] + "/{sample}_{lane}_R1_001.fastq.gz",
        fq2=config["inputdir"] + "/{sample}_{lane}_R2_001.fastq.gz",
    output:
        fq1=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_R1.fastq",
        fq1s=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_R1_singletons.fastq",
        fq2=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_R2.fastq",
        fq2s=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_R2_singletons.fastq"
    params:
        path=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}"
    log:
        config["outdir"] + "/logs/002_trimming/{sample}/{sample}_{lane}.log"
    benchmark:
        config["outdir"] + "/benchmarks/002_trimming/{sample}/{sample}_{lane}.txt"
    shell:
        """
        # unzip the input files
        gunzip -c {input.fq1} > {input.fq1}.fastq
        gunzip -c {input.fq2} > {input.fq2}.fastq

        prinseq-lite.pl \
        -fastq {input.fq1}.fastq \
        -fastq2 {input.fq2}.fastq \
        -out_good {params.path} \
        -out_bad null \
        -trim_qual_right 20 \
        -trim_qual_left 20 \
        -trim_qual_window 5 \
        -min_len 35 \
        &> {log}

        # rename output files to suit the rest of the pipeline
        mv {params.path}_1.fastq {output.fq1}
        mv {params.path}_1_singletons.fastq {output.fq1s}
        mv {params.path}_2.fastq {output.fq2}
        mv {params.path}_2_singletons.fastq {output.fq2s}
        """

rule trimming_fp:
    message: 
        "Trimming and removing adapters from sample {wildcards.sample}_{lane}"
    conda: 
        "icc_02_trimming"
    input:
        fq1=config["inputdir"] + "/{sample}_{lane}_R1_001.fastq.gz",
        fq2=config["inputdir"] + "/{sample}_{lane}_R2_001.fastq.gz",
    output:
        fq1=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_R1.fastq.gz",
        fq2=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_R2.fastq.gz",
        report=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_report.html",
        json=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_report.json"
    threads:
        config["threads_low"]
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
        -i {input.fq1} \
        -I {input.fq2} \
        -j {output.json} \
        -o {output.fq1} \
        -O {output.fq2} \
        -h {output.report} \
        -w {threads} \
        --length_required {params.min_length} \
        --cut_window_size {params.window_size} \
        &> {log}
        """


        