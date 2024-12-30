# A rule to trim and remove adapters from the raw data using prinseq

rule trimming:
    message: 
        "Trimming and removing adapters from sample {wildcards.sample}_{lane}"
    conda: 
        "icc_02_trimming"
    input:
        fw=config["inputdir"] + "/{sample}_{lane}_R1_001.fastq.gz",
        rv=config["inputdir"] + "/{sample}_{lane}_R2_001.fastq.gz",
    output:
        fw=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_R1.fastq",
        fws=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_R1_singletons.fastq",
        rv=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_R2.fastq",
        rvs=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_R2_singletons.fastq"
    threads:
        config["threads"]
    params:
        path=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}"
    log:
        config["outdir"] + "/logs/002_trimming/{sample}/{sample}_{lane}.log"
    benchmark:
        config["outdir"] + "/benchmarks/002_trimming/{sample}/{sample}_{lane}.txt"
    shell:
        """
        # unzip the input files
        gunzip -c {input.fw} > {input.fw}.fastq
        gunzip -c {input.rv} > {input.rv}.fastq

        prinseq-lite.pl \
        -fastq {input.fw}.fastq \
        -fastq2 {input.rv}.fastq \
        -out_good {params.path} \
        -out_bad null \
        -trim_qual_right 20 \
        -trim_qual_left 20 \
        -trim_qual_window 5 \
        -min_len 35 \
        &> {log}

        # rename output files to suit the rest of the pipeline
        mv {params.path}_1.fastq {output.fw}
        mv {params.path}_1_singletons.fastq {output.fws}
        mv {params.path}_2.fastq {output.rv}
        mv {params.path}_2_singletons.fastq {output.rvs}
        """