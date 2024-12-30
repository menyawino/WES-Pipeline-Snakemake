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