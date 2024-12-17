# A rule to trim and remove adapters from the raw data using prinseq

rule trimming:
    message: 
        "Trimming and removing adapters from sample {wildcards.sample}_{lane}"
    conda: 
        "icc_02_trimming"
    input:
        fw=lambda wildcards: config["inputdir"] + "/samples/{sample}_{lane}_R1_001.fastq.gz".format(sample=wildcards.sample, lane=wildcards.lane),
        rv=lambda wildcards: config["inputdir"] + "/samples/{sample}_{lane}_R2_001.fastq.gz".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        fw=lambda wildcards: config["outdir"] + "/analysis/002_trimming/{sample}/{sample}_{lane}_R1_trimmed.fastq.gz".format(sample=wildcards.sample, lane=wildcards.lane),
        rv=lambda wildcards: config["outdir"] + "/analysis/002_trimming/{sample}/{sample}_{lane}_R2_trimmed.fastq.gz".format(sample=wildcards.sample, lane=wildcards.lane)
    threads: 
        config["threads"]
    log:
        lambda wildcards: config["outdir"] + "/logs/002_trimming/{sample}/{sample}_{lane}.log".format(sample=wildcards.sample, lane=wildcards.lane)
    benchmark:
        lambda wildcards: config["outdir"] + "/benchmarks/002_trimming/{sample}/{sample}_{lane}.txt".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        prinseq-lite.pl \
        -fastq {input.fw} \
        -fastq2 {input.rv} \
        -out_good {output.fw} \
        -out_bad null \
        -trim_qual_right 20 \
        -trim_qual_left 20 \
        -trim_qual_window 5 \
        -min_len 35 \
        > {log} 2>&1
        """