# A rule to generate summary reports

rule generate_summary:
    message:
        "Generating summary report for sample {wildcards.sample}_{lane}"
    input:
        annotated_vcf=lambda wildcards: config["outdir"] + "/analysis/007_annotation/{sample}/{sample}_{lane}/{sample}_{lane}.annotated.vcf".format(sample=wildcards.sample, lane=wildcards.lane)
    output:
        summary=lambda wildcards: config["outdir"] + "/analysis/008_summary/{sample}/{sample}_{lane}/{sample}_{lane}_summary.txt".format(sample=wildcards.sample, lane=wildcards.lane)
    conda:
        "envs/008_summary.yml"
    threads:
        config["threads"]
    log:
        lambda wildcards: config["outdir"] + "/logs/008_summary/{sample}/{sample}_{lane}_summary.log".format(sample=wildcards.sample, lane=wildcards.lane)
    benchmark:
        lambda wildcards: config["outdir"] + "/benchmarks/008_summary/{sample}/{sample}_{lane}_summary.txt".format(sample=wildcards.sample, lane=wildcards.lane)
    shell:
        """
        python generate_summary.py \
        --input {input.annotated_vcf} \
        --output {output.summary} \
        > {log} 2>&1
        """
