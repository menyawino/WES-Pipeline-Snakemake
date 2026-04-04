# A rule to generate summary reports

rule generate_summary:
    message:
        "Generating summary report for sample {wildcards.sample}_{lane}"
    input:
        annotated_vcf=config["outdir"] + "/analysis/007_annotation/{sample}/{sample}_{lane}/{sample}_{lane}.annotated.vcf"
    output:
        summary=config["outdir"] + "/analysis/008_summary/{sample}/{sample}_{lane}/{sample}_{lane}_summary.txt"
    conda:
        "envs/008_summary.yml"
    threads:
        config["threads"]
    log:
        config["outdir"] + "/logs/008_summary/{sample}/{sample}_{lane}_summary.log"
    benchmark:
        config["outdir"] + "/benchmarks/008_summary/{sample}/{sample}_{lane}_summary.txt"
    shell:
        """
        python generate_summary.py \
        --input {input.annotated_vcf} \
        --output {output.summary} \
        > {log} 2>&1
        """
