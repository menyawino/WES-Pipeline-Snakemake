# Ensembl VEP & GeneBe ACMG 2015 Unified Variant Annotation Rules

rule vep_genebe_annotate_variants:
    message:
        "Annotating {wildcards.caller} variants and classifying ACMG guidelines via Ensembl VEP & GeneBe plugin for sample {wildcards.sample}"
    input:
        snp_vcf=rules.filter_snps.output.filtered_snp_vcf,
        indel_vcf=rules.filter_indels.output.filtered_indel_vcf
    output:
        vep_vcf=config["outdir"] + "/analysis/007_annotation/{sample}.{caller}.vep_annotated.vcf",
        acmg_tsv=config["outdir"] + "/analysis/007_annotation/{sample}.{caller}.acmg_variants.tsv"
    conda:
        "../envs/009_annotation.yml"
    threads:
        config.get("threads_mid", 8)
    log:
        config["outdir"] + "/logs/007_annotation/{sample}_{caller}_vep_genebe_annotation.log"
    benchmark:
        config["outdir"] + "/benchmarks/007_annotation/{sample}_{caller}_vep_genebe_annotation.txt"
    run:
        if config["vep"].get("mode", "online") == "online":
            shell(
                """
                python3 workflow/scripts/vep_online_annotator.py \
                --snp-vcf "{input.snp_vcf}" \
                --indel-vcf "{input.indel_vcf}" \
                --sample-name "{wildcards.sample}" \
                --output-vcf "{output.vep_vcf}" \
                --output-tsv "{output.acmg_tsv}" \
                > "{log}" 2>&1
                """
            )
        else:
            shell(
                """
                mkdir -p resources/vep_cache
                sample_safe=$(basename "{wildcards.sample}")
                tmp_vcf="/dev/shm/${sample_safe}_{wildcards.caller}_combined.vcf.gz"

                # First, combine the snp and indel VCFs
                bcftools concat -a "{input.snp_vcf}" "{input.indel_vcf}" -O z -o "$tmp_vcf"
                
                # Run offline VEP
                vep -i "$tmp_vcf" \
                -o "{output.vep_vcf}" \
                --offline --cache --dir_cache resources/vep_cache \
                --species homo_sapiens --assembly GRCh38 \
                --vcf --force_overwrite --fork {threads} \
                > "{log}" 2>&1
                
                # Touch empty TSV since offline VEP doesn't run GeneBe ACMG
                touch "{output.acmg_tsv}"
                rm -f "$tmp_vcf"
                """
            )

rule aggregate_acmg_annotations:
    message:
        "Aggregating cohort-wide ACMG/AMP clinical variant classifications for {wildcards.caller}"
    input:
        tsvs=lambda wildcards: expand(config["outdir"] + "/analysis/007_annotation/{sample}." + wildcards.caller + ".acmg_variants.tsv", sample=sample_filename)
    output:
        cohort_report=config["outdir"] + "/analysis/007_annotation/cohort_{caller}_acmg_report.md",
        cohort_table=config["outdir"] + "/analysis/007_annotation/cohort_{caller}_acmg_summary.tsv",
        cohort_json=config["outdir"] + "/analysis/007_annotation/cohort_{caller}_acmg_summary.json"
    conda:
        "../envs/009_annotation.yml"
    log:
        config["outdir"] + "/logs/007_annotation/cohort_{caller}_acmg_summary.log"
    benchmark:
        config["outdir"] + "/benchmarks/007_annotation/cohort_{caller}_acmg_summary.txt"
    shell:
        """
        python3 workflow/scripts/aggregate_acmg.py \
        --inputs {input.tsvs} \
        --report-md "{output.cohort_report}" \
        --summary-tsv "{output.cohort_table}" \
        --summary-json "{output.cohort_json}" \
        > "{log}" 2>&1
        """