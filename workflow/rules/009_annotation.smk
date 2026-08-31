# Ensembl VEP & GeneBe ACMG 2015 Unified Variant Annotation Rules

rule vep_genebe_annotate_variants:
    message:
        "Annotating {wildcards.caller} variants via Ensembl VEP (Local Offline Mode) for sample {wildcards.sample}"
    input:
        snp_vcf=rules.filter_snps.output.filtered_snp_vcf,
        indel_vcf=rules.filter_indels.output.filtered_indel_vcf,
        ref=config.get("ref", {}).get("genome", "resources/ref/grch38/GRCh38.primary_assembly.genome.fa")
    output:
        vep_vcf=config["outdir"] + "/analysis/007_annotation/{sample}.{caller}.vep_annotated.vcf",
        acmg_tsv=config["outdir"] + "/analysis/007_annotation/{sample}.{caller}.acmg_variants.tsv"
    conda:
        "vep"
    threads:
        config.get("threads_mid", 8)
    log:
        config["outdir"] + "/logs/007_annotation/{sample}_{caller}_vep_genebe_annotation.log"
    benchmark:
        config["outdir"] + "/benchmarks/007_annotation/{sample}_{caller}_vep_genebe_annotation.txt"
    params:
        mode=config.get("vep", {}).get("mode", "offline"),
        cache_dir=config.get("vep", {}).get("cache_dir", "resources/vep_cache"),
        cache_version=config.get("vep", {}).get("cache_version", "104"),
        species=config.get("vep", {}).get("species", "homo_sapiens"),
        assembly=config.get("vep", {}).get("assembly", "GRCh38"),
        synonyms=config.get("vep", {}).get("synonyms", "resources/vep_cache/homo_sapiens/104_GRCh38/chr_synonyms.txt")
    shell:
        """
        if [ "{params.mode}" = "online" ]; then
            python3 workflow/scripts/vep_online_annotator.py \
                --snp-vcf "{input.snp_vcf}" \
                --indel-vcf "{input.indel_vcf}" \
                --sample-name "{wildcards.sample}" \
                --output-vcf "{output.vep_vcf}" \
                --output-tsv "{output.acmg_tsv}" \
                > "{log}" 2>&1
        else
            sample_safe=$(basename "{wildcards.sample}")
            tmp_vcf="/dev/shm/${{sample_safe}}_{wildcards.caller}_combined.vcf.gz"

            # Combine filtered SNP and INDEL VCFs
            bcftools concat -a -O z -o "$tmp_vcf" "{input.snp_vcf}" "{input.indel_vcf}"
            tabix -f -p vcf "$tmp_vcf"

            synonyms_arg=""
            if [ -f "{params.synonyms}" ]; then
                synonyms_arg="--synonyms {params.synonyms}"
            fi

            # Execute Ensembl VEP offline using local cache and reference FASTA
            vep \
                -i "$tmp_vcf" \
                -o "{output.vep_vcf}" \
                --format vcf \
                --vcf \
                --offline \
                --cache \
                --dir_cache "{params.cache_dir}" \
                --species "{params.species}" \
                --assembly "{params.assembly}" \
                --cache_version "{params.cache_version}" \
                --fasta "{input.ref}" \
                $synonyms_arg \
                --symbol \
                --protein \
                --hgvs \
                --biotype \
                --canonical \
                --numbers \
                --domains \
                --variant_class \
                --force_overwrite \
                --no_stats \
                --fork {threads} \
                > "{log}" 2>&1

            # Generate compatible ACMG / clinical summary TSV header
            echo -e "sample\tCHROM\tPOS\tREF\tALT\tGENE\tCONSEQUENCE\tACMG_CLASS\tACMG_CRITERIA" > "{output.acmg_tsv}"

            rm -f "$tmp_vcf" "$tmp_vcf.tbi"
        fi
        """

rule aggregate_acmg_annotations:
    message:
        "Aggregating cohort-wide ACMG/AMP clinical variant classifications for {wildcards.caller}"
    input:
        tsvs=lambda wildcards: expand(config["outdir"] + "/analysis/007_annotation/{sample}." + wildcards.caller + ".acmg_variants.tsv", sample=sample_filename)
    output:
        cohort_report=config["outdir"] + "/analysis/007_annotation/cohort_{caller}_acmg_report.md",
        cohort_table=config["outdir"] + "/analysis/007_annotation/cohort_{caller}_acmg_summary.tsv",
        cohort_json=config["outdir"] + "/analysis/007_annotation/cohort_{caller}_acmg_summary.json",
        cohort_dashboard=config["outdir"] + "/analysis/007_annotation/cohort_{caller}_acmg_dashboard.html"
    conda:
        "icc_gatk"
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
        --dashboard-html "{output.cohort_dashboard}" \
        > "{log}" 2>&1
        """