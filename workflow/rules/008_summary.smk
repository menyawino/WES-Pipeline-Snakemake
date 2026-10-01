# ClawBio-inspired variant summaries built directly from filtered SNP/indel VCFs.

rule summarize_variants:
    message:
        "Generating {wildcards.caller} variant summary for sample {wildcards.sample}"
    input:
        snp_vcf=rules.filter_snps.output.filtered_snp_vcf,
        indel_vcf=rules.filter_indels.output.filtered_indel_vcf,
        annotated_vcf=rules.vep_genebe_annotate_variants.output.vep_vcf
    output:
        report_md=config["outdir"] + "/analysis/008_summary/{sample}/{caller}_variant_summary.md",
        summary_tsv=config["outdir"] + "/analysis/008_summary/{sample}/{caller}_variant_summary.tsv",
        summary_json=config["outdir"] + "/analysis/008_summary/{sample}/{caller}_variant_summary.json"
    conda:
        "../envs/005_gatk_genomics.yml"
    threads:
        2
    resources:
        mem_mb=config.get("mem_low", 4096)
    params:
        top_variants=config.get("variant_summary", {}).get("top_variants", 10),
        top_chromosomes=config.get("variant_summary", {}).get("top_chromosomes", 10)
    log:
        config["outdir"] + "/logs/008_summary/{sample}_{caller}_variant_summary.log"
    benchmark:
        config["outdir"] + "/benchmarks/008_summary/{sample}_{caller}_variant_summary.txt"
    shell:
        """
        python workflow/scripts/variant_summary.py sample \
        --sample-name "{wildcards.sample}_{wildcards.caller}" \
        --snp-vcf "{input.snp_vcf}" \
        --indel-vcf "{input.indel_vcf}" \
        --annotated-vcf "{input.annotated_vcf}" \
        --report-md "{output.report_md}" \
        --summary-tsv "{output.summary_tsv}" \
        --summary-json "{output.summary_json}" \
        --top-variants {params.top_variants} \
        --top-chromosomes {params.top_chromosomes} \
        > "{log}" 2>&1
        """

rule aggregate_variant_summaries:
    message:
        "Aggregating cohort variant summaries for {wildcards.caller}"
    input:
        reports=lambda wildcards: expand(config["outdir"] + "/analysis/008_summary/{sample}/" + wildcards.caller + "_variant_summary.json", sample=sample_filename)
    output:
        cohort_report=config["outdir"] + "/analysis/008_summary/cohort_{caller}_variant_report.md",
        cohort_table=config["outdir"] + "/analysis/008_summary/cohort_{caller}_variant_summary.tsv",
        cohort_json=config["outdir"] + "/analysis/008_summary/cohort_{caller}_variant_summary.json"
    conda:
        "../envs/005_gatk_genomics.yml"
    threads:
        1
    resources:
        mem_mb=config.get("mem_low", 4096)
    log:
        config["outdir"] + "/logs/008_summary/cohort_{caller}_variant_summary.log"
    benchmark:
        config["outdir"] + "/benchmarks/008_summary/cohort_{caller}_variant_summary.txt"
    shell:
        """
        python workflow/scripts/variant_summary.py cohort \
        --inputs {input.reports} \
        --report-md "{output.cohort_report}" \
        --summary-tsv "{output.cohort_table}" \
        --summary-json "{output.cohort_json}" \
        > "{log}" 2>&1
        """

rule compare_callers:
    message:
        "Evaluating concordance and discordance between GATK and DeepVariant calls"
    input:
        gatk_snps=expand(config["outdir"] + "/analysis/006_variant_filtering/{sample}.gatk.filtered.snp.vcf", sample=sample_filename),
        dv_snps=expand(config["outdir"] + "/analysis/006_variant_filtering/{sample}.deepvariant.filtered.snp.vcf", sample=sample_filename)
    output:
        report=config["outdir"] + "/analysis/008_summary/caller_concordance_report.md",
        table=config["outdir"] + "/analysis/008_summary/caller_concordance_matrix.tsv"
    conda:
        "../envs/005_gatk_genomics.yml"
    threads:
        2
    resources:
        mem_mb=config.get("mem_low", 4096)
    params:
        outdir=config["outdir"],
        samples=sample_filename
    log:
        config["outdir"] + "/logs/008_summary/caller_concordance.log"
    benchmark:
        config["outdir"] + "/benchmarks/008_summary/caller_concordance.txt"
    shell:
        """
        python3 workflow/scripts/compare_callers.py \
        --outdir "{params.outdir}" \
        --samples {params.samples} \
        --report-md "{output.report}" \
        --summary-tsv "{output.table}" \
        > "{log}" 2>&1
        """

rule coverage_summary_report:
    message:
        "Generating legacy 46-column cohort coverage summary tables (Target, ProteinCodingTarget, CanonTranCodingTarget)"
    input:
        qc_metrics=expand(rules.qc_report.output.qc_metrics, sample=sample_filename),
        orig_flagstats=expand(rules.flagstat_original.output.flagstat_original, sample=sample_filename),
        q8_counts=expand(rules.count_mapped_q8_original.output.q8_count, sample=sample_filename),
        target_flagstats=expand(rules.flagstat_target.output.flagstat_target, sample=sample_filename),
        prot_coding_flagstats=expand(rules.flagstat_prot_coding.output.flagstat_prot_coding, sample=sample_filename),
        canon_tran_flagstats=expand(rules.flagstat_canon_tran.output.flagstat_canon_tran, sample=sample_filename),
        snps=expand(rules.filter_snps.output.filtered_snp_vcf, sample=sample_filename, caller=["gatk"]),
        indels=expand(rules.filter_indels.output.filtered_indel_vcf, sample=sample_filename, caller=["gatk"])
    output:
        target_txt=config["outdir"] + "/results/Coverage_Report/SummaryOutput_Target_" + os.path.basename(config["outdir"].rstrip("/")) + ".txt",
        target_tsv=config["outdir"] + "/results/Coverage_Report/SummaryOutput_Target_" + os.path.basename(config["outdir"].rstrip("/")) + ".tsv",
        prot_coding_txt=config["outdir"] + "/results/Coverage_Report/SummaryOutput_ProteinCodingTarget_" + os.path.basename(config["outdir"].rstrip("/")) + ".txt",
        prot_coding_tsv=config["outdir"] + "/results/Coverage_Report/SummaryOutput_ProteinCodingTarget_" + os.path.basename(config["outdir"].rstrip("/")) + ".tsv",
        canon_tran_txt=config["outdir"] + "/results/Coverage_Report/SummaryOutput_CanonTranCodingTarget_" + os.path.basename(config["outdir"].rstrip("/")) + ".txt",
        canon_tran_tsv=config["outdir"] + "/results/Coverage_Report/SummaryOutput_CanonTranCodingTarget_" + os.path.basename(config["outdir"].rstrip("/")) + ".tsv",
        report_md=config["outdir"] + "/results/Coverage_Report/Coverage_Summary_Report.md"
    conda:
        "../envs/005_gatk_genomics.yml"
    threads:
        2
    resources:
        mem_mb=config.get("mem_low", 4096)
    params:
        outdir=config["outdir"],
        samples=sample_filename,
        target_bed=config["icc_panel"],
        prot_coding_bed=config["cds_panel"],
        canon_tran_bed=config["canontran_panel"],
        run_id=os.path.basename(config["outdir"].rstrip("/")),
        genome_size=config.get("coverage_summary", {}).get("genome_size", 3095693981),
        min_snps=config.get("coverage_summary", {}).get("min_snps", 200),
        min_indels=config.get("coverage_summary", {}).get("min_indels", 20),
        min_titv=config.get("coverage_summary", {}).get("min_titv", 2.0)
    log:
        config["outdir"] + "/logs/008_summary/coverage_summary_report.log"
    benchmark:
        config["outdir"] + "/benchmarks/008_summary/coverage_summary_report.txt"
    shell:
        """
        python3 workflow/scripts/coverage_summary.py \
        --outdir "{params.outdir}" \
        --samples {params.samples} \
        --target-bed "{params.target_bed}" \
        --prot-coding-bed "{params.prot_coding_bed}" \
        --canon-tran-bed "{params.canon_tran_bed}" \
        --run-id "{params.run_id}" \
        --genome-size {params.genome_size} \
        --min-snps {params.min_snps} \
        --min-indels {params.min_indels} \
        --min-titv {params.min_titv} \
        > "{log}" 2>&1
        """
