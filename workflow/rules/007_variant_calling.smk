num_chunks = int(config.get("intervals_scatter_count", 20))
scatter_chunks = [f"{i:02d}" for i in range(num_chunks)]

rule split_target_intervals:
    message:
        "Splitting target intervals into balanced chunks for parallel variant calling"
    input:
        target=config["icc_panel"],
        ref_fai=config["reference_genome"] + ".fai"
    output:
        chunks=expand(config["outdir"] + "/analysis/000_ref/intervals/chunk_{chunk}.bed", chunk=scatter_chunks)
    conda:
        "icc_gatk"
    log:
        config["outdir"] + "/logs/000_ref/split_target_intervals.log"
    shell:
        """
        python3 workflow/scripts/split_bed.py "{input.target}" "{input.ref_fai}" {output.chunks} > "{log}" 2>&1
        """

rule haplotypecaller_chunk:
    message:
        "Calling variants with HaplotypeCaller for sample {wildcards.sample} (chunk {wildcards.chunk})"
    input:
        bam=rules.filter_bam_target.output.bam_target,
        bai=rules.filter_bam_target.output.bai_target,
        interval=config["outdir"] + "/analysis/000_ref/intervals/chunk_{chunk}.bed"
    output:
        gvcf_chunk=temp(config["outdir"] + "/analysis/005_variant_calling/{sample}/chunks/{chunk}.g.vcf.gz"),
        gvcf_tbi=temp(config["outdir"] + "/analysis/005_variant_calling/{sample}/chunks/{chunk}.g.vcf.gz.tbi")
    conda:
        "icc_gatk"
    threads:
        config.get("threads_low", 4)
    resources:
        mem_mb=config.get("mem_low", 8192),
        tmpdir=config.get("tmpdir", "/tmp")
    params:
        ref=config["reference_genome"]
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_haplotypecaller_chunk_{chunk}.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_haplotypecaller_chunk_{chunk}.txt"
    shell:
        """
        mkdir -p "{resources.tmpdir}"
        gatk --java-options "-XX:+UseParallelGC -Xmx{resources.mem_mb}m" HaplotypeCaller \
        -R "{params.ref}" \
        -I "{input.bam}" \
        -O "{output.gvcf_chunk}" \
        -ERC GVCF \
        -L "{input.interval}" \
        --interval-padding 100 \
        --native-pair-hmm-threads {threads} \
        --tmp-dir "{resources.tmpdir}" \
        &> "{log}"
        """

rule gather_gvcfs:
    message:
        "Gathering interval GVCFs for sample {wildcards.sample}"
    input:
        gvcfs=expand(config["outdir"] + "/analysis/005_variant_calling/{{sample}}/chunks/{chunk}.g.vcf.gz", chunk=scatter_chunks)
    output:
        gvcf=temp(config["outdir"] + "/analysis/005_variant_calling/{sample}.haplotypecaller.g.vcf.gz"),
        gvcf_tbi=temp(config["outdir"] + "/analysis/005_variant_calling/{sample}.haplotypecaller.g.vcf.gz.tbi")
    conda:
        "icc_gatk"
    threads:
        config.get("threads_mid", 8)
    resources:
        mem_mb=config.get("mem_mid", 16384),
        tmpdir=config.get("tmpdir", "/tmp")
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_gather_gvcfs.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_gather_gvcfs.txt"
    shell:
        """
        mkdir -p "{resources.tmpdir}"
        inputs_args=""
        for f in {input.gvcfs}; do
            inputs_args="$inputs_args -I $f"
        done
        gatk --java-options "-XX:+UseParallelGC -Xmx{resources.mem_mb}m" GatherVcfs \
        $inputs_args \
        -O "{output.gvcf}" \
        --TMP_DIR "{resources.tmpdir}" \
        > "{log}" 2>&1

        tabix -p vcf "{output.gvcf}" >> "{log}" 2>&1
        """

rule genotype_gvcfs:
    message:
        "Genotyping GVCFs for sample {wildcards.sample}"
    input:
        gvcf=rules.gather_gvcfs.output.gvcf,
        gvcf_tbi=rules.gather_gvcfs.output.gvcf_tbi
    output:
        vcf=temp(config["outdir"] + "/analysis/005_variant_calling/{sample}.genotyped.vcf")
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    resources:
        mem_mb=config.get("mem_mid", 16384),
        tmpdir=config.get("tmpdir", "/tmp")
    params:
        ref=config["reference_genome"],
        target=config["icc_panel"],
        dbsnp=config["dbsnp"]
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_genotype_gvcfs.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_genotype_gvcfs.txt"
    shell:
        """
        mkdir -p "{resources.tmpdir}"
        gatk --java-options "-XX:+UseParallelGC -Xmx{resources.mem_mb}m" GenotypeGVCFs \
        -R "{params.ref}" \
        -V "{input.gvcf}" \
        -O "{output.vcf}" \
        -G StandardAnnotation \
        --intervals "{params.target}" \
        --interval-padding 100 \
        --create-output-variant-index true \
        --dbsnp "{params.dbsnp}" \
        --tmp-dir "{resources.tmpdir}" \
        &> "{log}"
        """

rule split_vcfs:
    message:
        "Splitting VCFs into SNPs and Indels for sample {wildcards.sample}"
    input:
        vcf=rules.genotype_gvcfs.output.vcf
    output:
        snp_vcf=temp(config["outdir"] + "/analysis/005_variant_calling/{sample}.genotyped.snp.vcf"),
        indel_vcf=temp(config["outdir"] + "/analysis/005_variant_calling/{sample}.genotyped.indel.vcf")
    conda:
        "icc_gatk"
    threads:
        config["threads_mid"]
    resources:
        mem_mb=config.get("mem_mid", 16384),
        tmpdir=config.get("tmpdir", "/tmp")
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_split_vcfs.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_split_vcfs.txt"
    shell:
        """
        mkdir -p "{resources.tmpdir}"
        gatk --java-options "-XX:+UseParallelGC -Xmx{resources.mem_mb}m" SplitVcfs \
        -I "{input.vcf}" \
        --SNP_OUTPUT "{output.snp_vcf}" \
        --INDEL_OUTPUT "{output.indel_vcf}" \
        --STRICT false \
        --TMP_DIR "{resources.tmpdir}" \
        &> "{log}"
        """
