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
        gvcf_chunk=temp("/dev/shm/wes_pipeline_{sample}_{chunk}.g.vcf.gz"),
        gvcf_tbi=temp("/dev/shm/wes_pipeline_{sample}_{chunk}.g.vcf.gz.tbi")
    container:
        "docker://broadinstitute/gatk:4.4.0.0"
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
        gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads={threads} -XX:ConcGCThreads={threads} -Xmx{resources.mem_mb}m" HaplotypeCaller \
        -R "{params.ref}" \
        -I "{input.bam}" \
        -O "{output.gvcf_chunk}" \
        -ERC GVCF \
        -L "{input.interval}" \
        --native-pair-hmm-threads {threads} \
        --smith-waterman FASTEST_AVAILABLE \
        --tmp-dir "/dev/shm" \
        &> "{log}"
        """

rule gather_gvcfs:
    message:
        "Gathering interval GVCFs for sample {wildcards.sample}"
    input:
        gvcfs=expand("/dev/shm/wes_pipeline_{{sample}}_{chunk}.g.vcf.gz", chunk=scatter_chunks),
        tbis=expand("/dev/shm/wes_pipeline_{{sample}}_{chunk}.g.vcf.gz.tbi", chunk=scatter_chunks)
    output:
        gvcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.gatk.g.vcf.gz",
        gvcf_tbi=config["outdir"] + "/analysis/005_variant_calling/{sample}.gatk.g.vcf.gz.tbi"
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
        inputs_args=""
        for f in {input.gvcfs}; do
            inputs_args="$inputs_args -I $f"
        done
        gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads={threads} -XX:ConcGCThreads={threads} -Xmx{resources.mem_mb}m" GatherVcfs \
        $inputs_args \
        -O "{output.gvcf}" \
        --TMP_DIR "/dev/shm" \
        > "{log}" 2>&1

        tabix -f -p vcf "{output.gvcf}" >> "{log}" 2>&1
        """

rule genotype_gvcfs:
    message:
        "Genotyping GATK GVCFs for sample {wildcards.sample}"
    input:
        gvcf=rules.gather_gvcfs.output.gvcf,
        gvcf_tbi=rules.gather_gvcfs.output.gvcf_tbi
    output:
        vcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.gatk.vcf.gz",
        tbi=config["outdir"] + "/analysis/005_variant_calling/{sample}.gatk.vcf.gz.tbi"
    container:
        "docker://broadinstitute/gatk:4.4.0.0"
    threads:
        config["threads_mid"]
    resources:
        mem_mb=config.get("mem_mid", 16384),
        tmpdir=config.get("tmpdir", "/tmp")
    params:
        ref="/dev/shm/wes_ref_grch38/GRCh38.primary_assembly.genome.fa",
        target=config["icc_panel"],
        dbsnp=config["dbsnp"]
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_genotype_gvcfs.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_genotype_gvcfs.txt"
    shell:
        """
        raw_vcf="{output.vcf}.tmp.vcf"
        gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads={threads} -XX:ConcGCThreads={threads} -Xmx{resources.mem_mb}m" GenotypeGVCFs \
        -R "{params.ref}" \
        -V "{input.gvcf}" \
        -O "$raw_vcf" \
        -G StandardAnnotation \
        --intervals "{params.target}" \
        --interval-padding 100 \
        --create-output-variant-index false \
        --dbsnp "{params.dbsnp}" \
        --tmp-dir "/dev/shm" \
        &> "{log}"

        bgzip -c "$raw_vcf" > "{output.vcf}" 2>> "{log}"
        tabix -f -p vcf "{output.vcf}" >> "{log}" 2>&1
        rm -f "$raw_vcf"
        """

rule deepvariant_call:
    message:
        "Calling variants with Google DeepVariant (WES model) for sample {wildcards.sample}"
    input:
        bam=rules.filter_bam_target.output.bam_target,
        bai=rules.filter_bam_target.output.bai_target,
        target=config["icc_panel"],
        ref_staged=rules.stage_ref_shm.output.staged_done
    output:
        vcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.deepvariant.vcf.gz",
        tbi=config["outdir"] + "/analysis/005_variant_calling/{sample}.deepvariant.vcf.gz.tbi",
        gvcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.deepvariant.g.vcf.gz"
    container:
        "docker://google/deepvariant:1.6.1"
    threads:
        config.get("threads_mid", 8)
    resources:
        mem_mb=config.get("mem_mid", 16384),
        tmpdir=config.get("tmpdir", "/tmp")
    params:
        ref="/dev/shm/wes_ref_grch38/GRCh38.primary_assembly.genome.fa",
        model_type=config.get("deepvariant", {}).get("model_type", "WES")
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_deepvariant.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_deepvariant.txt"
    shell:
        """
        mkdir -p "{resources.tmpdir}"
        dv_tmp="{resources.tmpdir}/dv_{wildcards.sample}"
        mkdir -p "$dv_tmp"

        if command -v /opt/deepvariant/bin/run_deepvariant >/dev/null 2>&1; then
            /opt/deepvariant/bin/run_deepvariant \
                --model_type="{params.model_type}" \
                --ref="{params.ref}" \
                --reads="{input.bam}" \
                --regions="{input.target}" \
                --output_vcf="{output.vcf}" \
                --output_gvcf="{output.gvcf}" \
                --num_shards={threads} \
                --intermediate_results_dir="$dv_tmp" \
                &> "{log}"
        else
            docker run --rm \
                -v /dev/shm:/dev/shm \
                -v /mnt/bucket:/mnt/bucket \
                -v /mnt/qnap-public:/mnt/qnap-public \
                google/deepvariant:1.6.1 \
                /opt/deepvariant/bin/run_deepvariant \
                --model_type="{params.model_type}" \
                --ref="{params.ref}" \
                --reads="{input.bam}" \
                --regions="{input.target}" \
                --output_vcf="{output.vcf}" \
                --output_gvcf="{output.gvcf}" \
                --num_shards={threads} \
                --intermediate_results_dir="$dv_tmp" \
                &> "{log}"
        fi

        tabix -f -p vcf "{output.vcf}" >> "{log}" 2>&1
        rm -rf "$dv_tmp"
        """

def get_raw_caller_vcf(wildcards):
    if wildcards.caller == "gatk":
        return rules.genotype_gvcfs.output.vcf
    elif wildcards.caller == "deepvariant":
        return rules.deepvariant_call.output.vcf
    else:
        return config["outdir"] + f"/analysis/005_variant_calling/{wildcards.sample}.{wildcards.caller}.vcf.gz"

rule split_vcfs:
    message:
        "Splitting {wildcards.caller} VCF into SNPs and Indels for sample {wildcards.sample}"
    input:
        vcf=get_raw_caller_vcf
    output:
        snp_vcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.{caller}.snp.vcf.gz",
        snp_tbi=config["outdir"] + "/analysis/005_variant_calling/{sample}.{caller}.snp.vcf.gz.tbi",
        indel_vcf=config["outdir"] + "/analysis/005_variant_calling/{sample}.{caller}.indel.vcf.gz",
        indel_tbi=config["outdir"] + "/analysis/005_variant_calling/{sample}.{caller}.indel.vcf.gz.tbi"
    conda:
        "icc_gatk"
    threads:
        config.get("threads_low", 4)
    resources:
        mem_mb=config.get("mem_low", 4096),
        tmpdir=config.get("tmpdir", "/tmp")
    log:
        config["outdir"] + "/logs/005_variant_calling/{sample}_{caller}_split_vcfs.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_variant_calling/{sample}_{caller}_split_vcfs.txt"
    shell:
        """
        # Split into SNPs
        bcftools view -v snps -O z -o "{output.snp_vcf}" "{input.vcf}" > "{log}" 2>&1
        tabix -f -p vcf "{output.snp_vcf}" >> "{log}" 2>&1

        # Split into Indels
        bcftools view -v indels,other -O z -o "{output.indel_vcf}" "{input.vcf}" >> "{log}" 2>&1
        tabix -f -p vcf "{output.indel_vcf}" >> "{log}" 2>&1
        """
