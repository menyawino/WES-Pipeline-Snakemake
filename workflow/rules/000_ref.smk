# Rule to download and index GRCh38 reference genome sequence if missing

ref_fasta = config.get("reference_genome", "resources/ref/grch38/GRCh38.primary_assembly.genome.fa")
ref_fai = ref_fasta + ".fai"
ref_dict = os.path.splitext(ref_fasta)[0] + ".dict"

rule download_grch38_reference:
    message:
        "Downloading and indexing GRCh38 reference genome FASTA"
    output:
        fasta=ref_fasta,
        fai=ref_fai
    params:
        genome=config.get("Genome", "grch38")
    log:
        config["outdir"] + "/logs/000_ref/download_grch38.log"
    benchmark:
        config["outdir"] + "/benchmarks/000_ref/download_grch38.txt"
    shell:
        """
        python3 workflow/scripts/download_ref.py {output.fasta} &> {log}
        """

rule stage_ref_shm:
    message:
        "Pre-loading reference genome and BWA-MEM2 index into RAM disk (/dev/shm) for high-speed alignment"
    input:
        ref=ref_fasta,
        bwt=ref_fasta + ".bwt.2bit.64",
        b0123=ref_fasta + ".0123",
        amb=ref_fasta + ".amb",
        ann=ref_fasta + ".ann",
        pac=ref_fasta + ".pac",
        fai=ref_fai,
        dict=ref_dict
    output:
        staged_done=touch(temp(config["outdir"] + "/analysis/000_ref/.staged_ref_shm.done"))
    shell:
        """
        mkdir -p /dev/shm/wes_ref_grch38
        for f in {input.ref} {input.bwt} {input.b0123} {input.amb} {input.ann} {input.pac} {input.fai} {input.dict}; do
            target="/dev/shm/wes_ref_grch38/$(basename "$f")"
            if [ ! -f "$target" ] || [ "$f" -nt "$target" ]; then
                cp -u "$f" "$target"
            fi
        done
        """
