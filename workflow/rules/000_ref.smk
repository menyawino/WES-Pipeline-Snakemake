# Rule to download and index GRCh38 reference genome sequence if missing

ref_fasta = config.get("reference_genome", "resources/ref/grch38/GRCh38.primary_assembly.genome.fa")
ref_fai = ref_fasta + ".fai"

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
