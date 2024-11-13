# Convert metadata table to JSON
# Then validate md5 checksums of fasta files


rule metadat2json:
    output:
        metadata_json="validation/metadata.json",
    params:
        metadata=config["metadata"],
    message:
        "[Input validation] Exporting metadata to JSON"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/metadat2json.log",
    script:
        "../scripts/metadata2json.py"


rule validate_checksum:
    input:
        metadata="validation/metadata.json",
    output:
        json="validation/checksum_status.json",
        tsv="validation/checksum_status.tsv",
    params:
        fasta_dir=config["fasta_dir"],
    message:
        "[Input validation] validating md5 checksums"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/validate_checksum.log",
    script:
        "../scripts/validate_checksums.py"


checkpoint copy_fasta_to_wdir:
    input:
        qc_pass="validation/checksum_status.json",
    output:
        outdir=directory("fastas/"),
    params:
        fasta_dir=config["fasta_dir"],
    message:
        "[Input validation] locally saving valid assemblies"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/copy_fasta_to_wdir.log",
    script:
        "../scripts/copy_fasta_to_wdir.py"
