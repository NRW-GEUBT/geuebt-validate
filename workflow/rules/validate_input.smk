# Validate user-submitted metadata with schema
# Then validate md5 checksums of fasta files


rule validate_metadata:
    output:
        json="validation/metadata_status.json",
        tsv="validation/metadata_status.tsv",
        metadata_json="validation/metadata.json",
    params:
        schema=f"{workflow.basedir}/schema/metadata.schema.json",
        metadata=config["metadata"],
    message:
        "[Input validation] validating metadata"
    conda:
        "../envs/pydantic.yaml"
    log:
        "logs/validate_metadata.log",
    script:
        "../scripts/validate_metadata.py"


rule validate_checksum:
    input:
        metadata_qc="validation/metadata_status.json",
    output:
        json="validation/checksum_status.json",
        tsv="validation/checksum_status.tsv",
    params:
        fasta_dir=config["fasta_dir"],
        metadata=config["metadata"],
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
