rule validate_metadata:
    input:
        schema=f"{workflow.basedir}/workflow/schema/metadata.schema.json",
        metadata=config['metadata'],
    output:
        json="validation/metadata_status.json",
        tsv="validation/metadata_status.tsv",
    message:
        "[Input validation] checking metadata"
    conda:
        "../envs/jsonschema.yaml"
    log:
        "logs/validate_metadata.log"
    script:
        "../scripts/validate_metadata.py"


rule validate_checksum:
    input:
        metadata=config['metadata'],
        metadata_qc="validation/metadata_status.tsv",
    output:
        json="validation/checksum_status.json",
        tsv="validation/checksum_status.tsv",
    params:
        fasta_dir=config['fasta_dir'],
    message:
        "[Input validation] exporting md5 checksums"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/export_md5_table.log"
    script:
        "../scripts/validate_checksums.py"

