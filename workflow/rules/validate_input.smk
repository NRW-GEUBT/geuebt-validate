# Validate user-submitted metadata with schema
# Then validate md5 checksums of fasta files


rule validate_metadata:
    output:
        json="validation/metadata_status.json",
        tsv="validation/metadata_status.tsv",
    params:
        schema=f"{workflow.basedir}/schema/metadata.schema.json",
        metadata=config['metadata'],
    message:
        "[Input validation] validating metadata"
    conda:
        "../envs/jsonschema.yaml"
    log:
        "logs/validate_metadata.log"
    script:
        "../scripts/validate_metadata.py"


rule validate_checksum:
    input:
        metadata_qc="validation/metadata_status.json",
    output:
        json="validation/checksum_status.json",
        tsv="validation/checksum_status.tsv",
    params:
        fasta_dir=config['fasta_dir'],
        metadata=config['metadata'],
    message:
        "[Input validation] validating md5 checksums"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/validate_checksum.log"
    script:
        "../scripts/validate_checksums.py"


checkpoint copy_fasta_to_wdir:
    input:
        qc_pass="validation/checksum_status.json",
    output:
        dir=directory("fastas/"),
    params:
        fasta_dir=config['fasta_dir'],
    message:
        "[Input validation] locally saving valid assemblies"
    log:
        "logs/copy_fasta_to_wdir.log"
    run:
        import os
        import json
        import shutil


        # load json
        with open(input.qc_pass, 'r') as f:
            qcp = json.load(f)
        # make output dir
        os.makedirs(output.dir, exist_ok=True)
        # Copy all selected fastas locally
        for k, v in qcp.items():
            if v['STATUS'] == 'PASS':
                shutil.copy(
                    os.path.join(params.fasta_dir, v['fasta_name']),
                    os.path.join(output.dir, f"{k}.fa")
                )

