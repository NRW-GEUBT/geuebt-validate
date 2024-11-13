# Gather all isolates sheets and sequences
# Post isolate_sheet and check status
# if response 422 -> parse error message andwrite to qc sheet
# if response 200 -> Pass
# wirte to QC sheet and post sequence
# copy passing sequences and isolate sheets to staging area to be passed on to CORE


rule post_data:
    input:
        isolate_sheets=aggregate_isolate_id,
        checksum_qc="validation/checksum_status.json",
    output:
        qc="staging/validation_status.json",
        isolate_sheet="staging/isolates_datasheet.json",
        isolate_sheets_dir=directory("staging/isolates_sheets"),
        fastas=directory("staging/fastas"),
    params:
        fasta_store="fastas",
        uri=config["mongo_uri"],
    message:
        "[Staging] Posting samples and sequences to database"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/post_data.log"
    script:
        "../scripts/post_data.py"

