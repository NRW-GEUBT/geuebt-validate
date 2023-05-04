# Assemblies are defined by checkpoint stage_fastas
# wildcard : qcpass_id
# For all passing samples, create a Json isolate sheet
# For all samples create an analysis report (json + tsv)


checkpoint stage_fastas:
    input:
        qc_pass="validation/assemblies_status.json",
    output:
        outdir=directory("staging/fastas"),
    params:
        fasta_dir=f"fastas/",
    message:
        "[Staging] Moving assemblies to staging area"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/stage_fastas.log",
    script:
        "../scripts/stage_fastas.py"


rule mlst:
    input:
        assembly="fastas/{qcpass_id}.fa",
    output:
        json="assembly_qc/mlst/{qcpass_id}.mlst.json",
    message:
        "[Staging][{wildcards.qcpass_id}] MLST scanning and sequence type"
    conda:
        "../envs/mlst.yaml"
    log:
        "logs/mlst_{qcpass_id}.log",
    shell:
        """
        exec 2> {log}
        mlst -json {output.json} {input.assembly}
        """


rule create_isolate_sheet:
    input:
        mlst="assembly_qc/mlst/{qcpass_id}.mlst.json",
        metadata="validation/metadata.json",
        assembly_qc="assembly_qc/summaries/{qcpass_id}.json",
    output:
        json="staging/isolates_sheets/{qcpass_id}.json",
    params:
        isolate_id="{qcpass_id}",
    message:
        "[Staging][{wildcards.qcpass_id}] Creating isolate datasheet"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/create_isolate_sheet_{qcpass_id}.log",
    script:
        "../scripts/create_isolate_sheet.py"


rule merge_isolate_sheets:
    input:
        qc_pass=aggregate_qcpass,
    output:
        json="staging/isolates_datasheet.json",
    message:
        "[Staging] Aggreagting isolate datasheets"
    log:
        "logs/merge_isolate_sheets.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/merge_isolate_sheets.py"


rule report_qc_status:
    input:
        metadata="validation/metadata_status.json",
        checksums="validation/checksum_status.json",
        assembly_qc="validation/assemblies_status.json",
    output:
        json="staging/validation_status.json",
    message:
        "[Staging] Creating validation status reports"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/report_qc_status.log"
    script:
        "../scripts/report_qc_status.py"
