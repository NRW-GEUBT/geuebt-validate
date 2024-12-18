# Assemblies are defined by checkpoint copy_fasta_to_wdir
# wildcard : metapass_id
# check assemblies individually and aggregate results in a JSON
# then validate qc values


import os


rule assembly_qc_quast:
    input:
        assembly="fastas/{isolate_id}.fa",
    output:
        report="assembly_qc/quast/{isolate_id}/transposed_report.tsv",
        outdir=directory("assembly_qc/quast/{isolate_id}"),
    params:
        min_contig_length=config["min_contig_length"],
    message:
        "[Assembly quality][{wildcards.isolate_id}] Calculating assembly metrics with QUAST"
    conda:
        "../envs/quast.yaml"
    log:
        "logs/assembly_qc_quast_{isolate_id}.log",
    threads: config["max_threads_per_job"]
    shell:
        """
        quast -o {output.outdir} -f \
            --min-contig {params.min_contig_length} \
            --threads {threads} \
            --no-html \
            {input.assembly} > {log} 2>&1
        """


rule assembly_qc_busco:
    input:
        assembly="fastas/{isolate_id}.fa",
    output:
        # Naming of outputs is dependent on the database name 
        # using a touch flag to ensure this rule is executed before
        # trying to run rule merge_metrics and specifying JSON output 
        # there with input funciton
        flag=touch("assembly_qc/busco/{isolate_id}/done.flag"),
        outdir=directory("assembly_qc/busco/{isolate_id}"),
    params:
        busco_db=os.path.expanduser(f"~/.nrw-geuebt/geuebt-validate-{version}/busco/bacteria_odb10"),
    message:
        "[Assembly quality][{wildcards.isolate_id}] Detecting and counting conserved genes with BUSCO"
    conda:
        "../envs/busco.yaml"
    log:
        "logs/assembly_qc_busco_{isolate_id}.log",
    threads: config["max_threads_per_job"]
    shell:
        """
        busco -f -i {input.assembly} \
            -l {params.busco_db} \
            -o {output.outdir} \
            -m genome \
            --offline \
            --cpu {threads} > {log} 2>&1
        """


rule assembly_qc_kraken2:
    input:
        assembly="fastas/{isolate_id}.fa",
    output:
        report="assembly_qc/kraken2/{isolate_id}.kreport",
        krout="assembly_qc/kraken2/{isolate_id}.kraken",
    params:
        kraken2_db=os.path.expanduser(f"~/.nrw-geuebt/geuebt-validate-{version}/kraken/kraken2_standard8"),
    message:
        "[Assembly quality][{wildcards.isolate_id}] Taxonomic classification of kmers with KRAKEN2"
    conda:
        "../envs/kraken2.yaml"
    log:
        "logs/assembly_qc_kraken2_{isolate_id}.log",
    threads: config["max_threads_per_job"]
    shell:
        """
        exec 2> {log}
        kraken2 --db {params.kraken2_db} \
            --threads {threads} \
            --output {output.krout} \
            --report {output.report} \
            --confidence 0 \
            {input.assembly}
        """


rule process_kraken:
    input:
        krout="assembly_qc/kraken2/{isolate_id}.kraken",
    output:
        json="assembly_qc/kraken2/{isolate_id}.kraken.json",
    params:
        taxdump=os.path.expanduser(f"~/.nrw-geuebt/geuebt-validate-{version}/taxdump/"),
    message:
        "[Assembly quality][{wildcards.isolate_id}] Calculating taxonomic ditribution of assembly"
    conda:
        "../envs/taxidtools.yaml"
    log:
        "logs/process_kraken_{isolate_id}.log",
    script:
        "../scripts/process_kraken.py"


rule merge_metrics:
    input:
        buscoflag="assembly_qc/busco/{isolate_id}/done.flag",
        quast="assembly_qc/quast/{isolate_id}/transposed_report.tsv",
        kraken="assembly_qc/kraken2/{isolate_id}.kraken.json",
        metadata="validation/metadata.json",
    output:
        json="assembly_qc/summaries/{isolate_id}.json",
    params:
        # Getting name for function in params as there is no output rule for it
        # see rule assembly_qc_busco
        busco=get_busco_out_name,
        isolate_id="{isolate_id}",
    message:
        "[Assembly quality][{wildcards.isolate_id}] Merging assembly QC metrics"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/merge_metrics_{isolate_id}.log",
    script:
        "../scripts/merge_metrics.py"


rule mlst:
    input:
        assembly="fastas/{isolate_id}.fa",
    output:
        json="assembly_qc/mlst/{isolate_id}.mlst.json",
    message:
        "[Staging][{wildcards.isolate_id}] MLST scanning and sequence type"
    conda:
        "../envs/mlst.yaml"
    log:
        "logs/mlst_{isolate_id}.log",
    shell:
        """
        exec 2> {log}
        mlst -json {output.json} {input.assembly} > {log}
        """


rule create_isolate_sheet:
    input:
        mlst="assembly_qc/mlst/{isolate_id}.mlst.json",
        metadata="validation/metadata.json",
        assembly_qc="assembly_qc/summaries/{isolate_id}.json",
    output:
        json="isolates_sheets/{isolate_id}.json",
    params:
        isolate_id="{isolate_id}",
    message:
        "[Staging][{wildcards.isolate_id}] Creating isolate datasheet"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/create_isolate_sheet_{isolate_id}.log",
    script:
        "../scripts/create_isolate_sheet.py"
