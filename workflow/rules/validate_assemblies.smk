# Assemblies are defined by checkpoint copy_fasta_to_wdir
# wildcard : metapass_id
# check assemblies individually and aggregate results in a JSON
# then validate qc values


rule assembly_qc_quast:
    input:
        assembly="fastas/{metapass_id}.fa",
    output:
        report="assembly_qc/quast/{metapass_id}/transposed_report.tsv",
        outdir=directory("assembly_qc/quast/{metapass_id}"),
    params:
        min_contig_length=config["min_contig_length"],
    message:
        "[Assembly quality][{wildcards.metapass_id}] Calculating assembly metrics with QUAST"
    conda:
        "../envs/quast.yaml"
    log:
        "logs/assembly_qc_quast_{metapass_id}.log",
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
        assembly="fastas/{metapass_id}.fa",
    output:
        # Naming of outputs is dependent on the database name 
        # using a touch flag to ensure this rule is executed before
        # trying to run rule merge_metrics and specifying JSON output 
        # there with input funciton
        flag=touch("assembly_qc/busco/{metapass_id}/done.flag"),
        outdir=directory("assembly_qc/busco/{metapass_id}"),
    params:
        busco_db=config["busco_db"],
    message:
        "[Assembly quality][{wildcards.metapass_id}] Detecting and counting conserved genes with BUSCO"
    conda:
        "../envs/busco.yaml"
    log:
        "logs/assembly_qc_busco_{metapass_id}.log",
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
        assembly="fastas/{metapass_id}.fa",
    output:
        report="assembly_qc/kraken2/{metapass_id}.kreport",
        krout="assembly_qc/kraken2/{metapass_id}.kraken",
    params:
        kraken2_db=config["kraken2_db"],
    message:
        "[Assembly quality][{wildcards.metapass_id}] Taxonomic classification of kmers with KRAKEN2"
    conda:
        "../envs/kraken2.yaml"
    log:
        "logs/assembly_qc_kraken2_{metapass_id}.log",
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
        krout="assembly_qc/kraken2/{metapass_id}.kraken",
    output:
        json="assembly_qc/kraken2/{metapass_id}.kraken.json",
    params:
        taxdump=config["taxdump"],
    message:
        "[Assembly quality][{wildcards.metapass_id}] Calculating taxonomic ditribution of assembly"
    conda:
        "../envs/taxidtools.yaml"
    log:
        "logs/process_kraken_{metapass_id}.log",
    script:
        "../scripts/process_kraken.py"


rule merge_metrics:
    input:
        buscoflag="assembly_qc/busco/{metapass_id}/done.flag",
        quast="assembly_qc/quast/{metapass_id}/transposed_report.tsv",
        kraken="assembly_qc/kraken2/{metapass_id}.kraken.json",
        metadata="validation/metadata.json",
    output:
        json="assembly_qc/summaries/{metapass_id}.json",
    params:
        # Getting name for function in params as there is no output rule for it
        # see rule assembly_qc_busco
        busco=get_busco_out_name,
        isolate_id="{metapass_id}",
    message:
        "[Assembly quality][{wildcards.metapass_id}] Merging assembly QC metrics"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/merge_metrics_{metapass_id}.log",
    script:
        "../scripts/merge_metrics.py"


rule aggregate_metrics:
    input:
        metrics=aggregate_metapass,
    output:
        mergedout="assembly_qc/assembly_metrics.json",
    message:
        "[Assembly quality] Merging assembly QC metrics"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/aggregate_metrics.log",
    script:
        "../scripts/aggregate_metrics.py"


rule validate_assembly_qc:
    input:
        metrics="assembly_qc/assembly_metrics.json",
    output:
        json="validation/assemblies_status.json",
        tsv="validation/assemblies_status.tsv",
    params:
        schema=f"{workflow.basedir}/schema/assembly_qc.schema.json",
    message:
        "[Assembly quality] Validating assemblies quality"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/validate_assembly_qc.log",
    script:
        "../scripts/validate_assemblies.py"
