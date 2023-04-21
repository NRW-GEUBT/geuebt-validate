# Assemblies are defined by checkpoint copy_fasta_to_wdir
# wildcard : metapass_id
# check assemblies individually and aggregate results in a JSON
# then validate qc values


rule assembly_qc_quast:
    input:
        assembly="fastas/{metapass_id}.fa",
    output:
        report="assembly_qc/quast/{metapass_id}/report.tsv",
    params:
        outdir="assembly_qc/quast/{metapass_id}",
        min_contig_length=config['min_contig_length'],
    message:
        "[Assembly quality][{wilcards.metapass_id}] Calculating assembly metrics"
    conda:
        "../envs/quast.yaml"
    log:
        "logs/assembly_qc_quast_{metapass_id}.log"
    threads:
        config['max_threads_per_job']
    shell:
        """
        exec 2> {log}
        quast -o {params.outdir} \
            --min-contig {params.min_contig_length} \
            --threads {threads} \
            --no-html \
            {input.assembly}
        """


rule assembly_qc_busco:
    input:
        assembly="fastas/{metapass_id}.fa",
    output:
        json="assembly_qc/busco/{metapass_id}/short_summary.{metapass_id}.json",
    params:
        outdir="assembly_qc/busco/{metapass_id}",
        busco_db=config['busco_db'],
    message:
        "[Assembly quality][{wilcards.metapass_id}] Deteting and counting conserved genes"
    conda:
        "../envs/busco.yaml"
    log:
        "logs/assembly_qc_busco_{metapass_id}.log"
    threads:
        config['max_threads_per_job']
    shell:
        """
        busco -i {input.assembly} \
            -l {params.busco_db} \
            -o {params.outdir} \
            -m genome \
            --cpu {threads}
        """


rule assembly_qc_kraken2:
    input:
        assembly="fastas/{metapass_id}.fa",
    output:
        report="assembly_qc/kraken2/{metapass_id}.kreport",
        krout="assembly_qc/kraken2/{metapass_id}.kraken",
    params:
        kraken2_db=config['kraken2_db'],
    message:
        "[Assembly quality][{wilcards.metapass_id}] Taxonomic classification of kmers"
    conda:
        "../envs/kraken2.yaml"
    log:
        "logs/assembly_qc_kraken2_{metapass_id}.log"
    threads:
        config['max_threads_per_job']
    shell:
        """
        exec 2> {log}
        kraken2 --db {params.kraken2_db} \
            --threads {threads} \
            --output {output.krout} \
            --report {output.report} \
            --confidence 0 
            {input.assembly}
        """


rule process_kraken:
    input:
        krout="assembly_qc/kraken2/{metapass_id}.kraken",
    output:
        json="assembly_qc/kraken2/{metapass_id}.kraken.json",
    params:
        taxdump=config['taxdump'],
    message:
        "[Assembly quality][{wilcards.metapass_id}] Calculating taxonomic ditribution of assembly"
    conda:
        "../envs/taxidtools.yaml"
    log:
        "logs/process_kraken_{metapass_id}.log"
    script:
        "../scripts/process_kraken.py"


rule mlst:
    input:
        assembly="fastas/{metapass_id}.fa",
    output:
        json="assembly_qc/mlst/{metapass_id}.mlst.json",
    message:
        "[Assembly quality][{wildcards.metapass_id}] MLST scanning and sequence type"
    conda:
        "../envs/mlst.yaml"
    log:
        "logs/mlst_{metapass_id}.log"
    shell:
        """
        exec 2> {log}
        mlst -json {output.json} {input.assembly}
        """


rule merge_metrics:
    input:
        busco="assembly_qc/busco/{metapass_id}/short_summary.{metapass_id}.json",
        quast="assembly_qc/quast/{metapass_id}/report.tsv",
        kraken2="assembly_qc/kraken2/{metapass_id}.kraken.json",
        mlst="assembly_qc/mlst/{metapass_id}.mlst.json",
    output:
        json="assembly_qc/summaries/{metapass_id}.json"
    message:
        "[Assembly quality][{wildcards.metapass_id}] Merging assembly QC metrics"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/merge_metrics_{metapass_id}.log"
    shell:
        """
        touch {output.json}
        """
        # QUAST report is a TSV!


rule aggregate_metrics:
    input:
        metrics=aggregate_metapass,
    output:
        metrics="assembly_qc/assembly_metrics.json",
    shell:
        """
        touch {output.metrics}
        """


# rule validate_assembly_qc:
    # Compare merge Json to schema and report deviations

