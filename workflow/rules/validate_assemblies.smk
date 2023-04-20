# Assemblies are defined by checkpoint copy_fasta_to_wdir
# check assemblies individually and aggregate results in a JSON
# then validate qc values


rule assembly_qc_quast:
    input:
        assembly="fastas/{isolate_id}.fa",
    output:
        report="assembly_qc/quast/{isolate_id}/report.tsv"
    params:
        outdir="assembly_qc/quast/{isolate_id},
        min_contig_length=config['min_contig_length']
    message:
        "[Assembly quality][{wilcards.isolate_id}] Calculating assembly metrics"
    conda:
        "../envs/quast.yaml"
    log:
        "logs/assembly_qc_quast_{isolate_id}.log"
    threads:
        config['max_threads_per_job']
    shell:
        """
        exec 2> {log}
        quast -o {params.outdir} \
            --min-contig {params.min_contig_length} \
            --threads {threads} \
            {input.assembly}
        """


rule assembly_qc_busco:
    input:
        assembly="fastas/{isolate_id}.fa",
    output:
        json="assembly_qc/busco/{isolate_id}/short_summary.{isolate_id}.json"
    params:
        outdir="assembly_qc/busco/{isolate_id},
        busco_db="/var/lib/busco/bacteria_odb10"
    message:
        "[Assembly quality][{wilcards.isolate_id}] Deteting and counting conserved genes"
    conda:
        "../envs/busco.yaml"
    log:
        "logs/assembly_qc_busco_{isolate_id}.log"
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
        assembly="fastas/{isolate_id}.fa",
    output:
        report="assembly_qc/kraken2/{isolate_id}.kreport",
        krout="assembly_qc/kraken2/{isolate_id}.kraken",
    params:
        kraken2_db="/var/lib/kraken2/minikraken2",
    message:
        "[Assembly quality][{wilcards.isolate_id}] Taxonomic classification of kmers"
    conda:
        "../envs/kraken2.yaml"
    log:
        "logs/assembly_qc_kraken2_{isolate_id}.log"
    threads:
        config['max_threads_per_job']
    shell:
        """
        exec 2> {log}
        kraken2 --db {params.kraken2_db} \
            --threads {threads} \
            --output {output.krout} \
            --report {output.report} \
            {input.assembly}
        """


rule aggregate_assembly_qc:
    input:
        quast=lambda w:aggregate_qcpass(w, "assembly_qc/quast/{isolate_id}/report.tsv")
    output:
        json="assembly_qc/assembly_metrics.json"
        tsv="assembly_qc/assembly_metrics.tsv"
    message:
        "[Assembly quality] Aggregating assembly metrics"
    log:
        "logs/aggregate_assembly_qc.log"
    run:
        inport json
        import pandas as pd


        # Read and concat all filesin a df
        df = pd.concat(
            (pd.read_csv(f, sep='\t', header=None).T for f in input.quast),
            ignore_index=True
        ).rename(columns={"Assembly": "isolate_id"})
        # Export to tsv
        df.to_csv(output.tsv, sep='\t', header=True, index=False)
        # Convert and export to JSON
        df_to_dict = json.loads(
            df.set_index('isolate_id').to_json(orient='index')
        )
        with open(output.json, 'w') as f:
            json.dump(df_to_dict, f, indent=4)


rule aggregate_kraken:
    # Aggregate rule for kraken output
    # Get cumlength of contigs per species/genus 
    # uses taxidTools






rule validate_assembly_qc:
