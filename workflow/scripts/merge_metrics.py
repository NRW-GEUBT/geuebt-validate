#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import json
import pandas as pd


def main(quast, kraken, metadata, jsonpath, busco, isolate_id):
    res = dict()
    # Metadata - access isolate via isolate_id
    with open(metadata, "r") as f:
        meta = json.load(f)
    res["isolate_id"] = str(meta[isolate_id]["isolate_id"])
    res["expect_species"] = str(meta[isolate_id]["organism"])
    res["seq_depth"] = float(meta[isolate_id]["seq_depth"])
    res["ref_coverage"] = float(meta[isolate_id]["ref_coverage"])
    res["q30"] = float(meta[isolate_id]["q30"])
    # QUAST - as TSV
    quast = pd.read_csv(quast, index_col=False, sep="\t")
    res["N50"] = int(quast["N50"].values[0])
    res["L50"] = int(quast["L50"].values[0])
    res["n_contigs_1kbp"] = int(quast["# contigs (>= 1000 bp)"].values[0])
    res["assembly_size"] = int(quast["Total length"].values[0])
    res["GC_perc"] = float(quast["GC (%)"].values[0])
    # BUSCO
    with open(busco, "r") as f:
        busco = json.load(f)
    res["orthologs_found"] = float(busco["results"]["Complete percentage"])
    res["duplicated_orthologs"] = float(busco["results"]["Multi copy percentage"])
    # Kraken
    with open(kraken, "r") as f:
        kraken = json.load(f)
    res["majority_genus"] = str(kraken["predicted_genus"])
    res["fraction_majority_genus"] = float(kraken["fraction_majority_genus"])
    res["majority_species"] = str(kraken["predicted_species"])
    res["fraction_majority_species"] = float(kraken["fraction_majority_species"])
    # Output
    with open(jsonpath, "w") as f:
        json.dump(res, f, indent=4)


if __name__ == '__main__':
    main(
        snakemake.input['quast'],
        snakemake.input['kraken'],
        snakemake.input['metadata'],
        snakemake.output['json'],
        snakemake.params['busco'],
        snakemake.params['isolate_id']
    )
