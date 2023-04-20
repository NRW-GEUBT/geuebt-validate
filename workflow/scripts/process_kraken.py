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
import taxidTools as txd


def get_ancestry(x, tax, ranks=["genus", "species"]):
    lin = tax.getAncestry(x[2])
    lin.filter(ranks)
    return lin


def main(kraken, taxdump, json):
    # Load taxonomy
    tax = txd.Taxonomy.from_taxdump(
        os.path.join(taxdump, "nodes.dmp"),
        os.path.join(taxdump, "rankedlineage.dmp")
    )
    # get at ranks
    krak = pd.read_csv(kraken, sep='\t', header=None)
    krak[["genus", "species"]] = krak.apply(
        lambda x: get_ancestry(x, tax), 
        axis=1) # does that even work???
    # get cumsum for spec and genus
    # keep top 2 for each
    # to json


if __name__ == '__main__':
    main(
        snakemake.input['krout'],
        snakemake.params['taxdump'],
        snakemake.output['json']
    )

