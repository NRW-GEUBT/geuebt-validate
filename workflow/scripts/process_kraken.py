#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Get frac majority species
# And predicted genus


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import os
import json
import pandas as pd
import taxidTools as txd


def get_ancestry(x, tax, ranks=['species', 'genus']):
    try:
        lin = tax.getAncestry(x['len'])
        lin.filter(ranks)
        return [n.name for n in lin]
    except KeyError:
        return ['unassigned']*len(ranks)


def rank_prop(df, rank, n=2):
    cumsum = df.loc[['len', rank]].groupby(
        rank
    ).sum()
    total_length = sum(cumsum['len'])
    return cumsum.apply(
        lambda x: pd.Series(
            [rank, x[rank], x['len']/total_length],
            index=['rank', 'name', 'proportion']
        ),
        axis=1
    ).sort_values(
        'proportion',
        ascending=False
    ).head(n)


def main(kraken, taxdump, json_path):
    # Load taxonomy
    tax = txd.Taxonomy.from_taxdump(
        os.path.join(taxdump, "nodes.dmp"),
        os.path.join(taxdump, "rankedlineage.dmp")
    )
    # get at ranks
    krak = pd.read_csv(
        kraken,
        sep='\t',
        header=None,
        names=[
            'flag',
            'id',
            'taxid',
            'len'
            'lca'
        ]
    )
    krak[['species', 'genus']] = krak.apply(
        lambda x: get_ancestry(x, tax),
        result_type='expand',
        axis=1)
    # Get predicted genus (majority taxid at genus level)
    pred_genus = rank_prop(
        krak, 'genus', 1
    ).loc['genus']
    # Get fraction majority species
    frac_maj_spec = rank_prop(
        krak, 'species', 1
    ).loc['proportion']
    d = {
        'predicted_genus': pred_genus,
        'fraction_majority_species': frac_maj_spec
    }
    with open(json_path, 'w') as f:
        json.dump(d, f, indent=4)


if __name__ == '__main__':
    main(
        snakemake.input['krout'],
        snakemake.params['taxdump'],
        snakemake.output['json']
    )
