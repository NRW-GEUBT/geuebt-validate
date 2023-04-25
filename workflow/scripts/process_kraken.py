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
        lin = tax.getAncestry(x['taxid'])
        lin.filter(ranks)
        return [n.name for n in lin]
    except KeyError:
        return ['unassigned']*len(ranks)


def rank_prop(df, rank, n=2):
    cumsum = df[['len', rank]].groupby(
        rank
    ).sum().reset_index()
    total_length = sum(cumsum['len'])
    cumsum['proportion'] = cumsum['len']/total_length
    cumsum['rank'] = rank
    return cumsum.rename(
        columns={rank: 'name'}
    ).sort_values(
        'proportion', ascending =False
    ).head(n)


def main(kraken, taxdump, json_path):
    # Load taxonomy
    taxdump = os.path.expanduser(taxdump)
    tax = txd.Taxonomy.from_taxdump(
        os.path.join(taxdump, "nodes.dmp"),
        os.path.join(taxdump, "rankedlineage.dmp")
    )
    # get at ranks
    krak = pd.read_csv(
        kraken,
        sep='\t',
        header=None,
        index_col=False,
        names=[
            'flag',
            'id',
            'taxid',
            'len',
            'lca'
        ]
    )
    krak[['species', 'genus']] = krak.apply(
        lambda x: get_ancestry(x, tax),
        result_type='expand',
        axis=1)
    # Get genus and species
    krak = krak[['len', 'species', 'genus']]
    # Species level
    species = rank_prop(
        krak, 'species', 1
    )
    pred_species = species['name'].values[0]
    frac_maj_species = species['proportion'].values[0]
    # Genus level
    genus = rank_prop(
        krak, 'genus', 1
    )
    pred_genus = genus['name'].values[0]
    frac_maj_genus = genus['proportion'].values[0]
    d = {
        'predicted_species': pred_species,
        'fraction_majority_species': frac_maj_species,
        'predicted_genus': pred_genus,
        'fraction_majority_genus': frac_maj_genus
    }
    with open(json_path, 'w') as f:
        json.dump(d, f, indent=4)


if __name__ == '__main__':
    main(
        snakemake.input['krout'],
        snakemake.params['taxdump'],
        snakemake.output['json']
    )
