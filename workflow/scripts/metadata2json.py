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
import numpy as np


def main(metadata, metadata_json):
    # load metadata as dataframe
    metatable = pd.read_csv(
        metadata,
        sep='\t',
        header=0,
        index_col=False
    )
    metatable = metatable.replace({np.nan: None})
    records = metatable.to_dict(orient='records')
    parsed_records = {r["isolate_id"]: r for r in records}

    # export valid metadata to Json file
    with open(metadata_json, 'w') as f:
        json.dump(parsed_records, f, indent=4)


if __name__ == '__main__':
    main(
        snakemake.params['metadata'],
        snakemake.output['metadata_json']
    )
