#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import json


def main(metrics, mergedout):
    merged = []
    for fp in metrics:
        with open(fp, "r") as f:
            merged.append(json.load(f))
    with open(mergedout, "w") as f:
        json.dump(merged, f, indent=4)


if __name__ == '__main__':
    main(
        snakemake.input['metrics'],
        snakemake.output['mergedout']
    )
