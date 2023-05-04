#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import json


def main(qc_pass, json_path):
    entries = []
    for sheet in qc_pass:
        with open(sheet, "r") as fp:
            entries.append(json.load(fp))
    with open(json_path, "w") as fp:
        json.dump(entries, fp, indent=4)


if __name__ == '__main__':
    main(
        snakemake.input['qc_pass'],
        snakemake.output['json']
    )
