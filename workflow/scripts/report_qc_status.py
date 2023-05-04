#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import json


def main(metadata_status, checksum_status, assembly_status, json_path):
    # Load Jsons
    with open(metadata_status, 'r') as fp:
        meta = json.load(fp)
    with open(checksum_status, 'r') as fp:
        checksum = json.load(fp)
    with open(assembly_status, 'r') as fp:
        ass = json.load(fp)
    # merge based on isolate_id
    dout = {}
    for k in meta.keys():
        dout[k] = meta[k]
        # Order matters here! Failed samples are not taken to the next step
        for d in [checksum, ass]: 
            try:
                dout[k]['STATUS'] = d[k]['STATUS']
                dout[k]['MESSAGES'].extend(d[k]['MESSAGES'])
            except KeyError:
                pass
    # dump JSON
    print(dout)
    with open(json_path, 'w') as fp:
        json.dump(dout, fp, indent=4)


if __name__ == '__main__':
    main(
        snakemake.input['metadata'],
        snakemake.input['checksums'],
        snakemake.input['assembly_qc'],
        snakemake.output['json']
    )
