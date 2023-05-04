#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import os
import json
import shutil


def main(qc_pass, fasta_dir, outdir):
    # load json
    with open(qc_pass, "r") as f:
        qcp = json.load(f)
        # make output dir
    os.makedirs(outdir, exist_ok=True)
    # Copy all selected fastas locally
    for k, v in qcp.items():
        if v["STATUS"] == "PASS":
            shutil.copy(
                os.path.join(fasta_dir, f"{k}.fa"),
                os.path.join(outdir, f"{k}.fa"),
            )


if __name__ == '__main__':
    main(
        snakemake.input['qc_pass'],
        snakemake.params['fasta_dir'],
        snakemake.output['outdir']
    )
