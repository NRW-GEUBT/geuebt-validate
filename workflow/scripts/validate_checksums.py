#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


sys.stderr = open(snakemake.log[0], "w")


import os
import hashlib
import pandas as pd


# TESTS
metadata="/home/debian/NGS/NRW-geübt/geuebt-validate/.tests/data/metadata_table_test.tsv"
metadata_qc= "/home/debian/NGS/NRW-geübt/geuebt-validate/.tests/validation_test/metadata_status.tsv"


def md5(fname):
    """
    Return the hex string representation of the md5 difest for a file
    From stackoverflow user @quantumSoup (2010-08-7)
    """
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.md5.update(chunk)
    return hash_md5.hexdigest()


def validate_fasta(x, fasta_prefix=None)
    """
    Check wether the file exists, calculate md5 and compare to user value
    """
    fasta = os.path.join(fasta_prefix, x['fasta_name'])
    if os.path.exists(fasta):
        md5hex = md5(fasta)
        if md5hex == x['fasta_md5']:
            status = 'PASS'
            mssg = ''
        else:
            status = 'FAIL'
            mssg = "f{x['fasta_name']} is corrupted"
    else:
        md5hex = pd.NA
        status = 'FAIL'
        mssg = f"{x['fasta_name']} does not exist"
    return pd.Series([md5hex, status, mssg])


def main(metadata, metadata_qc, fasta_dir, json, tsv):
    # Get valid entries as list
    pass_ids = pd.read_csv(
        metadata_qc,
        sep='\t',
    )
    pass_ids = pass_ids.loc[
        pass_ids['STATUS'] == 'PASS'
    ]['isolate_id'].to_list()
    # Extract fasta paths and md5
    fastas = pd.read_csv(
        metadata,
        sep='\t',
        usecols=['isolate_id', 'fasta_name', 'fasta_md5']
    )
    fastas = fastas.loc[fastas['isolate_id'].isin(pass_ids)]
    # Check that file exists and calculate md5
    fastas[['local_md5','STATUS','MESSAGES']] = fastas.apply(
        validate_fasta,
        fasta_prefix=fasta_dir,
        axis=1
    )
    # Export tsv
    fastas.to_csv(
        tsv_path,
        sep='\t',
        header=True,
        index=False
    )
    # Format and export JSON
    df_to_dict = set_index('isolate_id').to_json(orient='index')
    with open(json_path, 'w') as f:
        json.dump(df_to_dict, f, indent=4)


if __name__ == '__main__':
    main(
        snakemake.input['metadata'],
        snakemake.input['metadata_qc'],
        snakemake.params['fasta_dir'],
        snakemake.output['json'],
        snakemake.output['tsv']
    )
