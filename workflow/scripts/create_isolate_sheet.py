#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import json
from dateutil.parser import parse
from datetime import datetime


fmt_strs = ['%d.%m.%Y', '%Y-%m-%d']


def parse_dates(dt_str, fmt_strs=fmt_strs):
    for fmt in fmt_strs:
        try:
            return datetime.strptime(dt_str, fmt).strftime("%Y-%m-%d")
        except ValueError:
            pass
    # If not parsable, let the API reject it
    return dt_str


def main(mlst, metadata, assembly_qc, isolate_id, json_path):
    # load jsons
    with open(metadata, 'r') as fp:
        meta_dict = json.load(fp)[isolate_id]
    with open(assembly_qc, 'r') as fp:
        assembly_dict = json.load(fp)
    with open(mlst, 'r') as fp:
        mlst_dict = json.load(fp)[0]  # mlst res as list of dict
    # get relevant entries from metadata
    dictout = {key: meta_dict.get(key, None) for key in [
        "isolate_id",
        "sample_id",
        "alt_isolate_id",
        "organism",
        "third_party_owner",
        "sample_type",
        "fasta_name",
        "fasta_md5",
    ]}

    dictout["mlst"] ={
         key: mlst_dict.get(key, None) for key in ['sequence_type', 'scheme', 'alleles']
    }

    dictout["sample_info"] = {
        key: meta_dict.get(key, None) for key in [
            "isolation_org",
            "sequencing_org",
            "bioinformatics_org",
            "extraction_method",
            "library_kit",
            "sequencing_kit",
            "sequencing_instrument",
            "assembly_method",
        ]
    }

    dictout["epidata"] = {
        key: meta_dict.get(key, None) for key in [
            "collection_date",
            "isolation_date",
            "customer",
            "manufacturer",
            "collection_place",
            "collection_place_code",
            "collection_place_description",
            "collection_place_zip",
            "colleciton_place_city",
            "description",
            "manufacturer_type",
            "manufacturer_type_code",
            "manufacturer_type_description",
            "matrix",
            "matrix_code",
            "matrix_group_code",
            "processing",
            "processing_code",
            "collection_cause",
            "collection_cause_code",
            "analysis_cause",
            "analysis_cause_code",
            "control_program",
            "control_program_code",
            "lot_number",
        ]
    }
    dictout["qc_metrics"] = {
        key: assembly_dict[key] for key in [
            "seq_depth",
            "ref_coverage",
            "q30",
            "N50",
            "L50",
            "n_contigs_1kbp",
            "assembly_size",
            "GC_perc",
            "orthologs_found",
            "duplicated_orthologs",
            "majority_genus",
            "fraction_majority_genus",
            "majority_species",
            "fraction_majority_species",
        ]
    }

    # NB: fastas are renamed to match isolated_id
    dictout["fasta_name"] = f"{isolate_id}.fa"

    # Make sure all metadata except QC is parsed as string - important for (badly formatted) dates and the AVV/ADV codes!
    for k, v in dictout["epidata"].items():
        dictout["epidata"][k] = str(v) if v else None
    # And now format dates
    for field in ["collection_date", "isolation_date"]:
        if dictout["epidata"][field]:
            dictout["epidata"][field] = parse_dates(dictout["epidata"][field], fmt_strs)

    # export to json
    with open(json_path, 'w') as fp:
        json.dump(dictout, fp, indent=4)


if __name__ == '__main__':
    main(
        snakemake.input['mlst'],
        snakemake.input['metadata'],
        snakemake.input['assembly_qc'],
        snakemake.params['isolate_id'],
        snakemake.output['json']
    )
