#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import json


def main(mlst, metadata, assembly_qc, isolate_id, json_path):
    # load jsons
    with open(metadata, 'r') as fp:
        meta_dict = json.load(fp)
    with open(assembly_qc, 'r') as fp:
        assembly_dict = json.load(fp)
    with open(mlst, 'r') as fp:
        mlst_dict = json.load(fp)[0]  # mlst res as list of dict
    # get relevant entries from metadata
    dictout = meta_dict[isolate_id]
    # NB: fastas are renamed to match isolated_id
    dictout["fasta_name"] = f"{isolate_id}.fa"
    # reformat metadata
    dictout['isolation'] = {
        'org_name': dictout.pop('isolation_org'),
    }
    dictout['sequencing'] = {
        "org_name": dictout.pop('sequencing_org'),
        "extraction_method": dictout.pop('extraction_method'),
        "library_kit": dictout.pop('library_kit'),
        "sequencing_instrument": dictout.pop('sequencing_instrument'),
    }
    dictout['bioinformatics'] = {
        "org_name": dictout.pop('bioinformatics_org'),
        "assembly_method": dictout.pop('assembly_method'),
    }
    dictout['epidata'] = {
        "collection_date": dictout.pop('collection_date'),
        "collection_place": dictout.pop('collection_place'),
        "collection_place_code": dictout.pop('collection_place_code'),
        "cause": dictout.pop('collection_cause'),
        "customer": dictout.pop('customer'),
        "manufacturer": dictout.pop('manufacturer'),
        "description": dictout.pop('description'),
        "manufacturer_type": dictout.pop('manufacturer_type'),
        "manufacturer_type_code": dictout.pop('manufacturer_type_code'),
        "matrix": dictout.pop('matrix'),
        "matrix_code": dictout.pop('matrix_code'),
        "lot_number": dictout.pop('lot_number'),
    }
    # Add MLST infos
    dictout['mlst'] = {
        key: mlst_dict[key] for key in ['sequence_type', 'scheme', 'alleles']
    }
    # Add QC metrics
    dictout['qc_metrics'] = {
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
    # cleaning duplicated keys at top level
    for key in ["seq_depth", "ref_coverage", "q30"]:
        del dictout[key]
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
