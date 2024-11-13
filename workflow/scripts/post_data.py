#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import os
import shutil
import json
import requests
from urllib.parse import urljoin


def main(
        isolate_sheets,
        checksum_qc,
        qc_out,
        isolate_sheet_out,
        isolate_sheet_dir,
        fastas_dir,
        fasta_store,
        uri,
        ):
    
    # Create dirs
    os.makedirs(isolate_sheet_dir, exist_ok=True)
    os.makedirs(fastas_dir, exist_ok=True)

    # load chekcsum qc
    with open(checksum_qc, "r") as fi:
        checksums = json.load(fi)
    qc = {
        id: {"STATUS": record["STATUS"], "MESSAGES": record["MESSAGES"]} 
        for id, record in checksums.items()
    }
    
    # iterate over samples
    merged = []

    for isolate_file in isolate_sheets:
        # parse json
        with open(isolate_file, "r") as fi:
            data = json.load(fi)
        isolate_id = data["isolate_id"]
        
        # POST to server and get code
        response = requests.post(urljoin(uri, "isolates"), json=data)

        # On success
        if response.status_code == 200:
            qc[isolate_id]["MESSAGES"].append(response.json()["message"])
            
            # Add to merged
            merged.append(data)

            # POST fasta
            fasta_path = os.path.join(fasta_store, f"{isolate_id}.fa")
            with open(fasta_path, "r") as fasta:
                fdump = fasta.read()
            fasta_json = {
                "isolate_id": isolate_file,
                "sequence_type": "fasta",
                "sequence": fdump
            }
            fasta_response = requests.post(urljoin(uri, "sequences"), json=fasta_json)
            if fasta_response.status_code == 200:
                qc[isolate_id]["MESSAGES"].append(fasta_response.json()["message"])
            else:  # should never happen
                qc[isolate_id]["STATUS"] = "WARNING"
                qc[isolate_id]["MESSAGES"].append(
                    f"An unexpected error occured while posting sequence record."
                    f"Status: {fasta_response.status_code}."
                    f"Body: {json.dumps(fasta_response.text)}"
                )

            # copy to staging
            shutil.copy(isolate_file, isolate_sheet_dir)
            shutil.copy(fasta_path, fastas_dir)

        # On errors
        elif response.status_code == 422: 
            err_type = response.json()["detail"][0]["type"]
            err_msg = response.json()["detail"][0]["msg"]
            err_field = "/".join(response.json()["detail"][0]["loc"])
            nice_error = f"VALIDATION ERROR '{err_type}': {err_msg}; for field: '{err_field}'"
            qc[isolate_id]["STATUS"] = "FAIL"
            qc[isolate_id]["MESSAGES"].append(nice_error)

        else:
            qc[isolate_id]["STATUS"] = "FAIL"
            qc[isolate_id]["MESSAGES"].append(
                f"An unexpected error has occured, contact a geuebt admin."
                f"Status: {response.status_code}."
                f"Body: {response.text}"
                )
    
    # Write merged outputs
    with open(qc_out, "w") as fo:
        json.dump(qc, fo, indent = 4)
    with open(isolate_sheet_out, "w") as fo:
        json.dump(merged, fo, indent=4)


if __name__ == '__main__':
    main(
        snakemake.input["isolate_sheets"],
        snakemake.input["checksum_qc"],
        snakemake.output["qc"],
        snakemake.output["isolate_sheet"],
        snakemake.output["isolate_sheets_dir"],
        snakemake.output["fastas"],
        snakemake.params["fasta_store"],
        snakemake.params["uri"],
    )
