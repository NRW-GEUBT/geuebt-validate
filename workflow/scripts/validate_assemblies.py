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
from enum import Enum
from pydantic_core import PydanticCustomError
from pydantic import (
    BaseModel,
    ValidationError,
    model_validator,
    field_validator,
)


class OrganismEnum(str, Enum):
    """
    Define accepted organisms
    """
    listeria = 'Listeria monocytogenes'
    salmonella = 'Salmonella enterica'
    ecoli = 'Escherichia coli'
    campy = 'Campylobacter spp.'


class AssemblyQC(BaseModel):
    """
    Implements assembly QC thresholds from ASU L 00.00-183 2023:03
    """
    isolate_id: str
    expect_species: OrganismEnum
    n_contigs_1kbp: int
    assembly_size: int
    GC_perc: float
    orthologs_found: float
    duplicated_orthologs: float
    fraction_majority_genus: float
    majority_genus: str

    @field_validator('fraction_majority_genus')
    @classmethod
    def check_fractions(cls, v: float) -> float:
        if v < 0 or v > 1:
            raise ValueError("must be given as a fraction between 0 and 1")
        return v

    @field_validator('GC_perc', 'orthologs_found', 'duplicated_orthologs')
    @classmethod
    def check_percent(cls, v: float) -> float:
        if v < 0 or v > 100:
            raise ValueError("must be given as a percent between 0 and 100")
        return v

    @model_validator(mode='after')
    def check_species_specific_assembly_size(self) -> 'AssemblyQC':
        expect = self.expect_species
        assembly = self.assembly_size
        qcdict = {
            'Listeria monocytogenes': [2700000, 3200000],
            'Salmonella enterica': [4300000, 5200000],
            'Escherichia coli': [4500000, 5900000],
            'Campylobacter spp.': [1500000, 1900000],
        }
        if assembly < qcdict[expect][0] or assembly > qcdict[expect][1]:
            raise PydanticCustomError(
                "value_error",
                f"Value error: 'assembly_size' for '{expect}' must be between "
                f"{qcdict[expect][0]} and {qcdict[expect][1]}, got: {assembly}",
            )
        return self

    @model_validator(mode='after')
    def check_species_specific_orthologs_found(self) -> 'AssemblyQC':
        expect = self.expect_species
        ortho = self.orthologs_found
        qcdict = {
            'Listeria monocytogenes': 95,
            'Salmonella enterica': 95,
            'Escherichia coli': 95,
            'Campylobacter spp.': 80,
        }
        if ortho < qcdict[expect]:
            raise PydanticCustomError(
                "value_error",
                f"Value error: 'orthologs_found' for '{expect}' must be at least {qcdict[expect]}, got: {ortho}",
            )
        return self

    @model_validator(mode='after')
    def check_species_specific_duplicated_orthologs(self) -> 'AssemblyQC':
        expect = self.expect_species
        dup = self.duplicated_orthologs
        qcdict = {
            'Listeria monocytogenes': 5,
            'Salmonella enterica': 5,
            'Escherichia coli': 5,
            'Campylobacter spp.': 5,
        }
        if dup > qcdict[expect]:
            raise PydanticCustomError(
                "value_error",
                f"Value error: 'duplicated_orthologs' for '{expect}' must be at most {qcdict[expect]}, got: {dup}",
            )
        return self

    @model_validator(mode='after')
    def check_species_specific_fraction_majority_genus(self) -> 'AssemblyQC':
        expect = self.expect_species
        frac = self.fraction_majority_genus
        qcdict = {
            'Listeria monocytogenes': 0.95,
            'Salmonella enterica': 0.95,
            'Escherichia coli': 0.9,
            'Campylobacter spp.': 0.9,
        }
        if frac < qcdict[expect]:
            raise PydanticCustomError(
                "value_error",
                f"Value error: 'fraction_majority_genus' for '{expect}' must be at least {qcdict[expect]}, got: {frac}",
            )
        return self

    @model_validator(mode='after')
    def check_species_specific_majority_genus(self) -> 'AssemblyQC':
        expect = self.expect_species
        genus = self.majority_genus
        qcdict = {
            'Listeria monocytogenes': ["Listeria"],
            'Salmonella enterica': ["Salmonella"],
            'Escherichia coli': ["Escherichia", "Shigella"],
            'Campylobacter spp.': ["Campylobacter"],
        }
        if genus not in qcdict[expect]:
            raise PydanticCustomError(
                "value_error",
                f"Value error: 'majority_genus' for '{expect}' must be in {qcdict[expect]}, got: {genus}",
            )
        return self


def validate_record(record: dict, model: BaseModel):
    """
    Validate a single record against a data model

    Each entry of the record is checked against the schema. If errors are
    encountered, all o them are formattred in a user-friendly way.
    The output is a tuple whose first element is the check result and the
    second element is the list of formatted error messages.

    Parameters
    ----------
    record : dict
        A dictionnary of key-values generated by json.load.
    validator: dict
        A valid instance of jsonschema.protocols.Validator.

    Returns
    -------
    tuple:
        a tuple (STATUS, Errors) where status is either of PASS or FAIL and
        Errors us a list of all validation errors encountered.
    """
    messages = []
    try:
        m = model.model_validate(record)
        status = "PASS"
    except ValidationError as e:
        status = "FAIL"
        m = {}
        for error in e.errors():
            # Make error pretty
            if len(error['loc']) == 1:
                msg = f"Invalid value in field '{error['loc'][0]}'. "
            elif len(error['loc']) > 1:
                msg = f"Invalid value in field '{error['loc']}'. "
            else:
                msg = ''
            if type(error['input']) is dict:
                msg += f"{error['msg']}"
            else:
                msg += f"'{error['msg']}', got '{error['input']}' instead"
            messages.append(msg)
    return status, messages, m


def main(metrics, json_path, tsv_path):
    # load metrics
    with open(metrics, 'r') as fp:
        records = json.load(fp)

    # Validate each record and register validation error
    validation_status = {}
    for record in records:
        status, errors, parsed = validate_record(record, AssemblyQC)
        validation_status.update(
            {
                record['isolate_id']: {
                    "STATUS": status,
                    "MESSAGES": errors
                }
            }
        )

    # Export JSON
    with open(json_path, 'w') as f:
        json.dump(validation_status, f, indent=4)

    # Create DF and export to tsv
    dict_to_df = {
        k: [v['STATUS'], ';'.join(v['MESSAGES'])]
        for k, v in validation_status.items()
    }
    pd.DataFrame.from_dict(
        dict_to_df,
        orient='index',
        columns=['STATUS', 'MESSAGES']
    ).reset_index(
        names='isolate_id'
    ).to_csv(
        tsv_path,
        sep='\t',
        header=True,
        index=False
    )


if __name__ == '__main__':
    main(
        snakemake.input['metrics'],
        snakemake.output['json'],
        snakemake.output['tsv']
    )
