import os
import time
from snakemake.utils import validate


# Pipeline setup --------------------------------------
version = open(os.path.join(workflow.basedir, "..", "VERSION"), "r").read()
pipe_log = os.path.join(os.getcwd(), "PIPELINE_STATUS")


# Validating config ----------------------------------
validate(config, schema="../schema/config.schema.yaml")


# General puprose functions --------------------------
def get_local_time():
    return time.asctime(time.localtime(time.time()))


# Input functions ------------------------------------
def aggregate_qcpass(wildcards, pattern):
    checkpoint_output = checkpoints.copy_fasta_to_wdir.get(**wildcards).output[0]
    ids_map = glob_wildcards(
        os.path.join(
            checkpoint_output, "{isolate_id}.fa"
        )).isolate_id
    return expand(pattern, isolate_id=ids_map)