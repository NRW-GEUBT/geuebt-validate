import os
import sys

from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath
from json import load


sys.path.insert(0, os.path.dirname(__file__))


def test_validate_metadata():
    with TemporaryDirectory() as tmpdir:
        workdir = os.path.join(Path(tmpdir), "workdir")
        data_path = PurePosixPath(".tests/unit/validate_assemblies/data")
        expected_path = PurePosixPath(".tests/unit/validate_assemblies/expected")
        script_path = PurePosixPath(".tests/../workflow/scripts/validate_assemblies.py")
        schema_path = PurePosixPath(".tests/../workflow/schema/assembly_qc.schema.json")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(schema_path, workdir)
        shutil.copy(script_path, workdir)

        # dbg
        print("assemblies_status.json", file=sys.stderr)

        # run function
        sys.path.insert(0, workdir)
        from validate_assemblies import main
        main(
            schema=os.path.join(workdir, 'assembly_qc.schema.json'), 
            metrics=os.path.join(workdir, 'metrics_test.json'), 
            json_path=os.path.join(workdir, 'result.json'), 
            tsv_path=os.path.join(workdir, 'result.tsv'),
            # metadata_json=os.path.join(workdir, 'metadata.json')
        )

        # Compare resulting jsons as dict
        with open(os.path.join(workdir, 'result.json'), 'r') as res, open(os.path.join(expected_path, 'assemblies_status.json'), 'r') as expect:
            res_dict = load(res)
            exp_dict = load(expect)
            assert res_dict == exp_dict