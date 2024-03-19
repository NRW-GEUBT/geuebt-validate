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
        data_path = PurePosixPath(".tests/unit/validate_metadata/data")
        expected_path = PurePosixPath(".tests/unit/validate_metadata/expected")
        script_path = PurePosixPath(".tests/../workflow/scripts/validate_metadata.py")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(script_path, workdir)

        # run function
        sys.path.insert(0, workdir)
        from validate_metadata import main
        main(
            metadata=os.path.join(workdir, 'metadata_table_test.tsv'),
            json_path=os.path.join(workdir, 'result.json'),
            tsv_path=os.path.join(workdir, 'result.tsv'),
            metadata_json=os.path.join(workdir, 'metadata.json')
        )

        # Compare resulting jsons as dict
        with open(
            os.path.join(workdir, 'result.json'), 'r'
        ) as res, open(
            os.path.join(expected_path, 'metadata_status.json'), 'r'
        ) as expect:
            res_dict = load(res)
            exp_dict = load(expect)
            assert res_dict == exp_dict
