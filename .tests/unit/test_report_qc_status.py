import os
import sys

from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath
from json import load


sys.path.insert(0, os.path.dirname(__file__))


def test_report_qc_status():
    with TemporaryDirectory() as tmpdir:
        workdir = os.path.join(Path(tmpdir), "workdir")
        data_path = PurePosixPath(".tests/unit/report_qc_status/data")
        expected_path = PurePosixPath(".tests/unit/report_qc_status/expected")
        script_path = PurePosixPath(".tests/../workflow/scripts/report_qc_status.py")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(script_path, workdir)

        # run function
        sys.path.insert(0, workdir)
        from report_qc_status import main
        main(
            metadata_status=os.path.join(workdir, 'metadata_status.json'),
            checksum_status=os.path.join(workdir, 'checksum_status.json'),
            assembly_status=os.path.join(workdir, 'assemblies_status.json'),
            json_path=os.path.join(workdir, 'result.json'),
        )

        # Compare resulting jsons as dict
        with open(
            os.path.join(workdir, 'result.json'), 'r'
        ) as res, open(
            os.path.join(expected_path, 'result.json'), 'r'
        ) as expect:
            res_dict = load(res)
            exp_dict = load(expect)
            assert res_dict == exp_dict
