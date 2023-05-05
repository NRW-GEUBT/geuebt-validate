import os
import sys

from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath
from json import load


sys.path.insert(0, os.path.dirname(__file__))


def test_create_isolate_sheet():
    with TemporaryDirectory() as tmpdir:
        workdir = os.path.join(Path(tmpdir), "workdir")
        data_path = PurePosixPath(".tests/unit/create_isolate_sheet/data")
        expected_path = PurePosixPath(".tests/unit/create_isolate_sheet/expected")
        script_path = PurePosixPath(".tests/../workflow/scripts/create_isolate_sheet.py")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(script_path, workdir)

        # run function
        sys.path.insert(0, workdir)
        from create_isolate_sheet import main
        main(
            mlst=os.path.join(workdir, '2022-0232977-01.mlst.json'),
            metadata=os.path.join(workdir, 'metadata.json'),
            assembly_qc=os.path.join(workdir, '2022-0232977-01.qc_metrics.json'),
            isolate_id='2022-0232977-01',
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
