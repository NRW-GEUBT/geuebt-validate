import os
import sys

from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath
from json import load


sys.path.insert(0, os.path.dirname(__file__))


def test_merge_isolate_sheets():
    with TemporaryDirectory() as tmpdir:
        workdir = os.path.join(Path(tmpdir), "workdir")
        data_path = PurePosixPath(".tests/unit/merge_isolate_sheets/data")
        expected_path = PurePosixPath(".tests/unit/merge_isolate_sheets/expected")
        script_path = PurePosixPath(".tests/../workflow/scripts/merge_isolate_sheets.py")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(script_path, workdir)

        # run function
        sys.path.insert(0, workdir)
        from merge_isolate_sheets import main
        main(
            qc_pass=[
                os.path.join(workdir, 'test1.json'),
                os.path.join(workdir, 'test2.json'),
                os.path.join(workdir, 'test3.json')
            ],
            json_path=os.path.join(workdir, 'result.json'),
        )

        # Compare resulting jsons as dict
        with open(
            os.path.join(workdir, 'result.json'), 'r'
        ) as res, open(
            os.path.join(expected_path, 'result.json'), 'r'
        ) as expect:
            res_list = load(res)
            exp_list = load(expect)
            assert len(res_list) == len(exp_list)
            for res_dict, exp_dict in zip(res_list, exp_list):
                assert res_dict == exp_dict
