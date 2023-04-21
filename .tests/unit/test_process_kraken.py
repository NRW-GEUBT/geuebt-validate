import os
import sys

from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath
from json import load


sys.path.insert(0, os.path.dirname(__file__))


def test_process_kraken():
    with TemporaryDirectory() as tmpdir:
        workdir = os.path.join(Path(tmpdir), "workdir")
        data_path = PurePosixPath(".tests/unit/process_kraken/data")
        expected_path = PurePosixPath(".tests/unit/process_kraken/expected")
        script_path = PurePosixPath(".tests/../workflow/scripts/process_kraken.py")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(script_path, workdir)

        # dbg
        print("result.json", file=sys.stderr)

        # run function
        sys.path.insert(0, workdir)
        from process_kraken import main
        main(
            kraken=os.path.join(workdir, 'report.kraken'), 
            taxdump=os.path.join(workdir), 
            json_path=os.path.join(workdir, 'result.json'), 
        )

        # Compare resulting jsons as dict
        with open(os.path.join(workdir, 'result.json'), 'r') as res, open(os.path.join(expected_path, 'result.json'), 'r') as expect:
            res_dict = load(res)
            exp_dict = load(expect)
            assert res_dict == exp_dict
