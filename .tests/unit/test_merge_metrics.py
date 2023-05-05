import os
import sys

from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath
from json import load


sys.path.insert(0, os.path.dirname(__file__))


def test_merge_metrics():
    with TemporaryDirectory() as tmpdir:
        workdir = os.path.join(Path(tmpdir), "workdir")
        data_path = PurePosixPath(".tests/unit/merge_metrics/data")
        expected_path = PurePosixPath(".tests/unit/merge_metrics/expected")
        script_path = PurePosixPath(".tests/../workflow/scripts/merge_metrics.py")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(script_path, workdir)

        # run function
        sys.path.insert(0, workdir)
        from merge_metrics import main
        main(
            quast=os.path.join(workdir, 'transposed_report.tsv'),
            kraken=os.path.join(workdir, '2022-0232977-01.kraken.json'),
            metadata=os.path.join(workdir, 'metadata.json'),
            jsonpath=os.path.join(workdir, 'result.json'),
            busco=os.path.join(workdir, 'short_summary.json'),
            isolate_id='2022-0232977-01'
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
