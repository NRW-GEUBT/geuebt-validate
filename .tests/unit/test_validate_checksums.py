import os
import sys

from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath
from json import load


sys.path.insert(0, os.path.dirname(__file__))


def test_validate_checksums():
    with TemporaryDirectory() as tmpdir:
        workdir = os.path.join(Path(tmpdir), "workdir")
        data_path = PurePosixPath(".tests/unit/validate_checksums/data")
        expected_path = PurePosixPath(".tests/unit/validate_checksums/expected")
        script_path = PurePosixPath(".tests/../workflow/scripts/validate_checksums.py")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(script_path, workdir)

        # dbg
        print("fasta_checksums.json", file=sys.stderr)

        # run function
        sys.path.insert(0, workdir)
        from validate_checksums import main
        main(
            metadata=os.path.join(workdir, 'metadata_table_test.tsv'), 
            metadata_qc=os.path.join(workdir, 'metadata_status.json'), 
            fasta_dir=os.path.join(workdir, 'fastas'), 
            json_path=os.path.join(workdir, 'result.json'), 
            tsv_path=os.path.join(workdir, 'result.tsv')
        )

        # Compare resulting jsons as dict
        with open(os.path.join(workdir, 'result.json'), 'r') as res, open(os.path.join(expected_path, 'fasta_checksums.json'), 'r') as expect:
            res_dict = load(res)
            exp_dict = load(expect)
            assert res_dict == exp_dict
