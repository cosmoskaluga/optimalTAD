import configparser
import hashlib
import json
from argparse import Namespace

from optimalTAD import manifest


def test_run_manifest_records_parameters_versions_and_input_hashes(tmp_path):
    hic_path = tmp_path / "sample.cool"
    chipseq_path = tmp_path / "sample.bedgraph"
    hic_path.write_bytes(b"hic-data")
    chipseq_path.write_bytes(b"chipseq-data")

    args = Namespace(
        hic=[str(hic_path) + "::resolutions/20000"],
        chipseq=[str(chipseq_path)],
        output=str(tmp_path / "output"),
        seed=17,
        np=2,
    )
    cfg = configparser.ConfigParser()
    cfg.read_dict({"run": {"seed": "17", "resolution": "20000"}})

    run_manifest = manifest.create_run_manifest(args, cfg, seed=17)
    manifest_path = manifest.write_run_manifest(args.output, run_manifest)

    saved_manifest = json.loads(manifest_path.read_text())
    assert saved_manifest["seed"] == 17
    assert saved_manifest["status"] == "running"
    assert saved_manifest["parameters"]["arguments"]["np"] == 2
    assert saved_manifest["parameters"]["configuration"]["run"]["resolution"] == "20000"
    assert saved_manifest["inputs"]["hic"][0]["sha256"] == hashlib.sha256(
        b"hic-data"
    ).hexdigest()
    assert saved_manifest["inputs"]["chipseq"][0]["sha256"] == hashlib.sha256(
        b"chipseq-data"
    ).hexdigest()
    assert "numpy" in saved_manifest["package_versions"]
    assert saved_manifest["runtime"]["python_version"]


def test_mark_completed_sets_status_and_timestamp():
    run_manifest = {"status": "running", "completed_at_utc": None}
    manifest.mark_completed(run_manifest)

    assert run_manifest["status"] == "completed"
    assert run_manifest["completed_at_utc"].endswith("Z")
