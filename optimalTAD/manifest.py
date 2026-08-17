import hashlib
import json
import os
import platform
import sys
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path


MANIFEST_FILENAME = "run_manifest.json"
MANIFEST_SCHEMA_VERSION = 1

PACKAGE_DISTRIBUTIONS = {
    "optimalTAD": "optimalTAD",
    "numpy": "numpy",
    "scipy": "scipy",
    "pandas": "pandas",
    "h5py": "h5py",
    "matplotlib": "matplotlib",
    "seaborn": "seaborn",
    "cooler": "cooler",
    "pyBigWig": "pyBigWig",
    "cooltools": "cooltools",
    "bioframe": "bioframe",
    "numba": "numba",
}


def utc_timestamp():
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_file(path, chunk_size=1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as input_file:
        for chunk in iter(lambda: input_file.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _local_input_path(path):
    # Cooler URIs may include a group suffix such as
    # sample.mcool::resolutions/20000. Hash the containing file.
    filesystem_path = str(path).split("::", 1)[0]
    return Path(filesystem_path).expanduser().resolve(strict=True)


def collect_input_hashes(hic_files, chipseq_files):
    hash_cache = {}

    def describe_with_cache(path):
        local_path = _local_input_path(path)
        cache_key = str(local_path)
        if cache_key not in hash_cache:
            hash_cache[cache_key] = {
                "resolved_path": cache_key,
                "size_bytes": local_path.stat().st_size,
                "sha256": sha256_file(local_path),
            }
        return {"path": str(path), **hash_cache[cache_key]}

    return {
        "hic": [describe_with_cache(path) for path in hic_files],
        "chipseq": [describe_with_cache(path) for path in chipseq_files],
    }


def collect_package_versions():
    versions = {}
    for package_name, distribution_name in PACKAGE_DISTRIBUTIONS.items():
        try:
            versions[package_name] = metadata.version(distribution_name)
        except metadata.PackageNotFoundError:
            versions[package_name] = None
    return versions


def configuration_to_dict(cfg):
    return {
        section: dict(cfg.items(section))
        for section in cfg.sections()
    }


def arguments_to_dict(args):
    return {
        key: _json_safe(value)
        for key, value in vars(args).items()
    }


def _json_safe(value):
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if hasattr(value, "tolist"):
        return _json_safe(value.tolist())
    if hasattr(value, "item"):
        return _json_safe(value.item())
    return str(value)


def create_run_manifest(args, cfg, seed):
    created_at = utc_timestamp()
    return {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "status": "running",
        "created_at_utc": created_at,
        "completed_at_utc": None,
        "seed": seed,
        "parameters": {
            "arguments": arguments_to_dict(args),
            "configuration": configuration_to_dict(cfg),
        },
        "inputs": collect_input_hashes(args.hic, args.chipseq),
        "package_versions": collect_package_versions(),
        "runtime": {
            "python_version": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "platform": platform.platform(),
            "executable": sys.executable,
            "command": [str(argument) for argument in sys.argv],
            "working_directory": os.getcwd(),
        },
    }


def mark_completed(run_manifest):
    run_manifest["status"] = "completed"
    run_manifest["completed_at_utc"] = utc_timestamp()


def write_run_manifest(output_directory, run_manifest):
    output_path = Path(output_directory)
    output_path.mkdir(parents=True, exist_ok=True)
    manifest_path = output_path / MANIFEST_FILENAME
    temporary_path = output_path / (MANIFEST_FILENAME + ".tmp")

    with temporary_path.open("w", encoding="utf-8") as manifest_file:
        json.dump(run_manifest, manifest_file, indent=2, sort_keys=True)
        manifest_file.write("\n")

    os.replace(temporary_path, manifest_path)
    return manifest_path
