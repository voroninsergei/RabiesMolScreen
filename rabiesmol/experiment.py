
from __future__ import annotations
from pathlib import Path
from typing import Mapping, Optional, Iterable, Dict, Any
import sys, hashlib, json, platform, datetime, subprocess, os
import yaml

from .logging_config import get_logger

logger = get_logger(__name__)

def _sha256_of_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            h.update(chunk)
    return h.hexdigest()

def _pip_freeze() -> str:
    try:
        out = subprocess.check_output([sys.executable, "-m", "pip", "freeze"], text=True)
        return out.strip()
    except Exception as e:
        logger.warning("pip freeze failed: %s", e)
        return ""

def _git_info(dir_: Path) -> Dict[str, Any]:
    try:
        rev = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=str(dir_), text=True).strip()
        status = subprocess.check_output(["git", "status", "--porcelain"], cwd=str(dir_), text=True).strip()
        return {"commit": rev, "dirty": bool(status)}
    except Exception:
        return {}

def _load_yaml(path: Optional[Path]) -> Dict[str, Any]:
    if path and path.exists():
        with path.open("r", encoding="utf-8") as f:
            return yaml.safe_load(f) or {}
    return {}

def snapshot_experiment(
    out_dir: str,
    inputs: Optional[Iterable[str]] = None,
    config_path: Optional[str] = None,
    extra_params: Optional[Mapping[str, Any]] = None,
    code_dir: Optional[str] = None,
) -> str:
    """Create a reproducibility snapshot (manifest + hashes + env info).

    Parameters
    ----------
    out_dir : str
        Directory to write the snapshot into (creates `repro/` subfolder).
    inputs : Optional[Iterable[str]]
        Files or directories whose contents should be hashed.
    config_path : Optional[str]
        Path to YAML with run-time parameters (copied into snapshot).
    extra_params : Optional[Mapping[str, Any]]
        Any extra parameters to include (CLI args, overrides).
    code_dir : Optional[str]
        For Git metadata collection.

    Returns
    -------
    str
        Path to the created manifest.yaml
    """
    out = Path(out_dir) / "repro"
    out.mkdir(parents=True, exist_ok=True)

    # Hash inputs
    hashed: Dict[str, str] = {}
    if inputs:
        for item in inputs:
            p = Path(item)
            if p.is_dir():
                for fp in p.rglob("*"):
                    if fp.is_file():
                        hashed[str(fp)] = _sha256_of_file(fp)
            elif p.is_file():
                hashed[str(p)] = _sha256_of_file(p)

    # Config
    cfg = _load_yaml(Path(config_path)) if config_path else {}
    if config_path and Path(config_path).exists():
        dst = out / Path(config_path).name
        dst.write_text(Path(config_path).read_text(encoding="utf-8"), encoding="utf-8")

    # Env
    env = {
        "python": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "timestamp_utc": datetime.datetime.utcnow().isoformat() + "Z",
        "pip_freeze": _pip_freeze().splitlines(),
    }
    if code_dir:
        env["git"] = _git_info(Path(code_dir))

    manifest = {
        "inputs_sha256": hashed,
        "config": cfg,
        "extra_params": dict(extra_params) if extra_params else {},
        "environment": env,
    }

    manifest_path = out / "manifest.yaml"
    with manifest_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(manifest, f, sort_keys=False, allow_unicode=True)
    logger.info("Snapshot saved to %s", manifest_path)
    return str(manifest_path)


import platform, subprocess, hashlib, json, os

def write_manifest(output_dir: str, params: dict, inputs: list[str], outputs: list[str]):
    os.makedirs(output_dir, exist_ok=True)
    manifest = {
        "timestamp": datetime.datetime.utcnow().isoformat() + "Z",
        "platform": platform.platform(),
        "python": platform.python_version(),
        "git_hash": subprocess.getoutput("git rev-parse HEAD"),
        "pip_freeze": subprocess.getoutput("pip freeze"),
        "conda_list": subprocess.getoutput("conda list") if shutil.which("conda") else None,
        "params": params,
        "inputs": {p: sha256_of_file(p) for p in inputs if os.path.exists(p)},
        "outputs": {p: sha256_of_file(p) for p in outputs if os.path.exists(p)}
    }
    with open(os.path.join(output_dir, "manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

def sha256_of_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()
