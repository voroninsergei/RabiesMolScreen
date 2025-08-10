
from __future__ import annotations
from pathlib import Path
from typing import Optional, List, Dict
from .logging_config import get_logger
import pandas as pd
import subprocess

logger = get_logger(__name__)

def _parse_vina_log(log_file: Path) -> float:
    try:
        for line in log_file.read_text(encoding="utf-8", errors="ignore").splitlines():
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit():
                return float(parts[1])
    except Exception as e:
        logger.warning(f"Failed to parse {log_file}: {e}")
    return 0.0

def _run_engine(exe: str, receptor: Path, ligand: Path, out_dir: Path, exhaustiveness: int, seed: Optional[int]) -> float:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_pdbqt = out_dir / f"{ligand.stem}_{exe}.pdbqt"
    log_file = out_dir / f"{ligand.stem}_{exe}.log"
    cmd = [
        exe,
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--out", str(out_pdbqt),
        "--log", str(log_file),
        "--exhaustiveness", str(exhaustiveness),
    ]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    logger.debug(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    return _parse_vina_log(log_file)

def run_vina(receptor: Path, ligand: Path, out_dir: Path, exhaustiveness: int = 8, seed: Optional[int] = None) -> float:
    return _run_engine("vina", receptor, ligand, out_dir, exhaustiveness, seed)

def run_smina(receptor: Path, ligand: Path, out_dir: Path, exhaustiveness: int = 8, seed: Optional[int] = None) -> float:
    return _run_engine("smina", receptor, ligand, out_dir, exhaustiveness, seed)

def consensus_docking(receptor: Path, ligand: Path, out_dir: Path, exhaustiveness: int = 8, seed: Optional[int] = None) -> Dict[str, float]:
    v = run_vina(receptor, ligand, out_dir, exhaustiveness, seed)
    s = run_smina(receptor, ligand, out_dir, exhaustiveness, seed)
    return {"vina": v, "smina": s, "consensus": (v + s) / 2.0}

def run_docking_batch(protein: Path, ligands_dir: Path, out_dir: Path) -> pd.DataFrame:
    ligands = sorted([p for p in Path(ligands_dir).glob('*.sdf')] + [p for p in Path(ligands_dir).glob('*.pdbqt')])
    records = []
    for lig in ligands:
        score = run_vina(protein, lig, out_dir)
        records.append({"protein": protein.name, "ligand": lig.name, "vina_score": score})
    return pd.DataFrame.from_records(records)
