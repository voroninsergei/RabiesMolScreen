from pathlib import Path
from typing import Optional, List
from loguru import logger
import subprocess
import pandas as pd

def run_vina(protein: Path, ligand: Path, out_dir: Path, exhaustiveness: int = 8, seed: int = 42) -> float:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"{ligand.stem}_vs_{protein.stem}_vina.pdbqt"
    log_file = out_dir / f"{ligand.stem}_vs_{protein.stem}_vina.log"

    if out_file.exists():
        logger.info(f"Cached Vina result: {out_file}")
        return parse_vina_score(log_file)

    subprocess.run([
        "vina",
        "--receptor", str(protein),
        "--ligand", str(ligand),
        "--out", str(out_file),
        "--log", str(log_file),
        "--center_x", "0", "--center_y", "0", "--center_z", "0",
        "--size_x", "20", "--size_y", "20", "--size_z", "20",
        "--exhaustiveness", str(exhaustiveness),
        "--seed", str(seed)
    ], check=True)

    return parse_vina_score(log_file)

def run_smina(protein: Path, ligand: Path, out_dir: Path, exhaustiveness: int = 8, seed: int = 42) -> float:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"{ligand.stem}_vs_{protein.stem}_smina.pdbqt"
    log_file = out_dir / f"{ligand.stem}_vs_{protein.stem}_smina.log"

    if out_file.exists():
        logger.info(f"Cached smina result: {out_file}")
        return parse_vina_score(log_file)

    subprocess.run([
        "smina",
        "--receptor", str(protein),
        "--ligand", str(ligand),
        "--out", str(out_file),
        "--log", str(log_file),
        "--center_x", "0", "--center_y", "0", "--center_z", "0",
        "--size_x", "20", "--size_y", "20", "--size_z", "20",
        "--exhaustiveness", str(exhaustiveness),
        "--seed", str(seed)
    ], check=True)

    return parse_vina_score(log_file)

def parse_vina_score(log_file: Path) -> float:
    score = None
    with open(log_file, "r") as f:
        for line in f:
            if line.strip().startswith("1 "):  # First mode
                try:
                    score = float(line.split()[1])
                    break
                except:
                    pass
    if score is None:
        logger.warning(f"No score found in {log_file}")
        return 0.0
    return score

def consensus_docking(protein: Path, ligand: Path, out_dir: Path) -> dict:
    vina_score = run_vina(protein, ligand, out_dir)
    smina_score = run_smina(protein, ligand, out_dir)
    return {
        "vina_score": vina_score,
        "smina_score": smina_score,
        "consensus": (vina_score + smina_score) / 2
    }
