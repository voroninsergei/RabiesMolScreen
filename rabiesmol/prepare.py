from pathlib import Path
from typing import List
from loguru import logger
import subprocess
from rabiesmol.io import ensure_dir

# pH-dependent protonation with OpenBabel
def protonate_structure(input_file: Path, output_file: Path, ph: float = 7.4):
    cmd = [
        "obabel",
        str(input_file),
        "-O", str(output_file),
        "--pH", str(ph)
    ]
    subprocess.run(cmd, check=True)
    logger.info(f"Protonated {input_file} at pH {ph} -> {output_file}")

# Remove waters and ions, keep selected waters
def clean_protein(pdb_in: Path, pdb_out: Path, keep_waters: List[str] = None):
    keep_waters = keep_waters or []
    with open(pdb_in, "r") as f:
        lines = f.readlines()
    filtered = []
    for line in lines:
        if line.startswith("HETATM") and "HOH" in line and line[17:20] not in keep_waters:
            continue
        if line.startswith("HETATM") and any(ion in line for ion in ["NA", "CL", "MG", "CA", "K"]):
            continue
        filtered.append(line)
    with open(pdb_out, "w") as f:
        f.writelines(filtered)
    logger.info(f"Cleaned protein: {pdb_in} -> {pdb_out}")

def prepare_proteins(input_dir: Path, output_dir: Path, ph: float = 7.4, keep_waters: List[str] = None):
    ensure_dir(output_dir)
    for pdb_file in Path(input_dir).glob("*.pdb"):
        cleaned_file = output_dir / pdb_file.name.replace(".pdb", "_clean.pdb")
        protonated_file = output_dir / pdb_file.name.replace(".pdb", "_prepared.pdbqt")

        if protonated_file.exists():
            logger.info(f"Cached prepared protein: {protonated_file}")
            continue

        clean_protein(pdb_file, cleaned_file, keep_waters=keep_waters)
        protonate_structure(cleaned_file, cleaned_file, ph=ph)
        subprocess.run(["prepare_receptor", "-r", str(cleaned_file), "-o", str(protonated_file)], check=True)
        logger.info(f"Prepared protein: {protonated_file}")

def prepare_ligands(input_dir: Path, output_dir: Path, ph: float = 7.4):
    ensure_dir(output_dir)
    for mol_file in Path(input_dir).glob("*.sdf"):
        protonated_file = output_dir / mol_file.name.replace(".sdf", "_protonated.sdf")
        prepared_file = output_dir / mol_file.name.replace(".sdf", "_prepared.pdbqt")

        if prepared_file.exists():
            logger.info(f"Cached prepared ligand: {prepared_file}")
            continue

        protonate_structure(mol_file, protonated_file, ph=ph)
        subprocess.run(["prepare_ligand", "-l", str(protonated_file), "-o", str(prepared_file)], check=True)
        logger.info(f"Prepared ligand: {prepared_file}")
