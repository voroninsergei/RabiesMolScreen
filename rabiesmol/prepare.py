from __future__ import annotations
from pathlib import Path
from typing import List
from loguru import logger
import shutil
import subprocess
from .proc import run
from rabiesmol.io import ensure_dir

# --- Ligands ---
def protonate_structure(input_file: Path, output_file: Path, ph: float = 7.4):
    cmd = ["obabel", str(input_file), "-O", str(output_file), "--pH", str(ph)]
    logger.debug(f"Running: {' '.join(cmd)}")
    run(cmd)
    logger.info(f"Protonated {input_file} -> {output_file}")

def _smi_to_sdf(smi_file: Path, out_sdf: Path) -> Path:
    cmd = ["obabel", "-ismi", str(smi_file), "-O", str(out_sdf)]
    logger.debug(f"Running: {' '.join(cmd)}")
    try:
        run(cmd)
    except FileNotFoundError as e:
        logger.error("OpenBabel (`obabel`) is required to convert SMILES. Please install OpenBabel or provide SDF.")
        raise
    return out_sdf

def prepare_ligands(input_dir: Path, output_dir: Path, ph: float = 7.4):
    ensure_dir(output_dir)
    for mol_file in list(Path(input_dir).glob("*.sdf")) + list(Path(input_dir).glob("*.smi")):
        protonated_file = output_dir / mol_file.name.replace(mol_file.suffix, "_protonated.sdf")
        prepared_file = output_dir / mol_file.name.replace(mol_file.suffix, "_prepared.pdbqt")

        if prepared_file.exists():
            logger.info(f"Cached prepared ligand: {prepared_file}")
            continue

        source_file = mol_file
        if mol_file.suffix.lower() == '.smi':
            sdf_intermediate = output_dir / mol_file.name.replace('.smi', '.sdf')
            source_file = _smi_to_sdf(mol_file, sdf_intermediate)

        protonate_structure(source_file, protonated_file, ph=ph)
        logger.debug("Running: prepare_ligand")
        run(["prepare_ligand", "--ligand", str(protonated_file), "--out", str(prepared_file)])
        logger.info(f"Prepared ligand: {prepared_file}")

# --- Proteins ---
def clean_protein(protein_file: Path, out_file: Path) -> Path:
    """Clean protein (remove waters/ligands) and ensure PDBQT output.
    Prefer AutoDockTools `prepare_receptor`, otherwise try OpenBabel, otherwise copy.
    """
    out_file = Path(out_file)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    # Try AutoDockTools (if available)
    try:
        subprocess.run(["prepare_receptor", "-r", str(protein_file), "-o", str(out_file)], check=True)
        logger.info(f"Prepared receptor via AutoDockTools: {out_file}")
        return out_file
    except FileNotFoundError:
        logger.debug("prepare_receptor not found; falling back to OpenBabel/copy")  # continue

    # Try OpenBabel to convert to pdbqt
    try:
        cmd = ["obabel", str(protein_file), "-O", str(out_file), "-xr", "-xh"]
        logger.debug(f"Running: {' '.join(cmd)}")
        run(cmd)
        logger.info(f"Prepared receptor via OpenBabel: {out_file}")
        return out_file
    except FileNotFoundError:
        logger.warning("OpenBabel not found; copying input file as-is.")

    shutil.copyfile(protein_file, out_file)
    return out_file

def prepare_proteins(input_dir: Path, output_dir: Path) -> None:
    ensure_dir(output_dir)
    proteins = list(Path(input_dir).glob("*.pdb")) + list(Path(input_dir).glob("*.pdbqt"))
    if not proteins:
        logger.warning(f"No proteins found in {input_dir} (expected *.pdb or *.pdbqt)")
    for p in proteins:
        out = output_dir / p.name.replace(p.suffix, "_prepared.pdbqt")
        if out.exists():
            logger.info(f"Cached prepared protein: {out}")
            continue
        clean_protein(p, out)
