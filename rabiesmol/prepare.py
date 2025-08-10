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
        "--pH", str(ph),
    ]
    logger.debug(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    logger.info(f"Protonated {input_file} -> {output_file}")

def _smi_to_sdf(smi_file: Path, out_sdf: Path) -> Path:
    """Convert .smi to .sdf using OpenBabel if available."""
    cmd = ["obabel", "-ismi", str(smi_file), "-O", str(out_sdf)]
    logger.debug(f"Running: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
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
        subprocess.run(["prepare_ligand", "-l", str(protonated_file), "-o", str(prepared_file)], check=True)
        logger.info(f"Prepared ligand: {prepared_file}")
