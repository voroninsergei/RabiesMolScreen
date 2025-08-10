
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
    logger.debug(f"Running: {'{'}' '.join(cmd){'}'}")
        subprocess.run(cmd, check=True)
    logger.info(f"Protonated {input_file} at pH {ph} -> {output_file}")

def prepare_proteins(input_dir: Path, output_dir: Path, ph: float = 7.4, keep_waters: List[str] | None = None):
    """Prepare receptors to PDBQT via external toolchain.

    In tests, subprocess is mocked to just create the outputs.
    """
    keep_waters = keep_waters or []
    ensure_dir(output_dir)
    for pdb in Path(input_dir).glob("*.pdb"):
        protonated = output_dir / pdb.name.replace(".pdb", "_protonated.pdb")
        prepared = output_dir / pdb.name.replace(".pdb", "_prepared.pdbqt")
        if prepared.exists():
            logger.info(f"Cached prepared protein: {prepared}")
            continue
        # Protonate and then call a hypothetical prepare_receptor tool
        protonate_structure(pdb, protonated, ph=ph)
        cmd = ["prepare_receptor", "-r", str(protonated), "-o", str(prepared)]
        if keep_waters:
            cmd += ["--keep-waters", ",".join(keep_waters)]
        logger.debug(f"Running: {'{'}' '.join(cmd){'}'}")
        subprocess.run(cmd, check=True)
        logger.info(f"Prepared receptor: {prepared}")

def prepare_ligands(input_dir: Path, output_dir: Path, ph: float = 7.4):
    ensure_dir(output_dir)
    for mol_file in Path(input_dir).glob("*.sdf"):
        protonated_file = output_dir / mol_file.name.replace(".sdf", "_protonated.sdf")
        prepared_file = output_dir / mol_file.name.replace(".sdf", "_prepared.pdbqt")

        if prepared_file.exists():
            logger.info(f"Cached prepared ligand: {prepared_file}")
            continue

        protonate_structure(mol_file, protonated_file, ph=ph)
        logger.debug("Running: prepare_ligand")
        subprocess.run(["prepare_ligand", "-l", str(protonated_file), "-o", str(prepared_file)], check=True)
        logger.info(f"Prepared ligand: {prepared_file}")
