#!/usr/bin/env python
from pathlib import Path
from rabiesmol.docking import run_docking
from loguru import logger

proteins = Path("data/prepared_proteins")
ligands = Path("data/prepared_ligands")
out_dir = Path("data/docking")

for protein in proteins.glob("*.pdbqt"):
    for ligand in ligands.glob("*.pdbqt"):
        run_docking(protein, ligand, out_dir)
logger.info("Docking complete.")
