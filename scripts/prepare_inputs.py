#!/usr/bin/env python
from pathlib import Path
from rabiesmol.prepare import prepare_proteins, prepare_ligands
from loguru import logger

proteins = Path("data/proteins")
ligands = Path("data/ligands")
prepare_proteins(proteins, Path("data/prepared_proteins"))
prepare_ligands(ligands, Path("data/prepared_ligands"))
logger.info("Preparation complete.")
