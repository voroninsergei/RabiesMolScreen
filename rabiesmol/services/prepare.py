from __future__ import annotations
import os, pathlib, shutil
from typing import Iterable, List
from ..domain.models import Ligand, Protein
from ..utils.logging import get_logger

log = get_logger(__name__)

def prepare_proteins(proteins_dir: str, out_dir: str) -> None:
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
    # For demo: pass-through copy of files
    for p in pathlib.Path(proteins_dir).glob("*.pdb*"):
        shutil.copy(p, pathlib.Path(out_dir) / p.name)
    log.info("Prepared proteins -> %s", out_dir)

def prepare_ligands(ligands_dir: str, out_dir: str) -> None:
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
    for p in pathlib.Path(ligands_dir).glob("*.pdbqt") | pathlib.Path(ligands_dir).glob("*.smi") | pathlib.Path(ligands_dir).glob("*.mol*") :
        shutil.copy(p, pathlib.Path(out_dir) / p.name)
    log.info("Prepared ligands -> %s", out_dir)
