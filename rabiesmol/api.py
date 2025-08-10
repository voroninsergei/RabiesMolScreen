
from __future__ import annotations
from pathlib import Path
from typing import List, Optional
from .docking_ext import dock_batch, dock_with_constraints

def dock(ligands: List[Path], proteins: List[Path], out_dir: Path, **kwargs):
    return dock_batch(proteins, ligands, out_dir, **kwargs)
