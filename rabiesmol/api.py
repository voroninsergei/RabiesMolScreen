from __future__ import annotations
from pathlib import Path
from typing import List
from .docking_ext import batch_dock, DockParams

def dock(ligands: List[Path], proteins: List[Path], out_dir: Path, **kwargs):
    params = DockParams(
        receptor_paths=proteins,
        seed=kwargs.get("seed"),
        restarts=kwargs.get("restarts", 1),
    )
    return batch_dock(
        ligands,
        params,
        out_dir=out_dir,
        threads=kwargs.get("threads", 1),
        cache_dir=kwargs.get("cache_dir"),
    )
