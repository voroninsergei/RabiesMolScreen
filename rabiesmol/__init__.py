# RabiesMol package


from pathlib import Path
from typing import Iterable, List, Dict, Tuple
from .docking_ext import batch_dock, DockParams, Constraint

def dock(
    ligands: Iterable[str],
    receptor_paths: Iterable[Path],
    threads: int = 1,
    seed: int | None = None,
    restarts: int = 1,
) -> Tuple[List[Dict], List[Dict]]:
    """High-level Python API for docking without CLI."""
    params = DockParams(receptor_paths=list(receptor_paths), seed=seed, restarts=restarts)
    return batch_dock(list(ligands), params, out_dir=Path("work/api"), threads=threads, cache_dir=Path("work/cache"))
