"""rabiesmol package.

Lightweight top-level package to keep imports fast.
CLI is available via ``rabiesmol.cli:app``.
Programmatic APIs live in submodules (e.g., ``rabiesmol.docking_ext``, ``rabiesmol.docking``).

This file intentionally avoids importing heavy modules at import time.
"""
from __future__ import annotations

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("rabiesmol")  # type: ignore[assignment]
except PackageNotFoundError:  # local / editable installs
    __version__ = "0.0.0"

__all__ = ["__version__", "dock"]

def dock(ligands, receptor_paths, threads: int = 1, seed: int | None = None, restarts: int = 1):
    """Convenience wrapper that lazy-imports the heavy docking backend.

    This avoids importing RDKit/Plotly/etc. on ``import rabiesmol``.

    Parameters
    ----------
    ligands:
        Iterable of ligand SMILES strings to dock.
    receptor_paths:
        Iterable of Paths to prepared receptor files.
    threads:
        Number of parallel threads for docking.
    seed:
        Optional random seed for deterministic docking.
    restarts:
        Number of restarts for the docking algorithm.
    """
    from pathlib import Path
    from typing import Iterable, List, Dict, Tuple
    # Lazy import here:
    from .docking_ext import batch_dock, DockParams  # type: ignore

    params = DockParams(receptors=list(receptor_paths), seed=seed, restarts=restarts)  # type: ignore[arg-type]
    return batch_dock(list(ligands), params)