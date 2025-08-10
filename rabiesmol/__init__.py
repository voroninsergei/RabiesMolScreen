"""rabiesmol package.

Lightweight top-level package to keep imports fast.
CLI is available via `rabiesmol.cli:app`.
Programmatic APIs live in submodules (e.g., `rabiesmol.docking_ext`, `rabiesmol.docking`).

This file intentionally avoids importing heavy modules at import time.
"""
from __future__ import annotations

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("rabiesmol")
except PackageNotFoundError:  # local / editable installs
    __version__ = "0.0.0"

__all__ = ["__version__", "dock"]

def dock(ligands, receptor_paths, threads: int = 1, seed: int | None = None, restarts: int = 1):
    """Convenience wrapper that lazy-imports the heavy docking backend."
    This avoids importing RDKit/Plotly/etc. on `import rabiesmol`.
    """
    from pathlib import Path
    from typing import Iterable, List, Dict, Tuple
    # Lazy import here:
    from .docking_ext import batch_dock, DockParams  # type: ignore

    params = DockParams(receptor_paths=list(receptor_paths), seed=seed, restarts=restarts)
    return batch_dock(list(ligands), params, out_dir=Path("work/api"), threads=threads, cache_dir=Path("work/cache"))
