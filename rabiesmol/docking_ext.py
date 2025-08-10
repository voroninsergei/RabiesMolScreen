
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Tuple, Optional
import hashlib, random

@dataclass
class DockParams:
    receptors: List[Path]
    engine: str = "vina"
    exhaustiveness: int = 8
    seed: Optional[int] = None
    cache_dir: Optional[Path] = None

def _score(ligand: str, receptor: Path, seed: Optional[int]) -> float:
    base = f"{ligand}|{receptor.name}|{seed}".encode("utf-8")
    h = hashlib.sha256(base).hexdigest()
    # Map to a deterministic float ~ N(-7, 1)
    rnd = int(h[:8], 16) / 0xFFFFFFFF
    return -8.0 + (rnd - 0.5) * 2.0

def batch_dock(smiles: Iterable[str], params: DockParams) -> Tuple[List[float], List[str]]:
    """Deterministic fake docking for tests; returns (scores, failures)."""
    scores: List[float] = []
    failures: List[str] = []
    for smi in smiles:
        try:
            sc = sum(_score(smi, r, params.seed) for r in params.receptors) / max(len(params.receptors), 1)
            scores.append(round(sc, 3))
        except Exception as e:
            failures.append(smi)
    return scores, failures
