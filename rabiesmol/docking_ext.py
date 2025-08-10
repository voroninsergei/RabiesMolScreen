from __future__ import annotations
from pathlib import Path
from typing import Iterable, Dict, Any, List, Tuple
from dataclasses import dataclass
import time, random, json, hashlib

from .backend import parallel_map, Limits, SimpleCache

@dataclass
class Constraint:
    type: str
    atom_indices: list[int] | None = None
    distance: float | None = None

@dataclass
class DockParams:
    receptor_paths: List[Path]
    seed: int = 0
    restarts: int = 1
    box_center: tuple[float, float, float] | None = None
    box_size: tuple[float, float, float] | None = None

def _key_fn(smiles: str, params: DockParams) -> str:
    """Cache key for a docking job."""
    payload = {
        "smi": smiles,
        "rec_sha": [Path(p).name for p in params.receptor_paths],
        "seed": params.seed,
        "restarts": params.restarts,
        "box_center": params.box_center,
        "box_size": params.box_size,
    }
    dump = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(dump.encode("utf-8")).hexdigest()

def _dock_one(smiles: str, params: DockParams) -> Dict[str, Any]:
    """Deterministic fake docking used for tests.

    We seed python's RNG using a hash over input, then generate plausible scores.
    """
    # Derive deterministic seed from global seed + ligand + receptors list
    h = int(hashlib.sha256((smiles + "|" + "|".join(sorted([Path(p).name for p in params.receptor_paths])) + f"|{params.seed}|{params.restarts}").encode("utf-8")).hexdigest(), 16)
    rng = random.Random(h)
    # Simulate some time cost
    time.sleep(0.001)
    vina = -3.0 - rng.random() * 9.0
    smina = vina + (rng.random() - 0.5) * 0.5
    consensus = (vina + smina) / 2.0
    return {"smiles": smiles, "vina_score": round(vina, 3), "smina_score": round(smina, 3), "consensus": round(consensus, 3)}

def batch_dock(ligands: Iterable[str], params: DockParams, out_dir: Path, threads: int = 2, cache_dir: Path | None = None) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    limits = Limits(threads=threads)
    cache = SimpleCache(cache_dir) if cache_dir else None

    results, failures = parallel_map(
        ligands,
        lambda smi: _dock_one(smi, params),
        limits=limits,
        cache=cache,
        key_fn=lambda smi: _key_fn(smi, params),
    )
    # Record failed jobs
    failed_registry: List[Dict[str, Any]] = []
    if failures:
        for it, err in failures:
            failed_registry.append({"ligand": it, "error": err})
        (out_dir / "failed_jobs.json").write_text(json.dumps(failed_registry, indent=2), encoding="utf-8")
    return results, failed_registry
