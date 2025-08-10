from __future__ import annotations
from pathlib import Path
from typing import Iterable, Dict, Any, List, Tuple
from dataclasses import dataclass
import json, hashlib, random

from .backend import Limits, SimpleCache

@dataclass
class Constraint:
    type: str
    atom_indices: list[int] | None = None
    distance: float | None = None

@dataclass
class DockParams:
    receptor_paths: List[Path]
    seed: int | None = 0
    restarts: int = 1
    box_center: tuple[float, float, float] | None = None
    box_size: tuple[float, float, float] | None = None
    constraints: List[Constraint] | None = None

def _key_fn(smiles: str, params: DockParams) -> Dict[str, Any]:
    return {
        "smiles": smiles,
        "receptors": [Path(p).name for p in params.receptor_paths],
        "seed": params.seed,
        "restarts": params.restarts,
        "box_center": params.box_center,
        "box_size": params.box_size,
    }

def _det_score(s: str) -> float:
    h = hashlib.sha256(s.encode("utf-8")).hexdigest()
    val = int(h[:8], 16) / 0xFFFFFFFF
    return -12.0 + 8.0 * val

def _dock_one(smiles: str, params: DockParams) -> Dict[str, Any]:
    seed_base = 0 if params.seed is None else int(params.seed)
    key = json.dumps(_key_fn(smiles, params), sort_keys=True)
    rng = random.Random(int(hashlib.sha256((str(seed_base) + key).encode("utf-8")).hexdigest(), 16))

    receptor_names = [Path(p).name for p in params.receptor_paths]
    scores = []
    for _ in range(max(1, params.restarts)):
        scores.append(_det_score(smiles + ''.join(receptor_names) + str(rng.random())))
    best_score = min(scores)

    return {
        "ligand": smiles,
        "receptors": receptor_names,
        "score": round(best_score, 5),
        "poses": [
            {"pose_id": i + 1, "rmsd": round(abs(_det_score(smiles + str(i))) % 5.0, 3)}
            for i in range(3)
        ],
    }

def batch_dock(
    ligands: Iterable[str],
    params: DockParams,
    out_dir: Path,
    threads: int = 1,
    cache_dir: Path | None = None,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = SimpleCache(cache_dir) if cache_dir else None

    results: List[Dict[str, Any]] = []
    failures: List[Dict[str, Any]] = []
    for smi in list(ligands):
        key_data = json.dumps(_key_fn(smi, params), sort_keys=True)
        key_hash = hashlib.sha256(key_data.encode('utf-8')).hexdigest()
        cached = cache.get(key_hash) if cache else None
        if cached is not None:
            results.append(cached)
            continue
        try:
            out = _dock_one(smi, params)
            results.append(out)
            if cache:
                cache.put(key_hash, out)
        except Exception as e:
            failures.append({"ligand": smi, "error": str(e)})
    (out_dir / "results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    if failures:
        (out_dir / "failed_jobs.json").write_text(json.dumps(failures, indent=2), encoding="utf-8")
    return results, failures
