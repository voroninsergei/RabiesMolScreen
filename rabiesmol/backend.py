
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import cpu_count
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Dict, Tuple
import hashlib, json

@dataclass
class Limits:
    threads: int = max(1, (cpu_count() or 2) - 1)
    max_mem_mb: int | None = None

def _hash_key(record: dict) -> str:
(record: dict) -> str:
    """Create a stable short hash for a mapping.

    We sort keys to ensure determinism and json-dump with separators for compactness.
    """
    dump = json.dumps(record, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(dump.encode("utf-8")).hexdigest()

class SimpleCache:
    """Very small JSON-file based cache.

    Values are stored one-per-file under a directory using the key hash.
    Thread-safe for this test environment.
    """
    def __init__(self, root: Path):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()

    def _path(self, key: str) -> Path:
        return self.root / f"{key}.json"

    def contains(self, key: str) -> bool:
        return self._path(key).exists()

    def get(self, key: str):
        p = self._path(key)
        if not p.exists():
            return None
        try:
            with p.open("r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            return None

    def put(self, key: str, value) -> None:
        p = self._path(key)
        tmp = p.with_suffix(".tmp")
        with self._lock:
            with tmp.open("w", encoding="utf-8") as f:
                json.dump(value, f, indent=2)
            tmp.replace(p)

def parallel_map(
    items: Iterable,
    fn: Callable[[Any], Any],
    limits: Limits | None = None,
    cache: SimpleCache | None = None,
    key_fn: Callable[[Any], str] | None = None,
) -> Tuple[list, list[Tuple[Any, str]]]:
    """Map *fn* across *items* using threads with optional caching.

    Returns (results, failures) where failures contains (item, error_str).
    Order of *results* follows completion order (not input order).
    """
    limits = limits or Limits()
    results: list = []
    failures: list[Tuple[Any, str]] = []
    if not items:
        return results, failures

    with ThreadPoolExecutor(max_workers=max(1, limits.threads)) as ex:
        futs = {}
        for it in items:
            key = key_fn(it) if key_fn else None
            if cache and key and cache.contains(key):
                out = cache.get(key)
                if out is not None:
                    results.append(out)
                    continue
            futs[ex.submit(fn, it)] = (it, key)
        for f in as_completed(futs):
            it, key = futs[f]
            try:
                out = f.result()
                results.append(out)
                if cache and key:
                    cache.put(key, out)
            except Exception as e:
                failures.append((it, str(e)))
    return results, failures
