
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Dict, Tuple
import hashlib, json, threading

@dataclass
class Limits:
    threads: int = max(1, (cpu_count() or 2) - 1)
    max_mem_mb: int | None = None
    pool: str = "thread"  # 'thread' or 'process'

class SimpleCache:
    def __init__(self, root: Path):
        self.root = Path(root); self.root.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()

    def _path(self, key: str) -> Path:
        return self.root / f"{key}.json"

    def get(self, key: str):
        p = self._path(key)
        if not p.exists(): return None
        try:
            with p.open("r", encoding="utf-8") as f:
                import json; return json.load(f)
        except Exception:
            return None

    def put(self, key: str, value) -> None:
        p = self._path(key)
        with self._lock:
            with p.with_suffix(".tmp").open("w", encoding="utf-8") as f:
                import json; json.dump(value, f, indent=2)
            p.with_suffix(".tmp").replace(p)

def _hash_key(record: dict) -> str:
    dump = json.dumps(record, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(dump.encode("utf-8")).hexdigest()

def parallel_map(
    items: Iterable,
    fn: Callable[[Any], Any],
    limits: Limits | None = None,
    cache: SimpleCache | None = None,
    key_fn: Callable[[Any], dict] | None = None,
) -> Tuple[list, list]:
    limits = limits or Limits()
    results: list = []
    failures: list[tuple[Any, str]] = []
    futs: dict = {}
    key_fn = key_fn or (lambda it: it if isinstance(it, dict) else {"it": str(it)})
    Exec = ProcessPoolExecutor if (limits.pool == "process") else ThreadPoolExecutor
    with Exec(max_workers=limits.threads) as ex:
        for it in items:
            key = _hash_key(key_fn(it)) if cache else None
            if cache and key:
                cached = cache.get(key)
                if cached is not None:
                    results.append(cached)
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
