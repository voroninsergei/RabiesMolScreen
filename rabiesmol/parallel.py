
from __future__ import annotations
import concurrent.futures, hashlib, json
from pathlib import Path
from typing import Callable, Any

def _hash_inputs(data: dict) -> str:
    s = json.dumps(data, sort_keys=True)
    return hashlib.sha256(s.encode()).hexdigest()

class ParallelBackend:
    def __init__(self, max_workers: int = 4, cache_dir: Path = Path(".cache_jobs")):
        self.max_workers = max_workers
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def run(self, tasks: list[dict], func: Callable[[dict], Any]) -> list[Any]:
        results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_workers) as ex:
            fut_to_task = {}
            for task in tasks:
                key = _hash_inputs(task)
                cache_file = self.cache_dir / f"{key}.json"
                if cache_file.exists():
                    results.append(json.loads(cache_file.read_text()))
                else:
                    fut = ex.submit(func, task)
                    fut_to_task[fut] = (task, cache_file)
            for fut in concurrent.futures.as_completed(fut_to_task):
                task, cache_file = fut_to_task[fut]
                try:
                    res = fut.result()
                    cache_file.write_text(json.dumps(res))
                    results.append(res)
                except Exception as e:
                    results.append({"error": str(e), "task": task})
        return results
