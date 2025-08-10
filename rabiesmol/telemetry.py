
from __future__ import annotations
import time, json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Callable, Any

@dataclass
class Pulse:
    stage: str
    start_ts: float
    end_ts: float | None = None
    cpu_s: float | None = None
    notes: str | None = None

def pulse(stage: str, work_dir: Path, fn: Callable[[], Any]) -> Any:
    work_dir.mkdir(parents=True, exist_ok=True)
    p = Pulse(stage=stage, start_ts=time.time())
    try:
        out = fn()
        return out
    finally:
        p.end_ts = time.time()
        timeline = work_dir / "timeline.jsonl"
        with timeline.open("a", encoding="utf-8") as f:
            f.write(json.dumps(asdict(p)) + "\n")
