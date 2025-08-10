
from __future__ import annotations
import time, json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Callable, Any, Dict

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

def emit_event(work_dir: Path, kind: str, payload: Dict[str, Any]) -> None:
    """Append a structured event to run_summary.jsonl."""
    work_dir.mkdir(parents=True, exist_ok=True)
    rec = {"ts": time.time(), "kind": kind, **payload}
    with (work_dir / "run_summary.jsonl").open("a", encoding="utf-8") as f:
        f.write(json.dumps(rec) + "\n")
