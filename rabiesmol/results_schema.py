
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import json
from typing import List, Dict

SCHEMA_VERSION: str = "1.0"

# Core columns expected in any results parquet
RESULT_COLUMNS: List[str] = [
    "id", "smiles", "vina_score",
    "receptor", "engine", "exhaustiveness",
    "seed", "runtime_sec", "status", "pocket_center_x", "pocket_center_y", "pocket_center_z"
]

def write_sidecar(parquet_path: Path) -> None:
    sidecar = parquet_path.with_suffix(parquet_path.suffix + ".schema.json")
    meta = {"schema_version": SCHEMA_VERSION, "columns": RESULT_COLUMNS}
    sidecar.write_text(json.dumps(meta, indent=2), encoding="utf-8")

def validate_results(parquet_path: Path) -> Dict[str, str]:
    """
    Ensure required columns exist and optional sidecar carries matching schema_version.
    Returns dict with 'schema_version'.
    """
    import pyarrow.parquet as pq
    p = Path(parquet_path)
    table = pq.read_table(p)
    cols = set(table.column_names)
    required = {"id", "smiles", "vina_score"}
    missing = sorted(list(required - cols))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    # sidecar check
    sidecar = p.with_suffix(p.suffix + ".schema.json")
    if sidecar.exists():
        try:
            sc = json.loads(sidecar.read_text(encoding="utf-8"))
            if sc.get("schema_version") != SCHEMA_VERSION:
                raise ValueError("Schema version mismatch")
        except Exception as e:
            raise
    return {"schema_version": SCHEMA_VERSION}
