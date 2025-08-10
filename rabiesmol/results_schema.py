from __future__ import annotations
from pathlib import Path
import json
import pyarrow.parquet as pq
from pydantic import BaseModel, Field
from .results import RESULT_COLUMNS, SCHEMA_VERSION

class ResultSchema(BaseModel):
    schema_version: str = Field(default=SCHEMA_VERSION)
    required_columns: list[str] = Field(default_factory=lambda: ["id", "smiles", "vina_score"])
    allowed_columns: list[str] = Field(default_factory=lambda: list(RESULT_COLUMNS))

def validate_results(parquet_path: str | Path) -> None:
    p = Path(parquet_path)
    meta = pq.read_metadata(str(p))
    cols = set(meta.schema.names)
    model = ResultSchema()
    missing = set(model.required_columns) - cols
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    extras = cols - set(model.allowed_columns)
    if extras:
        raise ValueError(f"Unexpected columns present: {sorted(extras)}")
    # Sidecar schema check if exists
    sidecar = p.with_suffix(p.suffix + ".schema.json")
    if sidecar.exists():
        try:
            side = json.loads(sidecar.read_text(encoding="utf-8"))
            if str(side.get("schema_version")) != model.schema_version:
                raise ValueError(f"Schema version mismatch: file={side.get('schema_version')} expected={model.schema_version}")
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid schema sidecar JSON: {e}") from e
