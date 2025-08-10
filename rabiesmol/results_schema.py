
from __future__ import annotations
import pyarrow.parquet as pq
from pydantic import BaseModel, Field
from typing import Literal

class ResultSchema(BaseModel):
    schema_version: Literal["1"] = "1"
    required_columns: list[str] = Field(default_factory=lambda: ["ligand_id","score","rmsd","pose_path"])

def validate_results(parquet_path: str) -> None:
    """Raise ValueError if parquet doesn't conform to expected schema."""
    meta = pq.read_metadata(parquet_path)
    cols = set(meta.schema.names)
    schema = ResultSchema()
    missing = set(schema.required_columns) - cols
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
