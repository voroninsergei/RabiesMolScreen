"""
Definitions of the standardized docking result schema for RabiesMol.

This module exposes constants describing the canonical column order and schema
version for result Parquet files produced by RabiesMol.  It also provides
helpers for writing a JSON sidecar alongside a Parquet file and validating
existing result files.

The sidecar stores the schema version and list of columns so that downstream
tools can verify compatibility without reading the entire Parquet dataset.
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict
import json

# Semantic version of the result schema.  Increment this when the
# expected columns or their semantics change.
SCHEMA_VERSION: str = "1.0"

# Core columns expected in any results Parquet.  The first three columns
# (id, smiles, vina_score) are mandatory and are validated in
# :func:`validate_results`.  Additional columns provide metadata about the
# receptor, docking engine and configuration used for the run.
RESULT_COLUMNS: List[str] = [
    "id",
    "smiles",
    "vina_score",
    "receptor",
    "engine",
    "exhaustiveness",
    "seed",
    "runtime_sec",
    "status",
    "pocket_center_x",
    "pocket_center_y",
    "pocket_center_z",
]

def write_sidecar(parquet_path: Path) -> None:
    """Write a sidecar JSON file containing schema metadata next to a Parquet file.

    Parameters
    ----------
    parquet_path:
        Path to the Parquet file.  The sidecar will be created at
        ``parquet_path`` with the suffix ``.schema.json``.
    """
    p = Path(parquet_path)
    sidecar = p.with_suffix(p.suffix + ".schema.json")
    meta = {"schema_version": SCHEMA_VERSION, "columns": RESULT_COLUMNS}
    sidecar.write_text(json.dumps(meta, indent=2), encoding="utf-8")

def validate_results(parquet_path: Path) -> Dict[str, str]:
    """Validate a results Parquet file against the expected schema.

    This function reads a Parquet file using :mod:`pyarrow` to check that
    required columns exist.  It also inspects the sidecar JSON (if present)
    to ensure the schema version matches the current :data:`SCHEMA_VERSION`.

    Parameters
    ----------
    parquet_path:
        Path to the Parquet file to validate.

    Returns
    -------
    Dict[str, str]
        A dictionary containing the `schema_version` when validation passes.

    Raises
    ------
    ValueError
        If required columns are missing or if the sidecar's schema_version does not match.
    """
    import pyarrow.parquet as pq  # Delayed import to avoid heavy dependency at module import

    p = Path(parquet_path)
    table = pq.read_table(p)
    cols = set(table.column_names)
    required = {"id", "smiles", "vina_score"}
    missing = sorted(list(required - cols))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    # Validate the schema version in the sidecar if it exists
    sidecar = p.with_suffix(p.suffix + ".schema.json")
    if sidecar.exists():
        try:
            sc = json.loads(sidecar.read_text(encoding="utf-8"))
            if sc.get("schema_version") != SCHEMA_VERSION:
                raise ValueError(
                    f"Schema version mismatch: expected {SCHEMA_VERSION}, found {sc.get('schema_version')}"
                )
        except Exception as e:
            # Re-raise any unexpected exceptions to aid debugging
            raise
    return {"schema_version": SCHEMA_VERSION}