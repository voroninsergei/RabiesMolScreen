from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Dict

@dataclass(frozen=True)
class Ligand:
    id: str
    smiles: Optional[str] = None
    path: Optional[str] = None  # path to file if present

@dataclass(frozen=True)
class Protein:
    id: str
    path: str

# Simple contract for docking result rows
RESULT_SCHEMA_VERSION = "1.0.0"
RESULT_COLUMNS = {
    "schema_version": "string",
    "ligand_id": "string",
    "protein_id": "string",
    "engine": "string",
    "pose_id": "int64",
    "score": "float64",
    "unit": "string",
    "ligand_smiles": "string",
    "metadata": "string",  # JSON-encoded small blob (engine parameters, etc.)
}
