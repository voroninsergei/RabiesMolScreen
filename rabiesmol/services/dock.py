from __future__ import annotations
import os, pathlib
from typing import Iterable, List, Dict
import pandas as pd
from ..ports.docking import DockingEngine
from ..io.readers import read_vina_csv
from ..io.contract import attach_schema_meta
from ..domain.models import RESULT_SCHEMA_VERSION

def run_docking(engine: DockingEngine, receptor: str, ligands_dir: str, out_dir: str) -> str:
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
    lig_files: List[str] = [str(p) for p in pathlib.Path(ligands_dir).glob("*.pdbqt")]
    csv_path = engine.dock(receptor, lig_files, out_dir)
    # Normalize to contract
    df = read_vina_csv(csv_path)
    df = df.assign(engine=engine.name, protein_id=pathlib.Path(receptor).stem, ligand_smiles="", metadata="{}")  # fill minimal
    df = attach_schema_meta(df)
    out_csv = os.path.join(out_dir, "results.csv")
    df.to_csv(out_csv, index=False)
    return out_csv
