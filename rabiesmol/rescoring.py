import json
from pathlib import Path
import pandas as pd
from loguru import logger

def rfscore_vs_stub(docking_results: pd.DataFrame) -> pd.Series:
    # Stub for RF-Score-VS (real model would require features)
    logger.warning("RF-Score-VS is stubbed. Returning -score as RF prediction.")
    return -docking_results["vina_score"]

def gnina_rescoring_stub(docking_results: pd.DataFrame) -> pd.Series:
    # Stub for GNINA rescoring
    logger.warning("GNINA rescoring is stubbed. Returning -consensus as GNINA prediction.")
    return -docking_results["consensus"]

def rescore(docking_results_csv: Path, out_csv: Path, experimental: bool = True):
    df = pd.read_csv(docking_results_csv)
    df["rfscore_vs"] = rfscore_vs_stub(df) if experimental else None
    df["gnina_score"] = gnina_rescoring_stub(df) if experimental else None
    df.to_csv(out_csv, index=False)
    logger.info(f"Rescoring complete -> {out_csv}")

from rabiesmol.results import to_parquet
from rabiesmol.validation import compute_admet, apply_filters
from rdkit import Chem


def rescore_to_parquet(docking_results_csv: Path, out_parquet: Path, experimental: bool = True):
    df = pd.read_csv(docking_results_csv)
    df['rfscore_vs'] = rfscore_vs_stub(df) if experimental else None
    df['gnina_score'] = gnina_rescoring_stub(df) if experimental else None
    to_parquet(df, out_parquet)


def _write_metadata(out_csv: Path, experimental: bool):
    meta = {
        "rescoring": {
            "experimental": experimental,
            "rfscore_vs": {"impl": "stub", "version": "0.1"} if experimental else None,
            "gnina": {"impl": "stub", "version": "0.1"} if experimental else None,
        }
    }
    out_csv.with_suffix(".meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

# ADMET columns
admet = compute_admet(df.get("smiles", []).tolist())
df["cns_mpo"] = [a for a, b in admet]
df["bbb_pred"] = [b for a, b in admet]
try:
    mols = [Chem.MolFromSmiles(s) for s in df.get("smiles", [])]
    df["pains_flag"] = apply_filters(mols)
except Exception:
    df["pains_flag"] = None
