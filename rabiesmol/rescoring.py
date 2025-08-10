
from __future__ import annotations
import json, time
from pathlib import Path
import pandas as pd
from loguru import logger
from .results import to_parquet

def rfscore_vs_stub(df: pd.DataFrame) -> pd.Series:
    logger.warning("RF-Score-VS is stubbed. Returning -vina_score as RF prediction.")
    return -df.get("vina_score", pd.Series([0.0] * len(df)))

def gnina_rescoring_stub(df: pd.DataFrame) -> pd.Series:
    logger.warning("GNINA rescoring is stubbed. Returning -consensus as GNINA prediction.")
    return -df.get("consensus", pd.Series([0.0] * len(df)))

def rescore(docking_results_csv: Path, out_csv: Path, experimental: bool = True) -> None:
    df = pd.read_csv(docking_results_csv)
    df["rfscore_vs"] = rfscore_vs_stub(df) if experimental else None
    df["gnina_score"] = gnina_rescoring_stub(df) if experimental else None
    df.to_csv(out_csv, index=False)
    logger.info(f"Rescoring complete -> {out_csv}")

def rescore_to_parquet(docking_results_csv: Path, out_parquet: Path, experimental: bool = True) -> None:
    df = pd.read_csv(docking_results_csv)
    df["rfscore_vs"] = rfscore_vs_stub(df) if experimental else None
    df["gnina_score"] = gnina_rescoring_stub(df) if experimental else None
    # Optional ADMET / PAINS with lazy RDKit imports
    try:
        from rabiesmol.validation import compute_admet, apply_filters, get_filter_alerts
        from rdkit import Chem  # type: ignore
        mols = [Chem.MolFromSmiles(s or "") for s in df.get("smiles", [])]
        admet = compute_admet([s or "" for s in df.get("smiles", [])])
        df["cns_mpo"] = [a for a, _ in admet]
        df["bbb_pred"] = [b for _, b in admet]
        df["pains_flag"] = apply_filters(mols)
        df["alerts"] = [json.dumps(a) for a in get_filter_alerts(mols)]
    except Exception as e:
        logger.warning(f"ADMET/PAINS skipped: {e}")
        # be explicit – no silent True/False
        df["pains_flag"] = None
        df["alerts"] = [json.dumps(["rdkit_error"])] * len(df)
        try:
            out_parquet.with_suffix(".meta.json").write_text(json.dumps({"admet": "skipped", "reason": str(e)}, indent=2), encoding="utf-8")
        except Exception:
            pass
    to_parquet(df, out_parquet)
    logger.info(f"Saved parquet -> {out_parquet}")
