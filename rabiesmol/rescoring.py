from __future__ import annotations
import json
from pathlib import Path
import pandas as pd
from loguru import logger
from rabiesmol.results import to_parquet

def rfscore_vs_stub(df: pd.DataFrame) -> pd.Series:
    logger.warning("RF-Score-VS is stubbed. Returning -vina_score as RF prediction.")
    return -df.get("vina_score", pd.Series([0.0] * len(df)))

def gnina_rescoring_stub(df: pd.DataFrame) -> pd.Series:
    logger.warning("GNINA rescoring is stubbed. Returning -consensus as GNINA prediction.")
    return -df.get("consensus", pd.Series([0.0] * len(df)))

def rescore(docking_results_csv: Path, out_csv: Path, experimental: bool = True):
    df = pd.read_csv(docking_results_csv)
    df["rfscore_vs"] = rfscore_vs_stub(df) if experimental else None
    df["gnina_score"] = gnina_rescoring_stub(df) if experimental else None
    df.to_csv(out_csv, index=False)
    logger.info(f"Rescoring complete -> {out_csv}")

def rescore_to_parquet(docking_results_csv: Path, out_parquet: Path, experimental: bool = True):
    df = pd.read_csv(docking_results_csv)
    df["rfscore_vs"] = rfscore_vs_stub(df) if experimental else None
    df["gnina_score"] = gnina_rescoring_stub(df) if experimental else None
    # ADMET & PAINS: подключаем опционально, чтобы не ломать лёгкий образ
    try:
        from rabiesmol.validation import compute_admet, apply_filters  # heavy; optional
        from rdkit import Chem  # type: ignore
        admet = compute_admet(df.get("smiles", []).tolist())
        df["cns_mpo"] = [a for a, b in admet]
        df["bbb_pred"] = [b for a, b in admet]
        mols = [Chem.MolFromSmiles(s or "") for s in df.get("smiles", [])]
        df["pains_flag"] = apply_filters(mols)
        try:
            from rabiesmol.validation import get_filter_alerts
            df["alerts"] = [json.dumps(a) for a in get_filter_alerts(mols)]
        except Exception:
            pass
    except Exception as e:
        logger.warning(f"ADMET/PAINS skipped: {e}")
        df["pains_flag"] = "skipped"
        try:
            out_parquet.with_suffix(".meta.json").write_text(json.dumps({"admet": "skipped"}, indent=2), encoding="utf-8")
        except Exception:
            pass
    to_parquet(df, out_parquet)
    logger.info(f"Saved parquet -> {out_parquet}")

# Backward-compatible CLI wrapper used in rabiesmol.cli
def rescore_cmd(in_csv: Path, out_csv: Path, experimental: bool = False) -> None:
    rescore(in_csv, out_csv, experimental=experimental)
    meta = {
        "rescoring": {
            "rfscore_vs": {"impl": "stub", "version": "0.1"},
            "gnina": {"impl": "stub", "version": "0.1"} if experimental else None,
        }
    }
    out_csv.with_suffix(".meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
