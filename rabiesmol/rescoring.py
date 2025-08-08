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

def rescore(docking_results_csv: Path, out_csv: Path):
    df = pd.read_csv(docking_results_csv)
    df["rfscore_vs"] = rfscore_vs_stub(df)
    df["gnina_score"] = gnina_rescoring_stub(df)
    df.to_csv(out_csv, index=False)
    logger.info(f"Rescoring complete -> {out_csv}")

from rabiesmol.results import to_parquet


def rescore_to_parquet(docking_results_csv: Path, out_parquet: Path):
    df = pd.read_csv(docking_results_csv)
    df['rfscore_vs'] = rfscore_vs_stub(df)
    df['gnina_score'] = gnina_rescoring_stub(df)
    to_parquet(df, out_parquet)
