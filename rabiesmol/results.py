from pathlib import Path
import pandas as pd
from loguru import logger

RESULT_COLUMNS = [
    "id",
    "smiles",
    "vina_score",
    "smina_score",
    "consensus",
    "rfscore_vs",
    "gnina_score",
    "cns_mpo",
    "bbb_pred",
    "pains_flag",
    "alerts"
]

def to_parquet(df: pd.DataFrame, out_path: Path):
    # ensure all expected columns exist
    for col in RESULT_COLUMNS:
        if col not in df.columns:
            df[col] = None
    df = df[RESULT_COLUMNS]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, index=False)
    logger.info(f"Saved standardized results -> {out_path}")
