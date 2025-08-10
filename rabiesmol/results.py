
import json
from pathlib import Path
import pandas as pd
from loguru import logger

SCHEMA_VERSION = "v1"
RESULT_COLUMNS = [
    "id", "smiles", "vina_score", "smina_score", "consensus",
    "rfscore_vs", "gnina_score", "cns_mpo", "bbb_pred", "pains_flag", "alerts"
]

def to_parquet(df: pd.DataFrame, out_path: Path):
    # ensure all expected columns exist
    for col in RESULT_COLUMNS:
        if col not in df.columns:
            df[col] = None
    df = df[RESULT_COLUMNS]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, index=False)
    # write schema sidecar
    meta = {"schema_version": SCHEMA_VERSION, "columns": RESULT_COLUMNS}
    out_path.with_suffix(out_path.suffix + ".schema.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    logger.info(f"Saved standardized results -> {out_path}")
