from __future__ import annotations
import os, pathlib, json
import pandas as pd

def rescore(csv_path: str, out_parquet: str) -> str:
    df = pd.read_csv(csv_path)
    # Cheap "rescoring": identity mapping, ensure presence of column 'score'
    df["rescored"] = df["score"]
    df.to_parquet(out_parquet, index=False)
    return out_parquet
