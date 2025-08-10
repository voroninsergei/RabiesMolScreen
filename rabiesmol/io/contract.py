from __future__ import annotations
import json
from typing import Iterable, Dict, Any, List
import pandas as pd

from ..domain.models import RESULT_COLUMNS, RESULT_SCHEMA_VERSION

def schema() -> Dict[str, str]:
    return dict(RESULT_COLUMNS)

def validate_df(df: pd.DataFrame) -> List[str]:
    problems: List[str] = []
    for col, dtype in RESULT_COLUMNS.items():
        if col not in df.columns:
            problems.append(f"missing column: {col}")
            continue
        if dtype == "string" and not pd.api.types.is_string_dtype(df[col]):
            problems.append(f"column {col} must be string")
        if dtype == "float64" and not pd.api.types.is_float_dtype(df[col]):
            problems.append(f"column {col} must be float")
        if dtype == "int64" and not pd.api.types.is_integer_dtype(df[col]):
            problems.append(f"column {col} must be int")
    extra = set(df.columns) - set(RESULT_COLUMNS)
    if extra:
        problems.append("unexpected columns: " + ", ".join(sorted(extra)))
    return problems

def attach_schema_meta(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["schema_version"] = RESULT_SCHEMA_VERSION
    return df

def write_parquet(df: pd.DataFrame, out_path: str) -> None:
    df = attach_schema_meta(df)
    problems = validate_df(df)
    if problems:
        raise ValueError("Result schema validation errors: \n- " + "\n- ".join(problems))
    df.to_parquet(out_path, index=False)  # requires pyarrow
