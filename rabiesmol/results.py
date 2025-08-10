"""
Functions for standardizing RabiesMol docking results.

This module exposes a helper to persist docking results to a Parquet file
alongside a JSON sidecar containing the schema version and column order.  The
schema definition itself is housed in :mod:`rabiesmol.results_schema` and
should not be duplicated here.  Import the canonical `RESULT_COLUMNS` and
`SCHEMA_VERSION` from that module to avoid inconsistencies.

Users should call :func:`to_parquet` with a pandas DataFrame containing at
least the required result columns.  Any missing optional columns will be
added with null values so that downstream tooling can rely on a consistent
schema.  A corresponding ``.schema.json`` file will be written via
``write_sidecar`` to store the schema metadata.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable
import pandas as pd
# Use loguru for structured logging if available, otherwise fall back to the
# standard library logger.  This ensures the function works even when
# `loguru` isn't installed (e.g., minimal testing environments).
try:
    from loguru import logger  # type: ignore
except ImportError:  # pragma: no cover - fallback used only if loguru missing
    import logging
    logger = logging.getLogger(__name__)

# Pull the schema definition and sidecar writer from a single source of truth.
from .results_schema import RESULT_COLUMNS, SCHEMA_VERSION, write_sidecar

__all__ = ["to_parquet", "RESULT_COLUMNS", "SCHEMA_VERSION"]

def to_parquet(df: pd.DataFrame, out_path: Path) -> None:
    """Persist a DataFrame of docking results to Parquet with schema metadata.

    Parameters
    ----------
    df:
        The input pandas DataFrame containing docking results.  Columns
        corresponding to the required schema (defined in
        :data:`rabiesmol.results_schema.RESULT_COLUMNS`) will be included in
        the output.  Missing columns will be filled with ``None`` values.

    out_path:
        Destination path for the Parquet file.  A JSON sidecar will be
        created at ``out_path`` with the suffix ``.schema.json`` via
        :func:`rabiesmol.results_schema.write_sidecar`.

    Notes
    -----
    This function deliberately avoids re-defining the schema version or
    column order to prevent divergence from :mod:`rabiesmol.results_schema`.
    """
    # Ensure all expected columns exist; add missing optional columns as None.
    for col in RESULT_COLUMNS:
        if col not in df.columns:
            df[col] = None
    # Reorder columns to match the canonical schema
    df_out = df[RESULT_COLUMNS]
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # Write the Parquet file without an index
    df_out.to_parquet(out_path, index=False)
    # Persist the accompanying schema sidecar
    write_sidecar(out_path)
    logger.info(f"Saved standardized results -> {out_path}")