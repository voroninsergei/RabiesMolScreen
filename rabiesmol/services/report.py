from __future__ import annotations
import pathlib, pandas as pd
from ..report.html import render_dashboard

def build_html_report(parquet_path: str, out_html: str, template_dir: str) -> None:
    df = pd.read_parquet(parquet_path)
    render_dashboard(df, out_html, template_dir=template_dir)
