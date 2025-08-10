from __future__ import annotations
import pathlib
from typing import Any, Dict
import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape

def render_dashboard(df: pd.DataFrame, out_html: str, template_dir: str) -> None:
    env = Environment(
        loader=FileSystemLoader(template_dir),
        autoescape=select_autoescape(),
    )
    tpl = env.get_template("dashboard_template.html")
    html = tpl.render(
        n_rows=len(df),
        preview=df.head(20).to_html(index=False),
    )
    pathlib.Path(out_html).write_text(html, encoding="utf-8")
