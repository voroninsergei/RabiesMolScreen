
from __future__ import annotations
from pathlib import Path
from typing import List, Dict
import json
import pandas as pd
from jinja2 import Template
import plotly.express as px

def generate_interactive_report(results: List[Dict], template_path: Path, out_html: Path) -> Path:
    df = pd.DataFrame(results)
    df = df.sort_values("score").reset_index(drop=True)
    fig = px.scatter(df, x=range(len(df)), y="score", hover_data=["ligand_id","receptor"])
    plot_html = fig.to_html(full_html=False, include_plotlyjs="cdn")
    tmpl = Template(Path(template_path).read_text(encoding="utf-8"))
    html = tmpl.render(table=df.head(20).to_html(index=False), plot=plot_html)
    out_html.write_text(html, encoding="utf-8")
    return out_html
