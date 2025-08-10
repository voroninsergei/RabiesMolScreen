
from __future__ import annotations
try:
    import plotly.express as px
except Exception:
    px = None  # optional
from jinja2 import Environment, FileSystemLoader
import pandas as pd
from pathlib import Path

def generate_dashboard(data: pd.DataFrame, out_html: Path, template_dir: Path, snapshot_link: str = "#") -> None:
    if px is None:
        raise RuntimeError("Plotly is not installed. Install optional extra: rabiesmol[viz]")
    fig = px.scatter(data, x="ligand_id", y="score", color="score")
    plot_html = fig.to_html(full_html=False, include_plotlyjs="cdn")
    env = Environment(loader=FileSystemLoader(str(template_dir)))
    tpl = env.get_template("dashboard_template.html")
    html = tpl.render(plot_html=plot_html, top_hits=data.nsmallest(10, "score").to_dict(orient="records"), snapshot_link=snapshot_link)
    out_html.write_text(html, encoding="utf-8")
