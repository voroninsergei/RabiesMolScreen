
from pathlib import Path
import json
from jinja2 import Template
import plotly.express as px
import pandas as pd

def build_report(results_json: Path, template_path: Path, out_html: Path) -> None:
    results = json.loads(results_json.read_text())
    df = pd.DataFrame(results).sort_values("score")
    fig = px.scatter(df, x="ligand_id", y="score")
    plot_html = fig.to_html(full_html=False, include_plotlyjs="cdn")
    tmpl = Template(template_path.read_text())
    html = tmpl.render(table=df.head(20).to_html(), plot=plot_html)
    out_html.write_text(html, encoding="utf-8")
