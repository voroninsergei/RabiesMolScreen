from __future__ import annotations
from pathlib import Path
import pandas as pd
from jinja2 import Template
from .logging_config import get_logger

logger = get_logger(__name__)

TEMPLATE_FILE = Path('dashboard_template.html')
DEFAULT_TEMPLATE = """
<!doctype html>
<html><head><meta charset='utf-8'><title>rabiesmol report</title></head>
<body>
<h1>rabiesmol: docking summary</h1>
<p>Total poses: {{ n }}</p>
{{ table }}
</body></html>
"""

def generate_report(scores_csv: Path, out_html: Path) -> Path:
    df = pd.read_csv(scores_csv)
    html_table = df.sort_values('vina_score').to_html(index=False)
    html = Template(TEMPLATE_FILE.read_text(encoding='utf-8') if TEMPLATE_FILE.exists() else DEFAULT_TEMPLATE).render(n=len(df), table=html_table)
    out_html.write_text(html, encoding='utf-8')
    return out_html
