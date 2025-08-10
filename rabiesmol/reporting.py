from __future__ import annotations
from pathlib import Path
import pandas as pd
from jinja2 import Template
from .logging_config import get_logger

logger = get_logger(__name__)

DEFAULT_TEMPLATE = """
<!doctype html>
<html>
<head>
<meta charset='utf-8'><title>rabiesmol report</title>
<link rel="stylesheet" href="https://cdn.datatables.net/2.0.3/css/dataTables.dataTables.min.css">
<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<script src="https://cdn.datatables.net/2.0.3/js/dataTables.min.js"></script>
<style> body{font-family: system-ui, sans-serif; padding: 1rem;} table{width:100%;} </style>
</head>
<body>
<h1>rabiesmol: docking summary</h1>
<p>Total records: {{ n }} | Snapshot: <a href="../snapshot/manifest.json">manifest.json</a></p>
<table id="hits" class="display">
<thead><tr><th>Ligand</th><th>Vina score</th></tr></thead>
<tbody>
{% for r in rows %}
<tr><td>{{ r['ligand'] }}</td><td>{{ '%.3f'|format(r['vina_score']) }}</td></tr>
{% endfor %}
</tbody>
</table>
<script>
new DataTable('#hits', {paging: true, searching: true, order:[[1,'asc']]});
</script>
</body></html>
"""

def generate_report(in_csv: Path, out_html: Path, template: str | None = None) -> None:
    df = pd.read_csv(in_csv)
    tpl = Template(template or DEFAULT_TEMPLATE)
    html = tpl.render(n=len(df), rows=df.to_dict(orient="records"))
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html, encoding="utf-8")
    logger.info(f"Report generated -> {out_html}")
