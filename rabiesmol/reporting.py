
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
{% if embed %}
<style>
table { border-collapse: collapse; width: 100%; font-family: ui-sans-serif, system-ui; }
th, td { border: 1px solid #e5e7eb; padding: 6px; font-size: 14px; }
th { background: #f3f4f6; text-align: left; }
</style>
{% else %}
<link rel="stylesheet" href="https://cdn.datatables.net/2.0.3/css/dataTables.dataTables.min.css">
<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<script src="https://cdn.datatables.net/2.0.3/js/dataTables.min.js"></script>
{% endif %}
</head>
<body>
<h2>RabiesMol — docking hits (n={{n}})</h2>
{% if embed %}
<table>
<thead><tr>{% for k in header %}<th>{{k}}</th>{% endfor %}</tr></thead>
<tbody>
{% for r in rows %}<tr>{% for k in header %}<td>{{r.get(k)}}</td>{% endfor %}</tr>{% endfor %}
</tbody></table>
{% else %}
<table id='t' class='display' style='width:100%'></table>
<script>
const rows = {{ rows|tojson }};
const header = Object.keys(rows[0] || {});
$(document).ready(function() {
  $('#t').DataTable({
    data: rows,
    columns: header.map(h => ({title: h, data: h})),
    pageLength: 25
  });
});
</script>
{% endif %}
</body>
</html>
"""

def generate_report(in_csv: Path, out_html: Path, template: str | None = None, embed_assets: bool = False) -> None:
    """Generate interactive (CDN) or fully-embedded HTML report."""
    in_csv = Path(in_csv)
    if in_csv.suffix == ".parquet":
        df = pd.read_parquet(in_csv)
    else:
        df = pd.read_csv(in_csv)
        # Save Parquet twin for convenience
        try:
            from .results import to_parquet
            to_parquet(df.copy(), in_csv.with_suffix('.parquet'))
        except Exception:
            pass
    tpl = Template(template or DEFAULT_TEMPLATE)
    header = list(df.columns)
    html = tpl.render(n=len(df), rows=df.to_dict(orient="records"), header=header, embed=embed_assets)
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html, encoding="utf-8")
    logger.info(f"Report generated -> {out_html}")
