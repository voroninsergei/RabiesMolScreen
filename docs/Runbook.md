
# Prod Runbook

This page documents a strict-mode run of the pipeline.

## Requirements
- Python 3.9+
- External tools: `vina` or `smina`, `obabel`

## Commands
```bash
# prepare demo inputs
python -m rabiesmol.cli prepare --threads 2

# docking with strict mode (no fallbacks)
python -m rabiesmol.cli dock data/protein/receptor.pdbqt data/prepared/ligands --threads 2 --pool thread

# rescore to parquet
python -m rabiesmol.cli rescore outputs/hits.csv --out-parquet outputs/hits.parquet

# generate air-gapped report
python -m rabiesmol.cli report outputs/hits.parquet --out-html reports/report.html --embed-assets
```
