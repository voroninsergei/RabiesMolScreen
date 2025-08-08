# Initial Infrastructure PR

This PR adds:
- Conda environment and Dockerfile (Vina, RDKit, OpenBabel)
- Makefile targets (setup, download, prepare, dock, rescore, report)
- `rabiesmol/` package (I/O, preparation, docking, rescoring, validation, pockets)
- GitHub Actions: CI (ruff/black/pytest/smoke) and Docker build to GHCR
- Documentation: README, docs/concept.md, style guide, templates
- Examples: `examples/ligands_mini.smi` (5 ligands) and `examples/ligands_20.smi` (20 ligands)
- Standardized results format (`results.parquet`) and reporting notebook

## How to run a quick demo (3–5 ligands)

```bash
mamba env create -f environment.yml
mamba activate rabiesmol

# fetch RABV L/P
python scripts/fetch_structures.py

# copy mini ligands into data/ligands
mkdir -p data/ligands
cp examples/ligands_mini.smi data/ligands/ligands.smi

# prepare, dock, rescore, report
make prepare
make dock
make rescore
make report
```

Checklist:
- [ ] CI green (lint + tests + smoke)
- [ ] Docker image published to GHCR
- [ ] Results generated and report produced


## 📷 Screenshots

**Dashboard Overview**  
![Dashboard Placeholder](docs/images/dashboard_placeholder.png)

**Trend Graph**  
![Trend Placeholder](docs/images/trend_placeholder.png)
