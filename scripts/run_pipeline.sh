#!/usr/bin/env bash
set -euo pipefail
python -m rabiesmol.cli.main prepare --proteins data/proteins --ligands examples/ligands --out data/prepared
python -m rabiesmol.cli.main dock --receptor data/prepared/proteins/receptor.pdbqt --ligands data/prepared/ligands --out outputs/run1 --engine vina
python -m rabiesmol.cli.main rescore --csv outputs/run1/results.csv --out-parquet outputs/run1/results.parquet
python -m rabiesmol.cli.main report --parquet outputs/run1/results.parquet --out reports/report.html
