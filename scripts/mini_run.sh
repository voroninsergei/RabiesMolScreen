#!/usr/bin/env bash
set -euo pipefail
python scripts/fetch_structures.py
mkdir -p data/ligands
cp examples/ligands_mini.smi data/ligands/ligands.smi
make prepare
make dock
make rescore
make report
echo 'Mini run finished.'
