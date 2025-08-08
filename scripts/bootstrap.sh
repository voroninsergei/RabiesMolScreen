#!/usr/bin/env bash
set -euo pipefail

echo ">> Checking Git LFS..."
if ! command -v git &>/dev/null; then echo "git not found"; exit 1; fi
if ! git lfs env &>/dev/null; then
  echo "WARNING: Git LFS not configured. Run: git lfs install"
else
  echo "Git LFS OK"
fi

echo ">> Installing pre-commit hooks..."
if command -v pre-commit &>/dev/null; then
  pre-commit install
else
  echo "pre-commit not found; skipping (hooks will not run)"
fi

echo ">> Creating/Updating conda env..."
if command -v mamba &>/dev/null; then
  mamba env create -f environment.yml || mamba env update -f environment.yml
  echo "Activate with: mamba activate rabiesmol"
else
  echo "mamba not found; please install micromamba/mamba or use Docker"
fi

echo ">> Running mini E2E..."
bash scripts/mini_run.sh || (echo "Mini run failed"; exit 1)

echo "Bootstrap complete."
