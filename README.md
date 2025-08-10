# RabiesMol Screen

<p align="left">
  <a href="https://github.com/voroninsergei/RabiesMolScreen/actions/workflows/ci.yml">
    <img alt="CI" src="https://github.com/voroninsergei/RabiesMolScreen/actions/workflows/ci.yml/badge.svg">
  </a>
  <a href="https://codecov.io/gh/voroninsergei/RabiesMolScreen">
    <img alt="Coverage" src="https://codecov.io/gh/voroninsergei/RabiesMolScreen/branch/main/graph/badge.svg">
  </a>
  <a href="https://github.com/voroninsergei/RabiesMolScreen/wiki">
    <img alt="Wiki" src="https://img.shields.io/badge/wiki-GitHub%20Wiki-blue">
  </a>
  <a href="https://github.com/voroninsergei/RabiesMolScreen/blob/main/LICENSE">
    <img alt="License" src="https://img.shields.io/github/license/voroninsergei/RabiesMolScreen">
  </a>
</p>

Мини‑пайплайн для подготовки белков/лигангов, виртуального докинга,
перескоринга и генерации простого HTML‑отчёта.

## Установка
```bash
pip install -r requirements.txt        # базовые зависимости
# внешние бинарники: vina/smina/obabel — установите в системе (apt/brew/conda)
```

## Быстрый старт
```bash
python -m rabiesmol.cli prepare --proteins data/proteins --ligands examples --seed 42
python -m rabiesmol.cli dock --receptor data/prepared/proteins/receptor.pdbqt --ligands data/prepared/ligands --out outputs/run1
python -m rabiesmol.cli rescore --csv outputs/run1/results.csv --out-parquet outputs/run1/results.parquet
python -m rabiesmol.cli report --parquet outputs/run1/results.parquet --out reports/report.html
```

## Структура проекта
```
rabiesmol/        # пакет с кодом (CLI, prepare/docking/rescoring/…)
data/             # входные данные (PDB, SMILES) — по необходимости
examples/         # маленькие демо-лиганды
reports/          # HTML‑отчёты
```

## Воспроизводимость
- RDKit загружается лениво: CLI работает даже без RDKit; ADMET/PAINS помечаются как *skipped*.
- Внешние инструменты (vina/smina/obabel) не включены в `requirements.txt` — ставьте через пакетный менеджер.

## Документация
- Wiki: https://github.com/voroninsergei/RabiesMolScreen/wiki
- Docs (шаблон для GitHub Pages в `docs/`).

## Лицензия и вклад
Лицензия: MIT (см. `LICENSE`). Вклад приветствуется — см. [`CONTRIBUTING.md`](CONTRIBUTING.md).
