# RabiesMol — скрининг лигандов (мини-пайплайн)
[![License](https://img.shields.io/github/license/voroninsergei/RabiesMolScreen)](https://github.com/voroninsergei/RabiesMolScreen/blob/main/LICENSE)
[![Docs](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://voroninsergei.github.io/RabiesMolScreen/)
[![Wiki](https://img.shields.io/badge/wiki-on%20GitHub-lightgrey)](https://github.com/voroninsergei/RabiesMolScreen/wiki)

[![CI](https://img.shields.io/github/actions/workflow/status/voroninsergei/RabiesMolScreen/ci.yml?branch=main&label=CI)](https://github.com/voroninsergei/RabiesMolScreen/actions)
[![codecov](https://img.shields.io/codecov/c/github/voroninsergei/RabiesMolScreen?label=coverage)](https://codecov.io/gh/voroninsergei/RabiesMolScreen)

> Небольшой учебный пайплайн для подготовки входов, докинга и простого рескоринга.
> Работает без RDKit (lazy‑import), а при наличии RDKit добавляет ADMET/PAINS.

---

## Методология (кратко)
1. **Подготовка белков/лигандов**: лёгкая чистка PDB, протонирование и конвертация в PDBQT.
2. **Докинг**: запуск внешней утилиты (например, *smina*/*vina*) и парсинг логов (устойчивый к различным форматам).
3. **Рескоринг**: объединение базовых метрик; при наличии RDKit — фильтры PAINS и расчёт простых ADMET‑признаков.
4. **Отчёт**: HTML‑дашборд и выгрузка результатов в Parquet (+ сайдкар `*.schema.json`).

```
flow: prepare -> dock -> rescore -> report
```

### Мини‑flowchart (ASCII)
```
inputs/               data/prepared/           outputs/
  proteins/  ----->     proteins/  ----------------->  docking/
  ligands/   ----->     ligands/   --\               \-> results.parquet
                                 (parallel)            \-> report.html
```

## Требования
- **Python**: 3.9+
- **Обязательно**: `pandas`, `pyarrow`, `loguru`
- **Опционально**: `rdkit` (для валидации/ADMET; не обязателен), внешние бинарники `smina`/`vina`, `obabel`

Установка зависимостей:
```bash
pip install -r requirements.txt
# dev-инструменты (pytest/black/ruff и т.п.)
pip install -r dev-requirements.txt
```

Conda‑окружение (опционально):
```bash
conda env create -f environment.yml
conda activate rabiesmol
```

## Примеры запуска
```bash
# 1) Подготовка входов
python -m rabiesmol.cli prepare --proteins data/proteins --ligands examples --out-proteins data/prepared/proteins --out-ligands data/prepared/ligands --seed 42 --threads 2

# 2) Докинг
python -m rabiesmol.cli dock data/prepared/ligands data/prepared/proteins --out outputs --cache-dir outputs/cache

# 3) Рескоринг
python -m rabiesmol.cli rescore outputs/docking/results.csv --out outputs/results.parquet

# 4) Отчёт
python -m rabiesmol.cli report outputs/results.parquet --out reports/report.html
```

## Структура проекта
```
rabiesmol/        # пакет с кодом (CLI, prepare/docking/rescoring/...)
data/             # входные данные (PDB, SMILES) — по необходимости
examples/         # маленькие демо-лиганды
reports/          # HTML-отчёты
```

## Замечания по воспроизводимости
- Внешние инструменты (vina/smina/obabel) не входят в `requirements.txt` — установите в системе/образе.
- Пайплайн устойчив к отсутствию RDKit: CLI работает, а ADMET помечается `"skipped"`.

## Лицензия и вклад
Лицензия: MIT (см. `LICENSE`). Вклад приветствуется — см. [`CONTRIBUTING.md`](CONTRIBUTING.md).
