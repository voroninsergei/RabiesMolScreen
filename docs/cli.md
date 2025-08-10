# CLI

## prepare
Подготовка входов (белки/лигандов).
```
python -m rabiesmol.cli prepare --proteins DIR --ligands DIR --out-proteins DIR --out-ligands DIR [--threads N] [--seed INT]
```

## dock
Запуск докинга.
```
python -m rabiesmol.cli dock PREP_LIGANDS PREP_PROTEINS --out OUTDIR [--cache-dir DIR]
```

## rescore
Рескоринг CSV результатов:
```
python -m rabiesmol.cli rescore results.csv --out results.parquet
```

## report
Генерация отчёта:
```
python -m rabiesmol.cli report results.parquet --out reports/report.html
```
