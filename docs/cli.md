# CLI

## prepare

Подготовка входных данных. Команда принимает путь к каталогу с белками (`--proteins`), путь к каталогу с лигандами (`--ligands`) и путь для сохранения подготовленных данных (`--out`). В каталоге `out` будут созданы поддиректории `proteins` и `ligands`.

```bash
python -m rabiesmol.cli.main prepare --proteins data/proteins --ligands data/ligands --out data/prepared
```

## dock

Запуск докинга. Требуется подготовленный рецептор (файл `.pdbqt`) и каталог с подготовленными лигандами. Результаты сохраняются в CSV.

```bash
python -m rabiesmol.cli.main dock --receptor data/prepared/proteins/receptor.pdbqt --ligands data/prepared/ligands --out outputs/run1 [--engine vina|smina]
```

## rescore

Рескоринг CSV результатов. Команда принимает путь к CSV‑файлу с результатами (`--csv`) и путь к выходному файлу Parquet (`--out-parquet`).

```bash
python -m rabiesmol.cli.main rescore-cmd --csv outputs/run1/results.csv --out-parquet outputs/run1/results.parquet
```

## report

Генерация HTML‑отчёта из Parquet. Необходимо указать путь к Parquet‑файлу и место сохранения HTML‑отчёта.

```bash
python -m rabiesmol.cli.main report --parquet outputs/run1/results.parquet --out reports/report.html
```
