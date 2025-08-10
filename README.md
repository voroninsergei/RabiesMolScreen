# RabiesMolScreen (arch-10/10 sample)

Мини‑пайплайн с правильной архитектурой: четкая граница с внешними бинарями,
плагиноподобные адаптеры и стандартизованный контракт результатов.

## CLI
```bash
python -m rabiesmol.cli.main prepare --proteins data/proteins --ligands examples/ligands --out data/prepared
python -m rabiesmol.cli.main dock --receptor data/prepared/proteins/receptor.pdbqt --ligands data/prepared/ligands --out outputs/run1 --engine vina
python -m rabiesmol.cli.main rescore-cmd --csv outputs/run1/results.csv --out-parquet outputs/run1/results.parquet
python -m rabiesmol.cli.main report --parquet outputs/run1/results.parquet --out reports/report.html
```

## Контракт результатов
См. `config/results.schema.json` и `rabiesmol/io/contract.py`.

## Заметки
- Для демонстрации `VinaEngine` генерирует фиктивные результаты; чтобы
  запускать реальную `vina/smina`, расширьте методы и используйте `runner.run(...)`.
