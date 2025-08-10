# Getting Started

## Установка
- `pip install -r requirements.txt`
- (опц.) `conda env create -f environment.yml && conda activate rabiesmol`

## Мини-пайплайн
```bash
python -m rabiesmol.cli prepare --proteins data/proteins --ligands examples --out-proteins data/prepared/proteins --out-ligands data/prepared/ligands
python -m rabiesmol.cli dock data/prepared/ligands data/prepared/proteins --out outputs
python -m rabiesmol.cli rescore outputs/docking/results.csv --out outputs/results.parquet
python -m rabiesmol.cli report outputs/results.parquet --out reports/report.html
```
