ENV_NAME = rabiesmol

.PHONY: setup download prepare dock rescore report clean

setup:
	mamba env create -f environment.yml || mamba env update -f environment.yml
	echo "Activate with: mamba activate $(ENV_NAME)"

download:
	python scripts/fetch_structures.py --output data/proteins

prepare:
	python scripts/prepare_inputs.py --proteins data/proteins --ligands data/ligands

dock:
	python scripts/run_docking.py --proteins data/proteins --ligands data/ligands --out data/docking

rescore:
	python scripts/rescore.py --input data/docking --out data/rescored

report:
	jupyter nbconvert --execute reports/top_hits.ipynb --to html --output-dir reports/html

clean:
	rm -rf data/docking/* data/rescored/* reports/html/*

# Save standardized parquet after rescoring
rescore-parquet:
	python -c "from pathlib import Path; from rabiesmol.rescoring import rescore_to_parquet; rescore_to_parquet(Path('data/docking/results.csv'), Path('data/results/results.parquet'))"
