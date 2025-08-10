.PHONY: prepare dock rescore report report-html validate-config doctor snapshot
prepare:
	python -m rabiesmol.cli prepare --proteins data/proteins --ligands examples --track
dock:
	python -m rabiesmol.cli dock data/protein/receptor.pdbqt data/prepared/ligands --out-csv outputs/scores.csv
rescore:
	python -m rabiesmol.cli rescore-cmd outputs/scores.csv --out-csv outputs/rescored.csv --experimental
report:
	python -m rabiesmol.cli report outputs/scores.csv --out-html reports/report.html
report-html: report
validate-config:
	python -m rabiesmol.cli validate-config config/defaults.yaml
doctor:
	python -m rabiesmol.cli doctor
snapshot:
	python -m rabiesmol.cli snapshot --note "local"
