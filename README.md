# RabiesMolScreen

Pipeline for molecular screening against RABV (Rabies Virus) targets.

## Overview

RabiesMolScreen is an open-source pipeline for virtual screening of small molecules against rabies virus proteins (e.g., L and P). It integrates molecular docking (AutoDock Vina), cheminformatics (RDKit), and optional ML rescoring.

## Quickstart (End-to-End Example)

```bash
# 1. Setup environment
mamba env create -f environment.yml
mamba activate rabiesmol

# 2. Download example proteins (stubbed)
make download

# 3. Prepare proteins and ligands
make prepare

# 4. Run docking
make dock

# 5. Rescore (if implemented)
make rescore

# 6. Generate report
make report
```

### Expected Output

- `data/proteins/` — input protein structures
- `data/ligands/` — input ligands
- `data/docking/` — docking results (.pdbqt)
- `data/results/` — rescoring results (if available)
- `reports/html/` — generated HTML reports

## Makefile Targets

| Target   | Description |
|----------|-------------|
| setup    | Install conda env |
| download | Fetch proteins (stub) |
| prepare  | Prepare proteins/ligands |
| dock     | Run docking |
| rescore  | Rescore results |
| report   | Generate HTML report |
| clean    | Clean intermediate files |

## Dependencies

See [`docs/concept.md`](docs/concept.md) for a roadmap and full dependency table.

## License

MIT


## Target Proteins: RABV L and P

The primary targets in this project are the **L** (RNA-dependent RNA polymerase) and **P** (phosphoprotein) proteins of the rabies virus (RABV).

### Data Sources

- **AlphaFold DB** — predicted protein structures for RABV proteins by UniProt IDs
- **Protein Data Bank (PDB)** — experimentally solved structures (if available)

Structures are stored with:
- **FASTA** sequences
- **PDB** coordinate files
- **metadata.json** containing download date, UniProt ID, and source

### Multiple Conformations

The pipeline supports:
- Different **domains** of the L protein (polymerase domains)
- L–P complexes for interface docking

See [`scripts/fetch_structures.py`](scripts/fetch_structures.py) for automated retrieval.


## Defaults & Reproducibility

Default decisions are tracked in `config/defaults.yaml` (pH, salts, grid box, exhaustiveness, seeds, structural waters). Edit this file to change project-wide defaults.

## Privacy & Licensing

All protein structures and compound examples are sourced from **public/open** databases (AlphaFold DB, UniProt; example ligands are simple SMILES provided here for demonstration). External tools/libraries (RDKit, OpenBabel, AutoDock Vina, fpocket, DeepChem, etc.) are used under their respective licenses; see their documentation for details.
