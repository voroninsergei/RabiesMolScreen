# RabiesMolScreen

RabiesMolScreen is an open-source initiative to develop therapeutic options for symptomatic rabies by identifying blood–brain barrier (BBB) penetrating small-molecule inhibitors targeting the rabies virus (RABV) L (polymerase) and P (phosphoprotein) proteins.

## Goal

To perform in-silico screening of drug-like molecules for their affinity to RABV targets and BBB permeability, followed by in-vitro validation.

## Background

Rabies is almost always fatal once clinical symptoms appear, and no effective treatment exists today. One of the major obstacles is the blood–brain barrier, which limits delivery of therapeutic molecules to infected neurons. This project focuses on discovering small molecules that can cross the BBB and disrupt viral replication.

## Repository Structure

- `scripts/` – Python tools for fetching protein structures, preparing ligands and running docking.
- `docs/` – background documents, concept note and team information.
- `data/` – (placeholder) input data such as PDB files and compound libraries.
- `LICENSE` – license file (MIT License).

## Getting Started

1. Clone this repository.
2. Install dependencies (e.g. RDKit, AutoDock Vina, DeepChem).
3. Use `scripts/fetch_structures.py` to download RABV protein structures.
4. Prepare ligand libraries and perform docking.

Contributions and feedback are welcome. Please see `docs/concept.md` for an overview of the project roadmap.
