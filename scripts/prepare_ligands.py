#!/usr/bin/env python3
"""
prepare_ligands.py

This script converts a file containing SMILES strings into 3D coordinates and
writes out PDB files. Conversion to PDBQT can be done with external tools
such as Open Babel or MGLTools. Each line of the input file should contain
at least a SMILES string; an optional second column can be used to name the
molecule.

Usage:
    python prepare_ligands.py smiles.txt -o output_dir

Requirements: RDKit installed in your Python environment.
"""
import os
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem

def generate_3d(mol: Chem.Mol) -> Chem.Mol:
    """Generates 3D coordinates and optimizes geometry for a molecule."""
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    return mol

def process_smiles_file(smiles_file: str, output_dir: str) -> None:
    """Reads SMILES from a file and writes PDB files to the output directory."""
    with open(smiles_file, "r", encoding="utf-8") as f:
        for idx, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            smiles = parts[0]
            name = parts[1] if len(parts) > 1 else f"ligand_{idx+1}"
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"Skipping invalid SMILES: {smiles}")
                continue
            mol3d = generate_3d(mol)
            pdb_path = os.path.join(output_dir, f"{name}.pdb")
            Chem.MolToPDBFile(mol3d, pdb_path)
            print(f"Saved {pdb_path}")

def main():
    parser = argparse.ArgumentParser(description="Convert SMILES to PDB using RDKit")
    parser.add_argument("smiles", help="Path to input SMILES file (one per line)")
    parser.add_argument("-o", "--out", default="ligands", help="Output directory")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    process_smiles_file(args.smiles, args.out)

if __name__ == "__main__":
    main()
