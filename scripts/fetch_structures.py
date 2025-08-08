#!/usr/bin/env python
import requests
from pathlib import Path
import json
from datetime import date
from loguru import logger

ALPHAFOLD_URL = "https://alphafold.ebi.ac.uk/files/{uniprot}-F1-model_v4.pdb"
FASTA_URL = "https://rest.uniprot.org/uniprotkb/{uniprot}.fasta"

# Example UniProt IDs for RABV L and P proteins
TARGETS = {
    "L": "P0C569",  # RABV L protein
    "P": "P0C568",  # RABV P protein
}

OUT_DIR = Path("data/proteins")

def download_file(url: str, dest: Path):
    r = requests.get(url)
    r.raise_for_status()
    dest.write_bytes(r.content)
    logger.info(f"Downloaded {url} -> {dest}")

def fetch_protein(uniprot_id: str, name: str):
    protein_dir = OUT_DIR / name
    protein_dir.mkdir(parents=True, exist_ok=True)

    pdb_file = protein_dir / f"{name}.pdb"
    fasta_file = protein_dir / f"{name}.fasta"
    meta_file = protein_dir / "metadata.json"

    # Download PDB from AlphaFold
    download_file(ALPHAFOLD_URL.format(uniprot=uniprot_id), pdb_file)

    # Download FASTA from UniProt
    download_file(FASTA_URL.format(uniprot=uniprot_id), fasta_file)

    # Metadata
    metadata = {
        "uniprot_id": uniprot_id,
        "name": name,
        "source": {
            "alphafold": ALPHAFOLD_URL.format(uniprot=uniprot_id),
            "fasta": FASTA_URL.format(uniprot=uniprot_id)
        },
        "download_date": str(date.today())
    }
    meta_file.write_text(json.dumps(metadata, indent=2))
    logger.info(f"Saved metadata for {name}")

if __name__ == "__main__":
    for name, uid in TARGETS.items():
        fetch_protein(uid, name)
