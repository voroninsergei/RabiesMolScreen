import requests
import pathlib
import argparse

# List of PDB IDs to fetch (example: Rabies virus L-polymerase and a related structure)
PDB_IDS = ["6UEB", "6LGX"]

def fetch(pdb_id: str, out_dir: pathlib.Path) -> None:
    """Download a PDB file from RCSB and save it to the output directory."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    out_path = out_dir / f"{pdb_id}.pdb"
    out_path.write_bytes(response.content)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch PDB structures for RabiesMolScreen")
    parser.add_argument("-o", "--out", default="structures", help="Output directory for downloaded PDB files")
    args = parser.parse_args()

    output_dir = pathlib.Path(args.out)
    output_dir.mkdir(parents=True, exist_ok=True)

    for pid in PDB_IDS:
        fetch(pid, output_dir)
        print(f"\u2713 {pid}.pdb saved")
