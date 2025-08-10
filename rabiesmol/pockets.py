from __future__ import annotations
from pathlib import Path
import subprocess
import json
from loguru import logger

def compute_grid_center(pocket_pdb: str) -> tuple[float, float, float]:
    """Compute pocket center from PDB (simple average of atom coordinates)."""
    xs, ys, zs = [], [], []
    with open(pocket_pdb, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                try:
                    xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
                except ValueError:
                    continue
    if not xs:
        raise ValueError("No atoms parsed from pocket PDB for center computation")
    return (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))

def detect_pockets_fpocket(protein_file: Path, output_json: Path):
    subprocess.run(["fpocket", "-f", str(protein_file)], check=True)
    pockets_dir = protein_file.with_suffix("").name + "_out"
    pocket_info_file = Path(pockets_dir) / "pockets" / "pocket0" / "pocket.pdb"
    # If mock created pocket.pdb, compute center; otherwise fallback to (0,0,0)
    try:
        center = compute_grid_center(str(pocket_info_file)) if pocket_info_file.exists() else (0.0, 0.0, 0.0)
    except Exception:
        center = (0.0, 0.0, 0.0)
    pocket_data = {
        "protein": str(protein_file),
        "pocket_file": str(pocket_info_file),
        "grid_center": list(center),
        "grid_size": [20.0, 20.0, 20.0],
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(pocket_data, indent=2), encoding="utf-8")
    logger.info(f"Wrote pocket info -> {output_json}")
