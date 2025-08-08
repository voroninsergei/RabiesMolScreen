from pathlib import Path
import subprocess
import json
from loguru import logger

def detect_pockets_fpocket(protein_file: Path, output_json: Path):
    subprocess.run(["fpocket", "-f", str(protein_file)], check=True)
    pockets_dir = protein_file.with_suffix("").name + "_out"
    pocket_info_file = Path(pockets_dir) / "pockets" / "pocket0" / "pocket.pdb"
    # Здесь можно парсить координаты для grid-box
    pocket_data = {
        "protein": str(protein_file),
        "pocket_file": str(pocket_info_file),
        "grid_center": [0.0, 0.0, 0.0],  # TODO: вычислить из pocket.pdb
        "grid_size": [20.0, 20.0, 20.0]
    }
    output_json.write_text(json.dumps(pocket_data, indent=2))
    logger.info(f"Pocket detection output saved to {output_json}")

def detect_pockets_dogsite(protein_file: Path, output_json: Path):
    logger.warning("DoGSite integration not implemented; requires web API or local tool.")
