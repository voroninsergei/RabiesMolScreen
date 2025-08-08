
from pathlib import Path
import json
from rabiesmol.pockets import detect_pockets_fpocket

def test_detect_pockets_fpocket(tmp_path: Path):
    protein = tmp_path / "prot.pdb"
    protein.write_text("ATOM")

    out_json = tmp_path / "pocket.json"
    detect_pockets_fpocket(protein, out_json)

    data = json.loads(out_json.read_text())
    assert "grid_center" in data and "grid_size" in data
