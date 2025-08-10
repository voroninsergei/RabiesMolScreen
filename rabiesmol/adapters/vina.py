from __future__ import annotations
import csv, json, os, pathlib
from typing import Iterable, Any
from ..ports.docking import DockingEngine
from ..utils.runner import which_or_raise, run
from ..utils.logging import get_logger

log = get_logger(__name__)

class VinaEngine(DockingEngine):
    name = "vina"

    def _ensure(self) -> str:
        return which_or_raise("vina")

    def prepare_receptor(self, protein, out_dir: str, **kwargs: Any) -> str:
        # Vina expects pdbqt; assume .pdbqt provided or converted outside
        # Here we are no-op; in real life, call obabel to convert
        out = os.path.join(out_dir, f"{protein.id}.pdbqt")
        if not os.path.exists(out):
            # Create a tiny stub to avoid failing in demonstration runs
            pathlib.Path(out).write_text("RECEPTOR STUB\n")
        return out

    def prepare_ligand(self, ligand, out_dir: str, **kwargs: Any) -> str:
        out = os.path.join(out_dir, f"{ligand.id}.pdbqt")
        if not os.path.exists(out):
            pathlib.Path(out).write_text("LIGAND STUB\n")
        return out

    def dock(self, receptor_file: str, ligands: Iterable[str], out_dir: str, **kwargs: Any) -> str:
        # For illustrative purposes, we do NOT call vina. We simulate results.
        # In production, you'd run: run([self._ensure(), '--receptor', receptor_file, ...])
        csv_path = os.path.join(out_dir, "vina_results.csv")
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["ligand_id", "pose_id", "score", "unit"]) 
            writer.writeheader()
            for i, lig_file in enumerate(ligands, start=1):
                lig_id = pathlib.Path(lig_file).stem
                writer.writerow({"ligand_id": lig_id, "pose_id": 1, "score": -7.5 - i*0.1, "unit": "kcal/mol"})
        return csv_path
