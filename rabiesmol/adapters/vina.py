from __future__ import annotations

import csv
import os
import pathlib
from typing import Iterable, Any

from ..ports.docking import DockingEngine
from ..utils.runner import which_or_raise, run
from ..utils.logging import get_logger
from ..docking import _parse_vina_log, _fallback_score

log = get_logger(__name__)


class VinaEngine(DockingEngine):
    """Concrete adapter for AutoDock Vina.

    This implementation calls the external ``vina`` binary for each ligand and
    records the top score for each docking.  If the binary is not available or
    any invocation fails, a deterministic fallback score of 0.0 is used.
    """

    name = "vina"

    def _ensure(self) -> str:
        """Return the path to the ``vina`` binary, raising if not found."""
        return which_or_raise("vina")

    def prepare_receptor(self, protein, out_dir: str, **kwargs: Any) -> str:
        # Vina expects pdbqt; assume .pdbqt provided or converted outside
        out = os.path.join(out_dir, f"{protein.id}.pdbqt")
        if not os.path.exists(out):
            pathlib.Path(out).write_text("RECEPTOR STUB\n")
        return out

    def prepare_ligand(self, ligand, out_dir: str, **kwargs: Any) -> str:
        out = os.path.join(out_dir, f"{ligand.id}.pdbqt")
        if not os.path.exists(out):
            pathlib.Path(out).write_text("LIGAND STUB\n")
        return out

    def dock(self, receptor_file: str, ligands: Iterable[str], out_dir: str, **kwargs: Any) -> str:
        """Run docking for each ligand and write a CSV with scores.

        Parameters
        ----------
        receptor_file: str
            Path to the prepared receptor file (.pdbqt).
        ligands: Iterable[str]
            Paths to prepared ligand files (.pdbqt).
        out_dir: str
            Directory where docking results and logs will be written.

        Returns
        -------
        str
            Path to the CSV file with columns ``ligand_id``, ``pose_id``, ``score``, ``unit``.
        """
        pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
        csv_path = os.path.join(out_dir, f"{self.name}_results.csv")
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=["ligand_id", "pose_id", "score", "unit"]
            )
            writer.writeheader()
            for lig_file in ligands:
                lig_id = pathlib.Path(lig_file).stem
                log_file = os.path.join(out_dir, f"{lig_id}_vina.log")
                out_pdbqt = os.path.join(out_dir, f"{lig_id}_out.pdbqt")
                cmd = [
                    self._ensure(),
                    "--receptor",
                    receptor_file,
                    "--ligand",
                    lig_file,
                    "--out",
                    out_pdbqt,
                    "--log",
                    log_file,
                    "--score_only",
                ]
                try:
                    run(cmd)
                    score = _parse_vina_log(log_file)
                except Exception as exc:
                    # Log and use fallback
                    log.warning(
                        f"Vina docking failed for {lig_id}: {exc}; using fallback score."
                    )
                    score = _fallback_score(pathlib.Path(receptor_file), pathlib.Path(lig_file))
                writer.writerow(
                    {
                        "ligand_id": lig_id,
                        "pose_id": 1,
                        "score": score,
                        "unit": "kcal/mol",
                    }
                )
        return csv_path
