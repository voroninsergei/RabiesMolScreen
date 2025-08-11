import csv
import os
from pathlib import Path

import pytest

from rabiesmol.adapters.vina import VinaEngine
from rabiesmol.utils.runner import ExternalToolError


def test_vina_dock_creates_csv(tmp_path, monkeypatch):
    eng = VinaEngine()
    # Force runner.run to raise to trigger fallback
    monkeypatch.setattr("rabiesmol.utils.runner.run", lambda cmd, env=None, timeout=None: (_ for _ in ()).throw(ExternalToolError("test")))
    receptor = tmp_path / "rec.pdbqt"
    receptor.write_text("RECEPTOR")
    lig_dir = tmp_path / "ligs"
    lig_dir.mkdir()
    lig1 = lig_dir / "lig1.pdbqt"
    lig1.write_text("LIGAND")
    out_dir = tmp_path / "out"
    result_csv = eng.dock(str(receptor), [str(lig1)], str(out_dir))
    assert os.path.isfile(result_csv)
    # Ensure CSV contains a score column
    with open(result_csv) as f:
        rows = list(csv.DictReader(f))
        assert rows and "score" in rows[0]
