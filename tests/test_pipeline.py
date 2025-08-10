import subprocess
import sys
from pathlib import Path
import pytest

@pytest.fixture(autouse=True)
def patch_subprocess_run(monkeypatch):
    def dummy_run(cmd, check=True, **kwargs):
        # Simulate external binary by creating output files if specified
        if "--out" in cmd:
            out_idx = cmd.index("--out") + 1
            Path(cmd[out_idx]).write_text("DUMMY OUTPUT")
        if "-o" in cmd:
            out_idx = cmd.index("-o") + 1
            Path(cmd[out_idx]).write_text("DUMMY OUTPUT")
        if "--log" in cmd:
            log_idx = cmd.index("--log") + 1
            Path(cmd[log_idx]).write_text("1       -7.5    0.0  0.0")  # vina/smina log format
        return subprocess.CompletedProcess(cmd, 0)
    monkeypatch.setattr(subprocess, "run", dummy_run)
    yield

def test_vina_and_smina(tmp_path):
    from rabiesmol.docking import run_vina, run_smina
    protein = tmp_path / "prot.pdbqt"
    ligand = tmp_path / "lig.pdbqt"
    protein.write_text("")
    ligand.write_text("")
    score_vina = run_vina(protein, ligand, tmp_path)
    score_smina = run_smina(protein, ligand, tmp_path)
    assert score_vina == -7.5
    assert score_smina == -7.5

def test_prepare(tmp_path):
    from rabiesmol.prepare import clean_protein
    pdb_in = tmp_path / "prot.pdb"
    pdb_out = tmp_path / "clean.pdb"
    pdb_in.write_text("HETATM      1  O   HOH A   1      11.000  11.000  11.000\n")
    clean_protein(pdb_in, pdb_out)
    assert pdb_out.exists()

def test_results_parquet(tmp_path):
    import pandas as pd
    from rabiesmol.results import to_parquet
    df = pd.DataFrame([{"id": "mol1", "smiles": "CCO", "vina_score": -7.5}])
    out_file = tmp_path / "results.parquet"
    to_parquet(df, out_file)
    assert out_file.exists()
