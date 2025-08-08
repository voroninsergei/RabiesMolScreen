
from pathlib import Path
from rabiesmol.docking import run_vina, run_smina, consensus_docking

def test_docking_and_consensus(tmp_path: Path):
    prot = tmp_path / "prot_prepared.pdbqt"
    lig = tmp_path / "lig_prepared.pdbqt"
    prot.write_text("RECEPTOR PDBQT")
    lig.write_text("LIGAND PDBQT")

    out_dir = tmp_path / "dock"
    vina = run_vina(prot, lig, out_dir, exhaustiveness=8, seed=42)
    smina = run_smina(prot, lig, out_dir, exhaustiveness=8, seed=42)
    res = consensus_docking(prot, lig, out_dir)

    assert vina != 0.0
    assert smina != 0.0
    assert "consensus" in res
