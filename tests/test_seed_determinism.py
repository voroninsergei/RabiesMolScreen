
from pathlib import Path
from rabiesmol.docking_ext import batch_dock, DockParams

def test_seed_determinism(tmp_path: Path):
    ligands = ["CCO", "CCN", "CCC"]
    recs = [tmp_path/"r1.pdb", tmp_path/"r2.pdb"]
    # create dummy receptor files
    for r in recs: r.write_text("RECEPTOR", encoding="utf-8")
    p = DockParams(receptor_paths=recs, seed=42, restarts=2)
    res1, _ = batch_dock(ligands, p, out_dir=tmp_path/"out1", threads=2, cache_dir=tmp_path/"cache")
    res2, _ = batch_dock(ligands, p, out_dir=tmp_path/"out2", threads=2, cache_dir=tmp_path/"cache")
    assert res1 == res2
