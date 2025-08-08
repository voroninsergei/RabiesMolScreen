
from pathlib import Path
from rabiesmol.prepare import prepare_proteins, prepare_ligands

def test_prepare_proteins_and_ligands(tmp_path: Path):
    proteins_in = tmp_path / "in_proteins"
    ligands_in = tmp_path / "in_ligands"
    proteins_in.mkdir()
    ligands_in.mkdir()

    # Create dummy PDB and SDF
    (proteins_in / "L.pdb").write_text("ATOM")
    (ligands_in / "lig1.sdf").write_text("$$$$")

    out_prot = tmp_path / "prepared_proteins"
    out_lig = tmp_path / "prepared_ligands"

    prepare_proteins(proteins_in, out_prot, ph=7.4, keep_waters=[])
    prepare_ligands(ligands_in, out_lig, ph=7.4)

    assert list(out_prot.glob("*_prepared.pdbqt"))
    assert list(out_lig.glob("*_prepared.pdbqt"))
