
from __future__ import annotations
from pathlib import Path
import typer
from .prepare import prepare_ligands, prepare_proteins

app = typer.Typer(no_args_is_help=True)

@app.command()
def prepare(
    proteins: Path = typer.Option(Path("data/proteins"), help="Input proteins dir"),
    ligands: Path = typer.Option(Path("examples"), help="Input ligands dir"),
    out_proteins: Path = typer.Option(Path("data/prepared/proteins"), help="Output proteins"),
    out_ligands: Path = typer.Option(Path("data/prepared/ligands"), help="Output ligands"),
    ph: float = typer.Option(7.4, help="Ligand protonation pH"),
    seed: int = typer.Option(42, help="Seed for deterministic prep"),
    track: bool = typer.Option(False, help="No-op flag for CI tracking"),
):
    # Minimal: run mocked external tools
    prepare_proteins(proteins, out_proteins, ph=ph, keep_waters=[])
    # If a .smi file exists in ligands dir, use the first one
    smi_files = list(Path(ligands).glob("*.smi")) + list(Path(ligands).glob("*.smiles"))
    if smi_files:
        prepare_ligands(smi_files[0], out_ligands)

def main():
    app()

if __name__ == "__main__":
    main()
