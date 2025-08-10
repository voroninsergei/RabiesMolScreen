from __future__ import annotations
import json, pathlib, typer
from typing import Optional
from ..utils.logging import setup_logging
from ..adapters.vina import VinaEngine
from ..adapters.smina import SminaEngine
from ..services.prepare import prepare_proteins, prepare_ligands
from ..services.dock import run_docking
from ..services.rescore import rescore
from ..services.report import build_html_report

app = typer.Typer(help="RabiesMolScreen CLI - prepare -> dock -> rescore -> report")  # noqa

@app.callback()
def _main(log_level: str = typer.Option("INFO", help="Logging level")) -> None:
    setup_logging(log_level)

@app.command()
def prepare(proteins: str = typer.Option(..., exists=True, help="Path to proteins (PDB/PDBQT)"),
            ligands: str = typer.Option(..., exists=True, help="Path to ligands"),
            out: str = typer.Option("data/prepared", help="Output dir for prepared data")) -> None:
    prepare_proteins(proteins, pathlib.Path(out)/"proteins")
    prepare_ligands(ligands, pathlib.Path(out)/"ligands")

@app.command()
def dock(receptor: str = typer.Option(..., exists=True, help="Receptor .pdbqt"),
         ligands: str = typer.Option(..., exists=True, help="Folder with ligands .pdbqt"),
         out: str = typer.Option("outputs/run1", help="Output dir"),
         engine: str = typer.Option("vina", help="Docking engine (vina/smina)")) -> None:
    eng = VinaEngine() if engine == "vina" else SminaEngine()
    csv_path = run_docking(eng, receptor, ligands, out)
    typer.echo(csv_path)

@app.command()
def rescore_cmd(csv: str = typer.Option(..., exists=True),
                out_parquet: str = typer.Option("outputs/run1/results.parquet")) -> None:
    path = rescore(csv, out_parquet)
    typer.echo(path)

@app.command()
def report(parquet: str = typer.Option(..., exists=True),
           out: str = typer.Option("reports/report.html")) -> None:
    template_dir = str(pathlib.Path(__file__).resolve().parent.parent / "report")
    build_html_report(parquet, out, template_dir)
    typer.echo(out)

if __name__ == "__main__":
    app()
