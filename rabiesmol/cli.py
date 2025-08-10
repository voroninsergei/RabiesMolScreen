
from __future__ import annotations
import typer
from typing import Optional
from pathlib import Path
from .logging_config import get_logger
from .prepare import prepare_proteins, prepare_ligands
from .docking import run_docking_batch
from .reporting import generate_report
from .validation import validate_and_cluster

app = typer.Typer(help="rabiesmol: Unified CLI for ligand prep, docking, and reporting.")
logger = get_logger(__name__)

@app.command()
def prepare(
    proteins: Path = typer.Option(None, help="Dir with input PDBs"),
    ligands: Path = typer.Option(None, help="Dir with input SDFs"),
    out_proteins: Path = typer.Option(Path("prepared/proteins"), help="Output dir for prepared receptors"),
    out_ligands: Path = typer.Option(Path("prepared/ligands"), help="Output dir for prepared ligands"),
    seed: int = typer.Option(0, "--seed", help="Random seed used for any randomized steps"),
    track: bool = typer.Option(False, "--track", help="Emit simple tracking output"),
):
    """Prepare receptors and ligands. In tests, external tools are mocked."""
    if proteins:
        prepare_proteins(proteins, out_proteins)
    if ligands:
        prepare_ligands(ligands, out_ligands)
    if track:
        logger.info(f"Tracking enabled; seed={seed}")
    return 0

@app.command()
def dock(
    protein: Path = typer.Argument(..., help="Input receptor PDBQT"),
    ligands_dir: Path = typer.Argument(..., help="Directory with ligands"),
    out_csv: Path = typer.Option(Path("out/scores.csv"), help="Where to save scores CSV"),
):
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df = run_docking_batch(protein, ligands_dir, out_csv.parent)
    df.to_csv(out_csv, index=False)
    logger.info(f"Docking completed -> {out_csv}"); return 0

@app.command()
def validate(
    input_csv: Path = typer.Argument(..., help="Docking results CSV"),
    out_csv: Path = typer.Argument(..., help="Where to write validated CSV"),
):
    validate_and_cluster(input_csv, out_csv)

@app.command()
def report(
    scores_csv: Path = typer.Argument(..., help="CSV produced by 'dock'"),
    out_html: Path = typer.Option(Path("reports/report.html"), help="Where to save the HTML report"),
):
    out_html.parent.mkdir(parents=True, exist_ok=True)
    generate_report(scores_csv, out_html)
    logger.info(f"Report saved -> {out_html}")

if __name__ == "__main__":
    app()


@app.command()
def rescore_cmd(
    docking_results_csv: Path = typer.Argument(..., help="CSV with docking results"),
    out_csv: Path = typer.Option(Path("out/rescored.csv"), help="Output CSV"),
    experimental: bool = typer.Option(False, "--experimental", help="Enable experimental ML rescoring (RF-Score-VS/GNINA stubs)"),
):
    from .rescoring import rescore
    rescore(docking_results_csv, out_csv, experimental=experimental)
    logger.info(f"Rescored -> {out_csv}")


@app.command()
def doctor():
    """Show environment diagnostics (binaries, CUDA, etc.)."""
    from .doctor import run_doctor
    import json
    info = run_doctor()
    typer.echo(json.dumps(info, indent=2))


@app.command("validate-config")
def validate_config(cfg_path: Path = typer.Argument(..., help="YAML config")):
    from .config_schema import validate_config as _validate
    from .config import import_config
    try:
        cfg = import_config(str(cfg_path))
        _validate(cfg)
    except Exception as e:
        typer.secho(f"Invalid config: {e}", fg=typer.colors.RED)
        raise typer.Exit(code=1)
    typer.secho("OK", fg=typer.colors.GREEN)


@app.command()
def snapshot(out_dir: Path = typer.Option(Path("snapshots"), help="Where to write snapshot"),
             note: str = typer.Option("", help="Optional note")):
    """Write a lightweight experiment manifest (inputs/outputs hashes)."""
    import json, platform, hashlib, os
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "platform": platform.platform(),
        "python": platform.python_version(),
        "note": note,
        "files": {}
    }
    for p in Path('.').glob('**/*'):
        if p.is_file() and p.stat().st_size < 5000000:
            h = hashlib.sha256(p.read_bytes()).hexdigest()[:16]
            manifest["files"][str(p)] = {"sha256": h, "size": p.stat().st_size}
    (out_dir/"manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    typer.echo(f"Snapshot -> {out_dir/ 'manifest.json'}")
