from __future__ import annotations
import json, hashlib, platform
from typing import Optional
from pathlib import Path
import typer
from .logging_config import get_logger

app = typer.Typer(help="rabiesmol: Unified CLI for ligand prep, docking, and reporting.")
logger = get_logger(__name__)

def _run_validation(*args, **kwargs):
    """Lazy-import validator to avoid pulling RDKit unless needed."""
    try:
        from .validation import validate_and_cluster  # type: ignore
    except Exception as e:
        typer.secho("Validation requires optional dependencies (RDKit). Install `rabiesmol[validation]`.\n" + str(e), fg=typer.colors.YELLOW)
        raise typer.Exit(code=2)
    return validate_and_cluster(*args, **kwargs)

@app.command()
def prepare(
    proteins: Path = typer.Option(Path("data/proteins"), help="Input proteins dir"),
    ligands: Path = typer.Option(Path("examples"), help="Input ligands dir"),
    out_proteins: Path = typer.Option(Path("data/prepared/proteins"), help="Output proteins dir"),
    out_ligands: Path = typer.Option(Path("data/prepared/ligands"), help="Output ligands dir"),
    ph: float = typer.Option(7.4, help="Ligand protonation pH"),
    threads: int = typer.Option(1, help="Parallel workers for prep"),
    seed: Optional[int] = typer.Option(None, help="Random seed (optional)") ,
    track: bool = typer.Option(False, help="No-op flag for CI tracking"),
):
    """Prepare inputs: proteins and ligands."""
    from .prepare import prepare_proteins, prepare_ligands
    logger.info(f"Using seed={seed}")
    prepare_proteins(proteins, out_proteins, ph=ph, keep_waters=[], threads=threads)
    prepare_ligands(ligands, out_ligands, ph=ph, threads=threads)

@app.command()
def dock(
    protein: Path,
    ligands_dir: Path,
    out_csv: Path = typer.Option(Path("outputs/scores.csv"), help="Output CSV with scores"),
    threads: int = typer.Option(1, help="Parallel workers"),
    cache_dir: Optional[Path] = typer.Option(Path("outputs/cache"), help="Cache directory for docking results"),
):
    from .docking import run_docking_batch
    df = run_docking_batch(protein, ligands_dir, out_dir=Path("outputs/poses"), threads=threads, cache_dir=cache_dir)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    typer.echo(f"Saved scores -> {out_csv}")

@app.command("rescore-cmd")
def rescore_cmd(
    in_csv: Path,
    out_csv: Path = typer.Option(Path("outputs/rescored.csv"), help="Output CSV"),
    experimental: bool = typer.Option(False, help="Enable experimental steps"),
):
    from .rescoring import rescore_cmd as _rescore
    _rescore(in_csv, out_csv, experimental=experimental)
    typer.echo(f"Rescored -> {out_csv}")

@app.command()
def report(in_csv: Path, out_html: Path = Path("reports/report.html")):
    from .reporting import generate_report
    generate_report(in_csv, out_html)

@app.command("validate-config")
def validate_config(config_yaml: Path = Path("config/defaults.yaml")):
    if not config_yaml.exists():
        typer.echo(f"Config not found: {config_yaml}")
        raise typer.Exit(code=2)
    import yaml
    from .config_schema import validate_config as _validate
    cfg = yaml.safe_load(config_yaml.read_text(encoding="utf-8")) or {}
    _validate(cfg)
    typer.echo(f"Config validated: {config_yaml}")

@app.command()
def doctor():
    from .doctor import run_doctor
    import json as _json
    info = run_doctor()
    typer.echo(_json.dumps(info, indent=2))

@app.command()
def snapshot(note: Optional[str] = typer.Option("local", help="Note for snapshot"), out_dir: Path = Path("snapshot")):
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "platform": platform.platform(),
        "python": platform.python_version(),
        "note": note,
        "files": {}
    }
    for p in Path('.').glob('**/*'):
        if p.is_file() and p.stat().st_size < 5_000_000:
            h = hashlib.sha256(p.read_bytes()).hexdigest()[:16]
            manifest["files"][str(p)] = {"sha256": h, "size": p.stat().st_size}
    (out_dir/"manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    typer.echo(f"Snapshot -> {out_dir/ 'manifest.json'}")

if __name__ == "__main__":
    app()
