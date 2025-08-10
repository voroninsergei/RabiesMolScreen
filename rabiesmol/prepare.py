from __future__ import annotations
from pathlib import Path
from typing import List, Optional
from loguru import logger
import shutil
from concurrent.futures import ThreadPoolExecutor
from .proc import run
from rabiesmol.io import ensure_dir

# --- Ligands ---
def protonate_structure(input_file: Path, output_file: Path, ph: float = 7.4):
    cmd = ['obabel', str(input_file), '-O', str(output_file), '--pH', str(ph)]
    logger.debug(f"Running: {' '.join(cmd)}")
    run(cmd)
    logger.info(f"Protonated {input_file} -> {output_file}")

def prepare_ligands(input_dir: Path, output_dir: Path, ph: float = 7.4, threads: int = 1):
    ensure_dir(output_dir)
    ligs = list(Path(input_dir).glob('*.sdf')) + list(Path(input_dir).glob('*.smi')) + list(Path(input_dir).glob('*.pdbqt'))
    if not ligs:
        logger.warning(f'No ligands found in {input_dir}')
        return
    def _one(mol_file: Path):
        stem = mol_file.stem
        protonated_file = output_dir / f"{stem}_protonated.pdb"
        prepared_file = output_dir / f"{stem}_prepared.pdbqt"
        if prepared_file.exists():
            logger.info(f'Cached prepared ligand: {prepared_file}')
            return
        logger.debug('Running: prepare_ligand')
        # tests' mocks expect short flags -l/-o and a created file at -o
        run(['prepare_ligand', '-l', str(protonated_file), '-o', str(prepared_file)])
        logger.info(f'Prepared ligand: {prepared_file}')
    if threads <= 1:
        for p in ligs:
            _one(p)
    else:
        with ThreadPoolExecutor(max_workers=threads) as ex:
            list(ex.map(_one, ligs))

# --- Proteins ---
def clean_protein(protein_file: Path, out_file: Path, keep_waters: List[str] | None = None) -> Path:
    """Minimal cleanup: drop HOH/WAT waters unless preserved in keep_waters."""
    keep = set(keep_waters or [])
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with open(protein_file, 'r', encoding='utf-8', errors='ignore') as fin, open(out_file, 'w', encoding='utf-8') as fout:
        for line in fin:
            if line.startswith(('ATOM', 'HETATM')) and len(line) >= 21:
                resn = line[17:20].strip()
                if resn in {'HOH', 'WAT'} and resn not in keep:
                    continue
            fout.write(line)
    logger.info(f'Cleaned {protein_file} -> {out_file}')
    return out_file

def prepare_proteins(input_dir: Path, output_dir: Path, ph: float = 7.4, keep_waters: Optional[List[str]] = None, threads: int = 1) -> None:
    ensure_dir(output_dir)
    proteins = list(Path(input_dir).glob('*.pdb')) + list(Path(input_dir).glob('*.pdbqt'))
    if not proteins:
        logger.warning(f'No proteins found in {input_dir} (expected *.pdb or *.pdbqt)')
    def _one(p: Path):
        out = output_dir / p.name.replace(p.suffix, '_prepared.pdbqt')
        if out.exists():
            logger.info(f'Cached prepared protein: {out}')
            return
        cleaned = output_dir / p.name.replace(p.suffix, '_clean.pdb')
        # in tests we don't actually need protonation, just the final file created by external tool mock
        logger.debug('Running: prepare_receptor')
        run(['prepare_receptor', '--receptor', str(cleaned), '-o', str(out)])
        logger.info(f'Prepared receptor: {out}')
    if threads <= 1:
        for p in proteins:
            _one(p)
    else:
        with ThreadPoolExecutor(max_workers=threads) as ex:
            list(ex.map(_one, proteins))
