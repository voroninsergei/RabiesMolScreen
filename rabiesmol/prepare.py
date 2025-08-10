
from __future__ import annotations
from pathlib import Path
from typing import Iterable, List, Optional
from .proc import run
from loguru import logger

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def prepare_ligands(input_smiles: Path, output_dir: Path) -> None:
    ensure_dir(output_dir)
    # For each line in SMILES file create a fake prepared ligand via external tool (mocked in tests)
    lines = [l.strip() for l in input_smiles.read_text(encoding="utf-8").splitlines() if l.strip()]
    for i, smi in enumerate(lines, 1):
        out = output_dir / f"lig_{i}.pdbqt"
        run(["prepare_ligand", "-l", smi, "-o", str(out)])

def prepare_proteins(input_dir: Path, output_dir: Path, ph: float = 7.4, keep_waters: Optional[List[str]] = None) -> None:
    ensure_dir(output_dir)
    for prot in list(Path(input_dir).glob("*.pdb")) + list(Path(input_dir).glob("*.pdbqt")):
        run(["fpocket", "-f", str(prot)])
        # The mock in tests will create pockets/pocket0/pocket.pdb etc.
        logger.debug(f"Prepared protein: {prot.name}")
