from __future__ import annotations
import os
from . import registry
from ..utils.runner import which_or_raise, run
from ..utils.logging import get_logger
from ..domain.models import Ligand, Protein
from typing import Optional

log = get_logger(__name__)

@registry.register_tool("obabel")
def ensure_obabel() -> str:
    """Return path to obabel binary or raise a clear error."""
    return which_or_raise("obabel")

def smiles_to_mol(smiles: str, out_file: str, obabel_path: Optional[str] = None) -> str:
    bin_path = obabel_path or ensure_obabel()
    run([bin_path, "-:'" + smiles + "'", "-O", out_file])
    if not os.path.exists(out_file):
        raise RuntimeError("obabel produced no output: " + out_file)
    return out_file
