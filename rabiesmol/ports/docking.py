from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Iterable, Dict, Any, Optional
from ..domain.models import Ligand, Protein

class DockingEngine(ABC):
    """Abstract boundary for any docking engine (vina, smina, autodock etc.)."""

    name: str = "engine"

    @abstractmethod
    def prepare_receptor(self, protein: Protein, out_dir: str, **kwargs: Any) -> str:
        """Make engine-specific receptor file (e.g. .pdbqt). Returns path."""

    @abstractmethod
    def prepare_ligand(self, ligand: Ligand, out_dir: str, **kwargs: Any) -> str:
        """Make engine-specific ligand file (e.g. .pdbqt). Returns path."""

    @abstractmethod
    def dock(self, receptor_file: str, ligands: Iterable[str], out_dir: str, **kwargs: Any) -> str:
        """Run docking. Returns path to CSV with columns: ligand_id, pose_id, score, unit."""

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} name={self.name!r}>"
