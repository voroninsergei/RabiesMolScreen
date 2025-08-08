from pydantic import BaseModel
from pathlib import Path

class PipelineConfig(BaseModel):
    proteins_dir: Path
    ligands_dir: Path
    docking_dir: Path
    prepared_proteins_dir: Path
    prepared_ligands_dir: Path
