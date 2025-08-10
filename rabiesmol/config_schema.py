
from __future__ import annotations
from pydantic import BaseModel, Field, field_validator
from pathlib import Path
from typing import Optional, List

class DockingConfig(BaseModel):
    grid_size: int = Field(20, ge=8, le=128)
    exhaustiveness: int = Field(8, ge=1, le=64)
    seed: Optional[int] = None
    threads: int = Field(1, ge=1, le=256)
    receptor: Path
    pocket_pdb: Optional[Path] = None

    @field_validator("receptor", mode="before")
    @classmethod
    def _check_receptor(cls, v):
        p = Path(v)
        if not p.exists():
            raise ValueError(f"Receptor file does not exist: {p}")
        return p

class Settings(BaseModel):
    output_dir: Path = Path("work")
    docking: DockingConfig

def validate_config(cfg: dict) -> Settings:
    return Settings.model_validate(cfg)
