
from __future__ import annotations
from pathlib import Path
import re
from typing import Optional
from .proc import run

def _parse_vina_log(log_file: Path) -> float:
    """
    Return the first best energy from Vina/Smina log.
    Robust to headers, tabs/spaces, locales.
    """
    txt = Path(log_file).read_text(encoding="utf-8", errors="ignore").splitlines()
    nums = []
    for line in txt:
        line = line.strip()
        if not line or not re.match(r"^\d+", line):
            continue
        parts = re.split(r"[\s,]+", line)
        if len(parts) >= 2:
            try:
                nums.append(float(parts[1]))
            except ValueError:
                continue
    if not nums:
        raise ValueError(f"No energies found in {log_file}")
    return float(nums[0])

def _fallback_score(receptor: Path, ligand: Path) -> float:
    # Simple deterministic fallback (worst score)
    return 0.0
