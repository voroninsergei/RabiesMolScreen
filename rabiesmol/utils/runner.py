from __future__ import annotations
import shutil
import subprocess
from dataclasses import dataclass
from typing import List, Dict, Optional
from .logging import get_logger

log = get_logger(__name__)

class ExternalToolError(RuntimeError):
    pass

@dataclass
class CmdResult:
    returncode: int
    stdout: str
    stderr: str

def which_or_raise(tool: str) -> str:
    path = shutil.which(tool)
    if not path:
        raise ExternalToolError(f"Required tool '{tool}' not found in PATH. "
                                "Install it via apt/brew/conda and try again.")
    return path

def run(cmd: List[str], env: Optional[Dict[str, str]] = None, timeout: Optional[int] = None) -> CmdResult:
    log.debug("Running external command: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=timeout)
    if proc.returncode != 0:
        raise ExternalToolError(f"Command failed ({proc.returncode}): {' '.join(cmd)}\n{proc.stderr}")
    return CmdResult(proc.returncode, proc.stdout, proc.stderr)
