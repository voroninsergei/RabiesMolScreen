from __future__ import annotations
from typing import List, Sequence
import subprocess
from .logging_config import get_logger

logger = get_logger(__name__)

class ExternalRunError(RuntimeError):
    def __init__(self, cmd: Sequence[str], returncode: int, stdout: str, stderr: str):
        super().__init__(f"Command failed ({returncode}): {' '.join(cmd)}\nSTDOUT:\n{stdout}\nSTDERR:\n{stderr}")
        self.cmd = list(cmd)
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr

class ExternalTimeout(ExternalRunError):
    pass

def run(cmd: List[str], timeout: int = 600) -> subprocess.CompletedProcess:
    cmd_str = ' '.join(cmd)
    logger.debug("Running: %s", cmd_str)
    try:
        res = subprocess.run(cmd, check=False, text=True, capture_output=True, timeout=timeout)
    except subprocess.TimeoutExpired as e:
        raise ExternalTimeout(cmd, -9, e.stdout or "", e.stderr or "") from e
    if res.returncode != 0:
        raise ExternalRunError(cmd, res.returncode, res.stdout or "", res.stderr or "")
    return res
