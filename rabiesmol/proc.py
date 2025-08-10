
from __future__ import annotations
from typing import Sequence, Optional
import subprocess, shlex
from loguru import logger

class ExternalTimeout(Exception):
    pass

def run(cmd: Sequence[str], timeout: Optional[int] = None, check: bool = True):
    logger.debug(f"Run: {' '.join(shlex.quote(c) for c in cmd)}")
    try:
        res = subprocess.run(list(cmd), capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired as e:
        logger.error(f"Timeout: {cmd}")
        raise ExternalTimeout(str(e)) from e
    if check and res.returncode != 0:
        logger.error(f"Command failed rc={res.returncode}: {res.stderr}")
        raise subprocess.CalledProcessError(res.returncode, cmd, res.stdout, res.stderr)
    return res
