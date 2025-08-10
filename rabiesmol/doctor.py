
from __future__ import annotations
import shutil, subprocess, json, platform, os
from typing import Dict

def _bin_info(name: str) -> Dict[str, str | None]:
    path = shutil.which(name)
    ver = None
    if path:
        try:
            out = subprocess.check_output([name, "--version"], text=True, stderr=subprocess.STDOUT, timeout=5)
            ver = out.splitlines()[0].strip()
        except Exception:
            ver = "unknown"
    return {"path": path, "version": ver}

def run_doctor() -> dict:
    info = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "bins": {
            "smina": _bin_info("smina"),
            "fpocket": _bin_info("fpocket"),
            "gnina": _bin_info("gnina"),
            
            "vina": _bin_info("vina"),
            "obabel": _bin_info("obabel"),
            "nvidia-smi": _bin_info("nvidia-smi"),
            "cuda": _bin_info("cuda"),
        },
        "env": {
            "CUDA_VISIBLE_DEVICES": os.environ.get("CUDA_VISIBLE_DEVICES")
        }
    }
    return info
