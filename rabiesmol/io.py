from pathlib import Path
from typing import Union
from loguru import logger
import shutil

def ensure_dir(path: Union[str, Path]) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Ensured directory exists: {p}")
    return p

def copy_file(src: Union[str, Path], dst: Union[str, Path]) -> None:
    shutil.copy2(src, dst)
    logger.debug(f"Copied {src} -> {dst}")
