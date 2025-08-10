from __future__ import annotations
import logging
from typing import Optional

def setup_logging(level: str = "INFO") -> None:
    numeric = getattr(logging, level.upper(), logging.INFO)
    logging.basicConfig(
        level=numeric,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )

def get_logger(name: Optional[str] = None) -> logging.Logger:
    return logging.getLogger(name or "rabiesmol")
