from __future__ import annotations
import csv, json, os, pathlib
from typing import Iterable, List, Dict
import pandas as pd

def read_vina_csv(path: str) -> pd.DataFrame:
    return pd.read_csv(path)
