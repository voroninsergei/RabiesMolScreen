
#!/usr/bin/env python
from __future__ import annotations
from pathlib import Path
from datetime import date
from typing import Optional
import json
import hashlib

try:
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry
except Exception:
    requests = None
import urllib.request as _urllib
from loguru import logger

ALPHAFOLD_URL = "https://alphafold.ebi.ac.uk/files/{uniprot}-F1-model_v4.pdb"
FASTA_URL = "https://rest.uniprot.org/uniprotkb/{uniprot}.fasta"

# Example UniProt IDs for RABV L and P proteins
TARGETS = {
    "L": "P0C569",  # RABV L protein
    "P": "P0C568",  # RABV P protein
}

def _session(total_retries: int = 3, backoff: float = 0.5) -> requests.Session:
    s = requests.Session()
    retry = Retry(
        total=total_retries,
        read=total_retries,
        connect=total_retries,
        backoff_factor=backoff,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET",),
        raise_on_status=False,
    )
    s.mount("https://", HTTPAdapter(max_retries=retry))
    s.mount("http://", HTTPAdapter(max_retries=retry))
    return s

def _download(url: str, dest: Path, timeout: float = 10.0, cache_dir: Optional[Path] = None) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    # Simple cache by URL hash
    if cache_dir:
        cache_dir.mkdir(parents=True, exist_ok=True)
        key = hashlib.sha256(url.encode("utf-8")).hexdigest()[:16]
        cache_path = cache_dir / key
        if cache_path.exists():
            dest.write_bytes(cache_path.read_bytes())
            logger.info(f"Loaded from cache -> {dest}")
            return dest
    with _session() as s:
        with s.get(url, timeout=timeout, stream=True) as r:
            r.raise_for_status()
            tmp = dest.with_suffix(dest.suffix + ".part")
            with open(tmp, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
            tmp.replace(dest)
    if cache_dir:
        (cache_dir / key).write_bytes(dest.read_bytes())
    logger.info(f"Downloaded {url} -> {dest}")
    return dest

def fetch_protein(uniprot_id: str, name: str, out_dir: Path = Path("data/structures")) -> dict:
    out_dir = Path(out_dir)
    pdb_file = out_dir / f"{name}.pdb"
    fasta_file = out_dir / f"{name}.fasta"
    meta_file = out_dir / f"{name}.json"
    cache = out_dir / ".cache"

    _download(ALPHAFOLD_URL.format(uniprot=uniprot_id), pdb_file, cache_dir=cache)
    _download(FASTA_URL.format(uniprot=uniprot_id), fasta_file, cache_dir=cache)

    metadata = {
        "uniprot_id": uniprot_id,
        "name": name,
        "source": {
            "alphafold": ALPHAFOLD_URL.format(uniprot=uniprot_id),
            "fasta": FASTA_URL.format(uniprot=uniprot_id)
        },
        "download_date": str(date.today())
    }
    meta_file.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    logger.info(f"Saved metadata for {name}")
    return metadata

if __name__ == "__main__":
    for name, uid in TARGETS.items():
        fetch_protein(uid, name)
