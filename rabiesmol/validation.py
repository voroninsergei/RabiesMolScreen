from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, DataStructs
from rdkit.Chem import FilterCatalog, FilterCatalogParams
from loguru import logger
import random

def load_filters() -> FilterCatalog.FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.LILLY)
    return FilterCatalog.FilterCatalog(params)

def apply_filters(mols):
    catalog = load_filters()
    flags = []
    for mol in mols:
        entry = catalog.GetFirstMatch(mol)
        flags.append(bool(entry))
    return flags

def cluster_ecfp4(mols, cutoff: float = 0.7):
    fps = [rdMolDescriptors.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in mols]
    clusters = []
    assigned = set()
    for i, fp in enumerate(fps):
        if i in assigned:
            continue
        cluster = [i]
        assigned.add(i)
        for j in range(i + 1, len(fps)):
            if j in assigned:
                continue
            sim = DataStructs.TanimotoSimilarity(fp, fps[j])
            if sim >= cutoff:
                cluster.append(j)
                assigned.add(j)
        clusters.append(cluster)
    return clusters

def validate_hits(input_csv: Path, out_csv: Path):
    df = pd.read_csv(input_csv)
    mols = [Chem.MolFromSmiles(smi) for smi in df["smiles"]]

    df["filter_flag"] = apply_filters(mols)

    clusters = cluster_ecfp4(mols)
    cluster_id = {}
    for cid, members in enumerate(clusters):
        for m in members:
            cluster_id[m] = cid
    df["cluster_id"] = [cluster_id[i] for i in range(len(mols))]

    df.to_csv(out_csv, index=False)
    logger.info(f"Validation complete -> {out_csv}")
