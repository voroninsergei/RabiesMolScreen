from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, DataStructs
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from .logging_config import get_logger

logger = get_logger(__name__)

def load_filters() -> FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.LILLY)
    return FilterCatalog(params)

def apply_filters(mols):
    cat = load_filters()
    flags = []
    for m in mols:
        try:
            flags.append(cat.HasMatch(m))
        except Exception:
            flags.append(True)
    return flags

def cluster_ecfp4(mols, threshold: float = 0.6):
    fps = [rdMolDescriptors.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in mols]
    clusters = []
    used = set()
    for i, fi in enumerate(fps):
        if i in used:
            continue
        cluster = [i]
        used.add(i)
        for j in range(i+1, len(fps)):
            if j in used:
                continue
            sim = DataStructs.TanimotoSimilarity(fi, fps[j])
            if sim >= threshold:
                cluster.append(j)
                used.add(j)
        clusters.append(cluster)
    return clusters

def validate_and_cluster(input_csv: Path, out_csv: Path):
    df = pd.read_csv(input_csv)
    mols = [Chem.MolFromSmiles(smi) for smi in df['smiles']]
    df['filter_flag'] = apply_filters(mols)
    df['alerts'] = ['' if ok == False else '' for ok in df['filter_flag']]
    clusters = cluster_ecfp4(mols)
    cluster_id = {}
    for cid, members in enumerate(clusters):
        for m in members:
            cluster_id[m] = cid
    df['cluster_id'] = [cluster_id[i] for i in range(len(mols))]
    df.to_csv(out_csv, index=False)
    logger.info(f"Validation complete -> {out_csv}")


def validate_hits(input_csv: Path, out_csv: Path):
    return validate_and_cluster(input_csv, out_csv)


def compute_admet(smiles_list):
    """Compute lightweight CNS MPO and BBB predictions for demonstration."""
    results = []
    from rdkit.Chem import Descriptors
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi) if smi else None
        if not mol:
            results.append((None, False)); continue
        logp = Descriptors.MolLogP(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        mw = Descriptors.MolWt(mol)
        cns_mpo = max(0.0, 6.0 - abs(2.0 - logp) - (tpsa/60.0))
        bbb = (tpsa < 90) and (mw < 500)
        results.append((cns_mpo, bbb))
    return results
