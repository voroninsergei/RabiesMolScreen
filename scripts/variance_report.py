#!/usr/bin/env python
import argparse
from pathlib import Path
import pandas as pd
from rabiesmol.docking import run_vina, run_smina
from loguru import logger

def main(out: Path):
    proteins = list(Path('data/prepared_proteins').glob('*.pdbqt'))
    ligands = list(Path('data/prepared_ligands').glob('*.pdbqt'))
    rows = []
    seeds = [42, 1337, 2024]
    exhaust = [8, 16]
    for p in proteins[:1]:
        for l in ligands[:5]:
            for s in seeds:
                for ex in exhaust:
                    vs = run_vina(p, l, Path('data/docking'), exhaustiveness=ex, seed=s)
                    ss = run_smina(p, l, Path('data/docking'), exhaustiveness=ex, seed=s)
                    rows.append({
                        'id': l.stem,
                        'vina_score': vs,
                        'smina_score': ss,
                        'consensus': (vs + ss)/2,
                        'seed': s,
                        'exhaustiveness': ex
                    })
    df = pd.DataFrame(rows)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    logger.info(f"Saved variance scores -> {out}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    main(args.out)
