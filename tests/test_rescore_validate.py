
from pathlib import Path
import pandas as pd
from rabiesmol.rescoring import rescore, rescore_to_parquet
from rabiesmol.validation import validate_hits
import pyarrow  # ensure parquet engine available

def test_rescore_validate_parquet(tmp_path: Path):
    # Create docking results csv
    csv = tmp_path / "results.csv"
    df = pd.DataFrame({
        "id": ["a","b"],
        "smiles": ["CCO", "c1ccccc1"],
        "vina_score": [-7.5, -6.0],
        "smina_score": [-7.4, -6.1],
        "consensus": [-7.45, -6.05]
    })
    df.to_csv(csv, index=False)

    out_csv = tmp_path / "rescored.csv"
    rescore(csv, out_csv)
    assert out_csv.exists()

    out_parquet = tmp_path / "results.parquet"
    rescore_to_parquet(csv, out_parquet)
    assert out_parquet.exists()

    # validation augments csv
    validate_out = tmp_path / "validated.csv"
    validate_hits(out_csv, validate_out)
    assert validate_out.exists()
