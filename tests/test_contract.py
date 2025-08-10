import pandas as pd
from rabiesmol.io.contract import validate_df, attach_schema_meta
from rabiesmol.domain.models import RESULT_COLUMNS

def test_contract_smoke():
    # minimal DF with required cols
    df = pd.DataFrame({
        "ligand_id":["L1"],
        "protein_id":["P1"],
        "engine":["vina"],
        "pose_id":[1],
        "score":[-7.5],
        "unit":["kcal/mol"],
        "ligand_smiles":[""],
        "metadata":["{}"],
    })
    df = attach_schema_meta(df)
    problems = validate_df(df)
    assert not problems, problems
