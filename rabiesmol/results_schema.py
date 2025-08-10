from __future__ import annotations
import pyarrow.parquet as pq
from pydantic import BaseModel
class ResultSchema(BaseModel):
 schema_version: str='v1'
 required_columns: list[str]=['id','smiles','vina_score']

def validate_results(parquet_path: str)->None:
 meta=pq.read_metadata(parquet_path)
 cols=set(meta.schema.names)
 missing=set(ResultSchema().required_columns)-cols
 if missing:
  raise ValueError(str(sorted(missing)))
