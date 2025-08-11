from rabiesmol.results_schema import SCHEMA_VERSION, RESULT_COLUMNS
from rabiesmol.domain.models import RESULT_SCHEMA_VERSION as MODEL_VERSION, RESULT_COLUMNS as MODEL_COLUMNS


def test_schema_version_consistency():
    assert SCHEMA_VERSION == MODEL_VERSION


def test_result_columns_consistency():
    assert set(RESULT_COLUMNS) == set(MODEL_COLUMNS.keys())
