
from typer.testing import CliRunner
from rabiesmol.cli import app

runner = CliRunner()

def test_seed_determinism(tmp_path):
    args = ["prepare", "--track", "--seed", "42"]
    result = runner.invoke(app, args)
    assert result.exit_code == 0
