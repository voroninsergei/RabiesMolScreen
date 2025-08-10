
from typer.testing import CliRunner
from rabiesmol.cli import app

runner = CliRunner()

def test_help_runs():
    r = runner.invoke(app, ["--help"])
    assert r.exit_code == 0
    assert "rabiesmol" in r.stdout.lower()

def test_validate_config_missing_file():
    r = runner.invoke(app, ["validate-config", "no_such.yaml"])
    assert r.exit_code != 0
