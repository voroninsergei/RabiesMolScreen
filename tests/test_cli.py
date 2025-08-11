from typer.testing import CliRunner

from rabiesmol.cli.main import app


def test_prepare_help():
    runner = CliRunner()
    result = runner.invoke(app, ["prepare", "--help"])
    assert result.exit_code == 0
    assert "proteins" in result.output
    assert "ligands" in result.output
