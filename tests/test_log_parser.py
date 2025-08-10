
from pathlib import Path
from rabiesmol.docking import _parse_vina_log

def test_parse_vina_log_basic(tmp_path: Path):
    log = tmp_path/"vina.log"
    log.write_text("-----\n1   -7.123   0.0  0.0  0.0\n2 -6.1\n", encoding="utf-8")
    assert abs(_parse_vina_log(log) + 7.123) < 1e-6

def test_parse_vina_log_commas(tmp_path: Path):
    log = tmp_path/"vina2.log"
    log.write_text("header\n1,-8.0,0,0,0\n", encoding="utf-8")
    assert _parse_vina_log(log) == -8.0
