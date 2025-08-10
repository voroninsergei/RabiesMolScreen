

from pathlib import Path
import yaml
from typing import Any

def export_config(config: dict[str, Any], out_file: str) -> None:
    """Export configuration dictionary to YAML file."""
    Path(out_file).write_text(yaml.safe_dump(config), encoding="utf-8")

def import_config(file_path: str) -> dict[str, Any]:
    """Load configuration from YAML file."""
    return yaml.safe_load(Path(file_path).read_text(encoding="utf-8"))
