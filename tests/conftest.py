
import os
import io
from pathlib import Path
import builtins
import types
import pytest

class FakeCompletedProcess:
    def __init__(self, returncode=0):
        self.returncode = returncode

@pytest.fixture(autouse=True)
def mock_subprocess_run(monkeypatch, tmp_path):
    def fake_run(cmd, check=False, **kwargs):
        # Normalize command (list or str)
        if isinstance(cmd, (list, tuple)):
            exe = str(cmd[0])
            args = list(cmd[1:])
        else:
            parts = str(cmd).split()
            exe, args = parts[0], parts[1:]

        def ensure_parent(p: Path):
            p.parent.mkdir(parents=True, exist_ok=True)

        # Simulate external tools by filename
        if "prepare_receptor" in exe:
            # expect "-r <in> -o <out>"
            out_path = Path(args[args.index("-o")+1])
            ensure_parent(out_path)
            out_path.write_text("RECEPTOR PDBQT")
            return FakeCompletedProcess(0)

        if "prepare_ligand" in exe:
            out_path = Path(args[args.index("-o")+1])
            ensure_parent(out_path)
            out_path.write_text("LIGAND PDBQT")
            return FakeCompletedProcess(0)

        if "obabel" in exe:
            # "-O <output>"
            out_path = Path(args[args.index("-O")+1])
            ensure_parent(out_path)
            out_path.write_text("OBABEL OUTPUT")
            return FakeCompletedProcess(0)

        if "vina" in exe or "smina" in exe:
            # write both out pdbqt and a log file if specified
            if "--out" in args:
                out_path = Path(args[args.index("--out")+1])
                ensure_parent(out_path)
                out_path.write_text("DOCKED POSE")
            if "--log" in args:
                log_path = Path(args[args.index("--log")+1])
                ensure_parent(log_path)
                # minimal log with first mode and score in 2nd column
                log_path.write_text("   1   -7.5   0.0   0.0\n")
            return FakeCompletedProcess(0)

        if "fpocket" in exe:
            # create expected output structure
            # <protein>_out/pockets/pocket0/pocket.pdb
            # We don't have the input file path from args easily, but fpocket uses "-f <file>"
            if "-f" in args:
                inp = Path(args[args.index("-f")+1])
                out_dir = inp.with_suffix("").name + "_out"
                pockets_dir = Path(out_dir) / "pockets" / "pocket0"
                pockets_dir.mkdir(parents=True, exist_ok=True)
                (pockets_dir / "pocket.pdb").write_text("POCKET")
            return FakeCompletedProcess(0)

        # default: pretend success
        return FakeCompletedProcess(0)

    monkeypatch.setattr("subprocess.run", fake_run)
    yield
