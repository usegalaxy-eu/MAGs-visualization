"""docs/cli.md is generated from the argument parser and must stay in sync with it."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_cli_reference_is_up_to_date():
    res = subprocess.run([sys.executable, "scripts/build_docs.py", "--check"], cwd=ROOT, capture_output=True,
                         text=True)
    assert res.returncode == 0, res.stdout + res.stderr


def test_gallery_covers_every_subcommand():
    sys.path.insert(0, str(ROOT))
    from mags_visualization.main import build_parser
    import argparse

    sys.path.insert(0, str(ROOT / "scripts"))
    from build_docs import EXAMPLES

    sub = next(a for a in build_parser()._actions if isinstance(a, argparse._SubParsersAction))
    documented = {args[0] for _, _, _, args in EXAMPLES}
    assert set(sub.choices) - {"all"} <= documented
