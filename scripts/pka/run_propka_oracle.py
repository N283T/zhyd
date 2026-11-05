#!/usr/bin/env python3
"""Run an installed PROPKA executable as an external oracle (never imported)."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--executable", default="propka3")
    args = parser.parse_args()
    executable = shutil.which(args.executable)
    if not executable:
        raise SystemExit(f"PROPKA executable not found: {args.executable}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for source in args.inputs:
        with tempfile.TemporaryDirectory(prefix="zhyd-propka-") as directory:
            work = Path(directory)
            local = work / source.name
            shutil.copy2(source, local)
            subprocess.run([executable, local.name], cwd=work, check=True)
            result = work / f"{local.stem}.pka"
            if not result.exists():
                raise SystemExit(f"PROPKA did not create expected output for {source}")
            shutil.copy2(result, args.output_dir / result.name)


if __name__ == "__main__":
    main()
