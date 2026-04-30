"""Tiny helper: convert a parquet file to CSV on stdout (or a path).

Used by prepare_app_data.R when no R parquet reader is available. We keep
the index column as a regular column named via the parquet metadata
(typically `patient_id`), so R can re-attach it on load.

Usage:
    python _pq2csv.py <input.parquet> [<output.csv>]
If <output.csv> is omitted, writes CSV to stdout.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: _pq2csv.py <input.parquet> [<output.csv>]", file=sys.stderr)
        return 2
    src = Path(sys.argv[1])
    if not src.exists():
        print(f"missing parquet: {src}", file=sys.stderr)
        return 1
    df = pd.read_parquet(src)
    df = df.reset_index()
    out = sys.argv[2] if len(sys.argv) >= 3 else "-"
    if out == "-":
        df.to_csv(sys.stdout, index=False)
    else:
        df.to_csv(out, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
