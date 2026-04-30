"""Make `src.*` imports work when scripts are run directly.

Each script does `from _bootstrap import setup`. Importing the file
inserts the repo root onto sys.path; calling setup() returns a
loaded config and parsed CLI args.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))


def parse_args(description: str) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--config",
        default=str(_ROOT / "config.yaml"),
        help="Path to config.yaml",
    )
    return parser.parse_args()


def setup(description: str):
    """Parse CLI, load config, ensure folders, set seed, return (cfg, args)."""
    args = parse_args(description)
    from src.utils import ensure_dirs, load_config, set_global_seed  # noqa: E402

    cfg = load_config(args.config)

    # Resolve relative paths against the config file's directory, so the
    # pipeline is invariant to where the user calls the script from.
    config_dir = Path(args.config).resolve().parent
    cfg["paths"] = {
        k: str((config_dir / Path(v).expanduser()).resolve()) if not Path(v).expanduser().is_absolute() else str(Path(v).expanduser().resolve())
        for k, v in cfg["paths"].items()
    }

    ensure_dirs(cfg)
    set_global_seed(int(cfg.get("seed", 42)))
    return cfg, args
