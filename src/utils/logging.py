"""Project-wide logger factory.

Logs go to stdout AND to logs/<name>.log so each script run leaves a
trail next to the data and metrics.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional

_FMT = "%(asctime)s | %(levelname)-7s | %(name)s | %(message)s"


def get_logger(
    name: str,
    log_dir: Optional[str | Path] = None,
    level: str = "INFO",
) -> logging.Logger:
    """Return a configured logger that writes to stdout and (optionally) a file."""
    logger = logging.getLogger(name)
    if logger.handlers:  # already configured by an earlier call this process
        return logger

    logger.setLevel(level)
    formatter = logging.Formatter(_FMT)

    stream = logging.StreamHandler(sys.stdout)
    stream.setFormatter(formatter)
    logger.addHandler(stream)

    if log_dir is not None:
        log_dir = Path(log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_dir / f"{name}.log")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    logger.propagate = False
    return logger
