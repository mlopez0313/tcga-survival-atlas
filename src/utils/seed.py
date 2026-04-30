"""Single entry point for reproducibility.

Seeds the python `random`, numpy, and torch (if installed) RNGs.
"""
from __future__ import annotations

import os
import random

import numpy as np


def set_global_seed(seed: int) -> None:
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch

        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        # Determinism trade-off: enabling here would slow training a lot.
        torch.backends.cudnn.deterministic = False
        torch.backends.cudnn.benchmark = True
    except ImportError:
        pass
