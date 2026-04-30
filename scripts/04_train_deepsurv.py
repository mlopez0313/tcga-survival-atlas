"""Train the DeepSurv MLP."""
from __future__ import annotations

import _bootstrap  # noqa: F401


def main() -> None:
    cfg, _ = _bootstrap.setup("Step 04: train DeepSurv")
    from src.modeling import train_deepsurv
    from src.utils.logging import get_logger

    log = get_logger("script_deepsurv", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    if not cfg["models"]["deepsurv"].get("enabled", True):
        log.info("DeepSurv disabled in config; skipping.")
        return
    train_deepsurv(cfg)


if __name__ == "__main__":
    main()
