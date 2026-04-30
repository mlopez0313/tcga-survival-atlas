"""Train the custom multi-branch survival model."""
from __future__ import annotations

import _bootstrap  # noqa: F401


def main() -> None:
    cfg, _ = _bootstrap.setup("Step 05: train MultiBranch survival model")
    from src.modeling import train_multibranch
    from src.utils.logging import get_logger

    log = get_logger("script_multibranch", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    if not cfg["models"]["multibranch"].get("enabled", True):
        log.info("MultiBranch disabled in config; skipping.")
        return
    train_multibranch(cfg)


if __name__ == "__main__":
    main()
