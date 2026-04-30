"""Train Cox elastic-net and RSF baselines on the processed data."""
from __future__ import annotations

import _bootstrap  # noqa: F401


def main() -> None:
    cfg, _ = _bootstrap.setup("Step 03: train Cox elastic-net + RSF baselines")
    from src.modeling import train_cox_elastic_net, train_rsf
    from src.utils.logging import get_logger

    log = get_logger("script_baselines", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))

    if cfg["models"]["cox_elastic_net"].get("enabled", True):
        log.info("=== Cox elastic-net ===")
        train_cox_elastic_net(cfg)
    else:
        log.info("Cox elastic-net disabled in config; skipping.")

    if cfg["models"]["rsf"].get("enabled", True):
        log.info("=== Random Survival Forest ===")
        train_rsf(cfg)
    else:
        log.info("RSF disabled in config; skipping.")


if __name__ == "__main__":
    main()
