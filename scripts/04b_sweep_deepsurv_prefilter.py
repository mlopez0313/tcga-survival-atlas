from __future__ import annotations

import json

import _bootstrap  # noqa: F401


def main() -> None:
    cfg, _ = _bootstrap.setup("Step 04b: sweep DeepSurv prefilter_top_k")
    from src.modeling import train_deepsurv
    from src.utils.logging import get_logger

    log = get_logger("deepsurv_sweep", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    base_spec = dict(cfg["models"]["deepsurv"])
    ks = [None, 250, 500, 1000, 1500, 2000]
    rows = []
    best = None
    best_score = float("-inf")

    for k in ks:
        cfg["models"]["deepsurv"] = dict(base_spec)
        cfg["models"]["deepsurv"]["prefilter_top_k"] = k
        log.info(f"Running DeepSurv with prefilter_top_k={k}")
        out = train_deepsurv(cfg)
        m = out["metrics"]
        row = {
            "prefilter_top_k": k,
            "c_index_train": m["c_index_train"],
            "c_index_test": m["c_index_test"],
            "gap_train_test": m["c_index_train"] - m["c_index_test"],
            "best_val_loss": m["best_val_loss"],
            "epochs_trained": m["epochs_trained"],
        }
        rows.append(row)
        if row["c_index_test"] > best_score:
            best_score = row["c_index_test"]
            best = row
        log.info(
            f"DeepSurv k={k} train={row['c_index_train']:.4f} "
            f"test={row['c_index_test']:.4f} gap={row['gap_train_test']:.4f}"
        )

    with open(f"{cfg['paths']['metrics']}/deepsurv_prefilter_sweep.json", "w") as fh:
        json.dump(rows, fh, indent=2)
    with open(f"{cfg['paths']['metrics']}/deepsurv_prefilter_best.json", "w") as fh:
        json.dump(best, fh, indent=2)

    log.info(f"Best DeepSurv prefilter result: {best}")


if __name__ == "__main__":
    main()
