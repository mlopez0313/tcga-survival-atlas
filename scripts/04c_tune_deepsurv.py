from __future__ import annotations

import itertools
import json

import _bootstrap  # noqa: F401


def main() -> None:
    cfg, _ = _bootstrap.setup("Step 04c: tune DeepSurv")
    from src.modeling import train_deepsurv
    from src.utils.logging import get_logger

    log = get_logger("deepsurv_tune", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    base_spec = dict(cfg["models"]["deepsurv"])
    base_spec["prefilter_top_k"] = 1500

    hidden_opts = [[128, 64], [64, 32], [32]]
    dropout_opts = [0.3, 0.4, 0.5]
    wd_opts = [1e-4, 5e-4, 1e-3]
    lr_opts = [1e-3, 3e-4, 1e-4]
    bn_opts = [True, False]

    combos = list(itertools.product(hidden_opts, dropout_opts, wd_opts, lr_opts, bn_opts))
    log.info(f"DeepSurv tuning over {len(combos)} combinations")

    rows = []
    best = None
    best_score = float("-inf")

    for i, (hidden, dropout, wd, lr, bn) in enumerate(combos, start=1):
        cfg["models"]["deepsurv"] = dict(base_spec)
        cfg["models"]["deepsurv"]["hidden_layers"] = hidden
        cfg["models"]["deepsurv"]["dropout"] = dropout
        cfg["models"]["deepsurv"]["weight_decay"] = wd
        cfg["models"]["deepsurv"]["lr"] = lr
        cfg["models"]["deepsurv"]["batch_norm"] = bn
        log.info(f"[{i}/{len(combos)}] hidden={hidden} dropout={dropout} wd={wd} lr={lr} bn={bn}")
        out = train_deepsurv(cfg)
        m = out["metrics"]
        row = {
            "hidden_layers": hidden,
            "dropout": dropout,
            "weight_decay": wd,
            "lr": lr,
            "batch_norm": bn,
            "prefilter_top_k": 1500,
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
            f"DeepSurv result hidden={hidden} dropout={dropout} wd={wd} lr={lr} bn={bn} "
            f"train={row['c_index_train']:.4f} test={row['c_index_test']:.4f} gap={row['gap_train_test']:.4f}"
        )

    with open(f"{cfg['paths']['metrics']}/deepsurv_tuning_results.json", "w") as fh:
        json.dump(rows, fh, indent=2)
    with open(f"{cfg['paths']['metrics']}/deepsurv_tuning_best.json", "w") as fh:
        json.dump(best, fh, indent=2)

    log.info(f"Best DeepSurv config: {best}")


if __name__ == "__main__":
    main()
