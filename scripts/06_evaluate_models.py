"""Aggregate metrics across all trained models into a single comparison table.

Reads `results/metrics/*.json` (one per model) and produces:
    - results/metrics/summary.csv
    - results/figures/cindex_comparison.png
"""
from __future__ import annotations

import json
from pathlib import Path

import _bootstrap  # noqa: F401


def main() -> None:
    cfg, _ = _bootstrap.setup("Step 06: aggregate model metrics")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    from src.utils.logging import get_logger

    log = get_logger("script_evaluate", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))

    metrics_dir = Path(cfg["paths"]["metrics"])
    fig_dir = Path(cfg["paths"]["figures"])

    rows = []
    canonical = {}
    for path in sorted(metrics_dir.glob("*.json")):
        if path.name.startswith("summary"):
            continue
        with open(path) as fh:
            data = json.load(fh)
        if not isinstance(data, dict):
            continue
        if "model" in data:
            canonical[data["model"]] = data

    # Prefer best-tuning artifacts when present.
    tuned_map = {
        "rsf_tuning_prefilter_best.json": {"model": "random_survival_forest", "feature_set": "baseline"},
        "deepsurv_tuning_best.json": {"model": "deepsurv", "feature_set": "baseline"},
        "multibranch_phase2_best.json": {"model": "multibranch"},
    }
    for fn, base in tuned_map.items():
        path = metrics_dir / fn
        if not path.exists():
            continue
        with open(path) as fh:
            data = json.load(fh)
        if not isinstance(data, dict):
            continue
        merged = dict(base)
        merged.update(data)
        canonical[base["model"]] = merged

    rows = list(canonical.values())

    if not rows:
        log.warning("No per-model metrics JSON files found. Run training scripts first.")
        return

    df = pd.DataFrame(rows)
    keep_cols = [c for c in [
        "model", "feature_set", "n_train", "n_test",
        "events_train", "events_test",
        "c_index_train", "c_index_test",
        "km_test_logrank_p", "km_test_logrank_stat",
        "alpha", "l1_ratio", "n_features", "branches",
        "best_val_loss", "best_val_cindex", "epochs_trained",
    ] if c in df.columns]
    df = df[keep_cols].sort_values("c_index_test", ascending=False)
    summary_path = metrics_dir / "summary.csv"
    df.to_csv(summary_path, index=False)
    log.info(f"Wrote summary -> {summary_path}")
    log.info("\n" + df.to_string(index=False))

    # Side-by-side bar of train vs test C-index
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    x = range(len(df))
    width = 0.4
    ax.bar([i - width / 2 for i in x], df["c_index_train"], width, label="train", color="#9ecae1")
    ax.bar([i + width / 2 for i in x], df["c_index_test"], width, label="test", color="#3182bd")
    ax.set_xticks(list(x))
    ax.set_xticklabels(df["model"], rotation=20, ha="right")
    ax.axhline(0.5, ls="--", color="gray", lw=0.8, label="random")
    ax.set_ylabel("Concordance index")
    ax.set_title("Survival model comparison")
    ax.set_ylim(0.4, 1.0)
    ax.legend()
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(fig_dir / "cindex_comparison.png", dpi=150)
    plt.close(fig)
    log.info(f"Wrote comparison figure -> {fig_dir / 'cindex_comparison.png'}")


if __name__ == "__main__":
    main()
