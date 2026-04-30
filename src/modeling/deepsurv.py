"""DeepSurv: an MLP trained with Cox partial-likelihood loss."""
from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from sklearn.model_selection import train_test_split
from torch import nn

from ..evaluation.metrics import concordance_index, summarize_metrics
from ..evaluation.survival_curves import plot_km_by_risk
from ..utils.io import save_pickle
from ..utils.logging import get_logger
from ._data import fit_scaler_train, load_processed, select_feature_set, split_xy
from ._torch_surv import EarlyStopper, cox_ph_loss, to_tensor


class DeepSurvNet(nn.Module):
    """A standard DeepSurv MLP: Linear -> (BN) -> ReLU -> Dropout, repeated."""

    def __init__(self, in_features: int, hidden: List[int],
                 dropout: float = 0.3, batch_norm: bool = True):
        super().__init__()
        layers: list[nn.Module] = []
        prev = in_features
        for h in hidden:
            layers.append(nn.Linear(prev, h))
            if batch_norm:
                layers.append(nn.BatchNorm1d(h))
            layers.append(nn.ReLU(inplace=True))
            if dropout > 0:
                layers.append(nn.Dropout(dropout))
            prev = h
        layers.append(nn.Linear(prev, 1))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x).squeeze(-1)


def _predict(model: nn.Module, X: np.ndarray, device: torch.device) -> np.ndarray:
    model.eval()
    with torch.no_grad():
        risk = model(to_tensor(X, device)).cpu().numpy()
    return risk


def train_deepsurv(cfg: Dict[str, Any]) -> Dict[str, Any]:
    log = get_logger("deepsurv", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    spec = cfg["models"]["deepsurv"]
    seed = int(cfg.get("seed", 42))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log.info(f"Using device: {device}")

    processed = load_processed(cfg)
    X = select_feature_set(processed, spec.get("feature_set", "baseline"))
    y = processed["y"]
    X_tr, y_tr, X_te, y_te = split_xy(X, y, processed["split"])
    prefilter_top_k = spec.get("prefilter_top_k")
    if prefilter_top_k is not None and int(prefilter_top_k) > 0 and X_tr.shape[1] > int(prefilter_top_k):
        keep = X_tr.var(axis=0).sort_values(ascending=False).head(int(prefilter_top_k)).index
        X_tr = X_tr.loc[:, keep]
        X_te = X_te.loc[:, keep]
        log.info(f"DeepSurv prefilter: keeping top-{int(prefilter_top_k)} variance features")
    X_tr_s, X_te_s, scaler = fit_scaler_train(X_tr, X_te)

    # Train/val split inside train (stratified on event)
    val_size = float(spec.get("val_size", 0.2))
    train_ids, val_ids = train_test_split(
        X_tr_s.index.tolist(), test_size=val_size, random_state=seed,
        stratify=y_tr["OS_event"].astype(int).values,
    )
    Xtr = X_tr_s.loc[train_ids].values.astype(np.float32)
    Xva = X_tr_s.loc[val_ids].values.astype(np.float32)
    ytr_t = y_tr.loc[train_ids, "OS_time"].values.astype(np.float32)
    ytr_e = y_tr.loc[train_ids, "OS_event"].values.astype(np.float32)
    yva_t = y_tr.loc[val_ids, "OS_time"].values.astype(np.float32)
    yva_e = y_tr.loc[val_ids, "OS_event"].values.astype(np.float32)

    log.info(f"DeepSurv shapes: train={Xtr.shape} val={Xva.shape} test={X_te_s.shape}")

    model = DeepSurvNet(
        in_features=Xtr.shape[1],
        hidden=list(spec.get("hidden_layers", [128, 64])),
        dropout=float(spec.get("dropout", 0.3)),
        batch_norm=bool(spec.get("batch_norm", True)),
    ).to(device)

    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=float(spec.get("lr", 1e-3)),
        weight_decay=float(spec.get("weight_decay", 1e-4)),
    )

    epochs = int(spec.get("epochs", 200))
    batch_size = int(spec.get("batch_size", 256))
    stopper = EarlyStopper(patience=int(spec.get("patience", 25)))

    Xtr_t = to_tensor(Xtr, device)
    ttr_t = to_tensor(ytr_t, device)
    etr_t = to_tensor(ytr_e, device)
    Xva_t = to_tensor(Xva, device)
    tva_t = to_tensor(yva_t, device)
    eva_t = to_tensor(yva_e, device)

    history = {"train_loss": [], "val_loss": [], "val_cindex": []}
    best_state = None
    best_val_loss = float("inf")

    n_train = Xtr.shape[0]
    use_minibatch = n_train > batch_size
    rng = np.random.default_rng(seed)

    for epoch in range(1, epochs + 1):
        model.train()
        if use_minibatch:
            idx = rng.permutation(n_train)
            losses = []
            for start in range(0, n_train, batch_size):
                sel = idx[start: start + batch_size]
                # Cox loss needs at least one event in the batch
                if etr_t[sel].sum().item() < 1:
                    continue
                optimizer.zero_grad()
                risk = model(Xtr_t[sel])
                loss = cox_ph_loss(risk, ttr_t[sel], etr_t[sel])
                loss.backward()
                optimizer.step()
                losses.append(loss.item())
            train_loss = float(np.mean(losses)) if losses else float("nan")
        else:
            optimizer.zero_grad()
            risk = model(Xtr_t)
            loss = cox_ph_loss(risk, ttr_t, etr_t)
            loss.backward()
            optimizer.step()
            train_loss = float(loss.item())

        # Validation
        model.eval()
        with torch.no_grad():
            val_risk = model(Xva_t)
            val_loss = float(cox_ph_loss(val_risk, tva_t, eva_t).item())
            val_cindex = concordance_index(yva_t, yva_e, val_risk.cpu().numpy())

        history["train_loss"].append(train_loss)
        history["val_loss"].append(val_loss)
        history["val_cindex"].append(val_cindex)

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state = deepcopy(model.state_dict())

        if epoch == 1 or epoch % 10 == 0:
            log.info(f"epoch {epoch:3d} | train_loss={train_loss:.4f} "
                     f"val_loss={val_loss:.4f} val_C={val_cindex:.4f}")

        if stopper.step(val_loss):
            log.info(f"Early stopping at epoch {epoch} (best val loss {best_val_loss:.4f})")
            break

    if best_state is not None:
        model.load_state_dict(best_state)

    risk_train = _predict(model, X_tr_s.values.astype(np.float32), device)
    risk_test = _predict(model, X_te_s.values.astype(np.float32), device)

    metrics = summarize_metrics(
        "deepsurv", y_tr, risk_train, y_te, risk_test,
        extra={"feature_set": spec.get("feature_set", "baseline"),
               "n_features": int(X_tr.shape[1]),
               "prefilter_top_k": None if prefilter_top_k is None else int(prefilter_top_k),
               "best_val_loss": best_val_loss,
               "epochs_trained": len(history["train_loss"])},
    )
    log.info(f"DeepSurv   C-index train={metrics['c_index_train']:.4f}  test={metrics['c_index_test']:.4f}")

    # Persist
    models_dir = Path(cfg["paths"]["models"])
    fig_dir = Path(cfg["paths"]["figures"])
    metrics_dir = Path(cfg["paths"]["metrics"])

    torch.save(model.state_dict(), models_dir / "deepsurv.pt")
    save_pickle(
        {"scaler": scaler, "feature_names": list(X_tr.columns),
         "config": spec, "feature_set": spec.get("feature_set", "baseline"),
         "in_features": Xtr.shape[1]},
        models_dir / "deepsurv_meta.pkl",
    )

    # Training curves
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].plot(history["train_loss"], label="train")
    axes[0].plot(history["val_loss"], label="val")
    axes[0].set_xlabel("epoch")
    axes[0].set_ylabel("Cox loss")
    axes[0].set_title("DeepSurv loss")
    axes[0].legend()
    axes[0].grid(alpha=0.25)
    axes[1].plot(history["val_cindex"], color="#2ca02c")
    axes[1].set_xlabel("epoch")
    axes[1].set_ylabel("Validation C-index")
    axes[1].set_title("DeepSurv val C-index")
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(fig_dir / "deepsurv_training.png", dpi=150)
    plt.close(fig)

    km = plot_km_by_risk(
        y_te["OS_time"].to_numpy(), y_te["OS_event"].to_numpy(),
        risk_test, "DeepSurv (test)", fig_dir / "km_deepsurv.png",
    )
    metrics.update({f"km_test_{k}": v for k, v in km.items()})

    with open(metrics_dir / "deepsurv.json", "w") as fh:
        json.dump(metrics, fh, indent=2)

    return {
        "model": model, "scaler": scaler, "metrics": metrics, "history": history,
        "risk_train": risk_train, "risk_test": risk_test,
        "X_train": X_tr_s, "X_test": X_te_s, "y_train": y_tr, "y_test": y_te,
    }
