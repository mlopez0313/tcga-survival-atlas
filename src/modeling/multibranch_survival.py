"""Multi-branch deep survival model.

Each modality is fed through its own MLP encoder, the embeddings are
concatenated, and a fusion MLP outputs a single risk score per patient.
The Cox partial-likelihood loss is then applied as in DeepSurv.

Modalities are optional. Any modality not present in the processed
data (or removed from `cfg.models.multibranch.modalities`) is skipped.
A ready-to-train model needs at least the clinical branch.
"""
from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from torch import nn

from ..evaluation.metrics import concordance_index, summarize_metrics
from ..evaluation.survival_curves import plot_km_by_risk
from ..utils.io import load_parquet, save_pickle
from ..utils.logging import get_logger
from ._data import load_processed
from ._torch_surv import EarlyStopper, cox_ph_loss, to_tensor

# Map config name -> processed parquet key
_MODALITY_TO_KEY = {
    "clinical": "X_clinical",
    "rna": "X_rna",
    "mutation": "X_mutation",
    "cnv": "X_cnv",
    "methylation": "X_methylation",
    "mirna": "X_mirna",
    "histology": "X_histology",  # only if a future preprocessing step makes it
}


def _make_branch(in_features: int, hidden: List[int], embedding: int,
                 dropout: float, batch_norm: bool = True) -> nn.Sequential:
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
    layers.append(nn.Linear(prev, embedding))
    layers.append(nn.ReLU(inplace=True))
    return nn.Sequential(*layers)


class MultiBranchSurvivalModel(nn.Module):
    """Per-modality encoders + fusion MLP -> scalar risk."""

    def __init__(self,
                 input_dims: Dict[str, int],
                 branch_hidden: List[int],
                 branch_embedding: int,
                 fusion_hidden: List[int],
                 dropout: float = 0.3,
                 batch_norm: bool = True):
        super().__init__()
        self.modality_names = list(input_dims.keys())
        self.encoders = nn.ModuleDict({
            name: _make_branch(d, branch_hidden, branch_embedding, dropout, batch_norm)
            for name, d in input_dims.items()
        })

        fusion_in = branch_embedding * len(input_dims)
        layers: list[nn.Module] = []
        prev = fusion_in
        for h in fusion_hidden:
            layers.append(nn.Linear(prev, h))
            if batch_norm:
                layers.append(nn.BatchNorm1d(h))
            layers.append(nn.ReLU(inplace=True))
            if dropout > 0:
                layers.append(nn.Dropout(dropout))
            prev = h
        layers.append(nn.Linear(prev, 1))
        self.fusion = nn.Sequential(*layers)

    def forward(self, batch: Dict[str, torch.Tensor]) -> torch.Tensor:
        embeddings = [self.encoders[name](batch[name]) for name in self.modality_names]
        z = torch.cat(embeddings, dim=1)
        return self.fusion(z).squeeze(-1)


def _gather_modalities(processed: Dict[str, Any], wanted: List[str], log) -> Dict[str, pd.DataFrame]:
    out: Dict[str, pd.DataFrame] = {}
    for name in wanted:
        key = _MODALITY_TO_KEY.get(name)
        if key and key in processed:
            out[name] = processed[key]
        else:
            log.warning(f"Modality '{name}' requested but missing from processed data; skipping.")
    return out


def _predict(model: MultiBranchSurvivalModel,
             tensors: Dict[str, torch.Tensor]) -> np.ndarray:
    model.eval()
    with torch.no_grad():
        return model(tensors).cpu().numpy()


def _build_tensors(modalities: Dict[str, pd.DataFrame],
                   ids: List[str],
                   device: torch.device) -> Dict[str, torch.Tensor]:
    return {
        name: to_tensor(df.loc[ids].values.astype(np.float32), device)
        for name, df in modalities.items()
    }


def train_multibranch(cfg: Dict[str, Any]) -> Dict[str, Any]:
    log = get_logger("multibranch", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    spec = cfg["models"]["multibranch"]
    seed = int(cfg.get("seed", 42))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log.info(f"Using device: {device}")

    processed = load_processed(cfg)
    y = processed["y"]
    split = processed["split"]

    modalities = _gather_modalities(processed, list(spec.get("modalities", [])), log)
    if "clinical" not in modalities:
        raise RuntimeError("Multi-branch model needs at least the clinical modality.")

    branch_prefilter = spec.get("branch_prefilter_top_k", {}) or {}
    for name, top_k in branch_prefilter.items():
        if name in modalities and top_k is not None and int(top_k) > 0 and modalities[name].shape[1] > int(top_k):
            keep = modalities[name].var(axis=0).sort_values(ascending=False).head(int(top_k)).index
            modalities[name] = modalities[name].loc[:, keep]
            log.info(f"Branch prefilter for {name}: keeping top-{int(top_k)} variance features")

    log.info(f"Active branches: {list(modalities.keys())}")

    # Patients present in every selected modality (should already be aligned)
    common_ids = list(y.index)
    for df in modalities.values():
        common_ids = [p for p in common_ids if p in df.index]

    train_ids = [p for p in split["train"] if p in common_ids]
    test_ids = [p for p in split["test"] if p in common_ids]

    val_size = float(spec.get("val_size", 0.2))
    inner_train_ids, val_ids = train_test_split(
        train_ids, test_size=val_size, random_state=seed,
        stratify=y.loc[train_ids, "OS_event"].astype(int).values,
    )

    # Per-branch scaler fit on inner-train only
    scalers: Dict[str, StandardScaler] = {}
    scaled: Dict[str, pd.DataFrame] = {}
    for name, df in modalities.items():
        scaler = StandardScaler()
        scaler.fit(df.loc[inner_train_ids].values)
        scaled[name] = pd.DataFrame(
            scaler.transform(df.values).astype(np.float32),
            index=df.index, columns=df.columns,
        )
        scalers[name] = scaler

    input_dims = {n: scaled[n].shape[1] for n in scaled}
    log.info(f"Input dims per branch: {input_dims}")

    model = MultiBranchSurvivalModel(
        input_dims=input_dims,
        branch_hidden=list(spec.get("branch_hidden", [64, 32])),
        branch_embedding=int(spec.get("branch_embedding", 16)),
        fusion_hidden=list(spec.get("fusion_hidden", [64, 32])),
        dropout=float(spec.get("dropout", 0.3)),
        batch_norm=bool(spec.get("batch_norm", True)),
    ).to(device)

    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=float(spec.get("lr", 1e-3)),
        weight_decay=float(spec.get("weight_decay", 1e-4)),
    )

    # Pre-build tensors
    tr_tensors = _build_tensors(scaled, inner_train_ids, device)
    va_tensors = _build_tensors(scaled, val_ids, device)
    te_tensors = _build_tensors(scaled, test_ids, device)
    full_train_tensors = _build_tensors(scaled, train_ids, device)

    ttr = to_tensor(y.loc[inner_train_ids, "OS_time"].values, device)
    etr = to_tensor(y.loc[inner_train_ids, "OS_event"].values, device)
    tva = to_tensor(y.loc[val_ids, "OS_time"].values, device)
    eva = to_tensor(y.loc[val_ids, "OS_event"].values, device)

    epochs = int(spec.get("epochs", 300))
    batch_size = int(spec.get("batch_size", 256))
    stopper = EarlyStopper(patience=int(spec.get("patience", 30)))

    history = {"train_loss": [], "val_loss": [], "val_cindex": []}
    best_state = None
    best_val_cindex = -np.inf

    n_train = len(inner_train_ids)
    use_minibatch = n_train > batch_size
    rng = np.random.default_rng(seed)

    for epoch in range(1, epochs + 1):
        model.train()
        if use_minibatch:
            idx = rng.permutation(n_train)
            losses = []
            for start in range(0, n_train, batch_size):
                sel = idx[start: start + batch_size]
                if etr[sel].sum().item() < 1:
                    continue
                optimizer.zero_grad()
                risk = model({k: v[sel] for k, v in tr_tensors.items()})
                loss = cox_ph_loss(risk, ttr[sel], etr[sel])
                loss.backward()
                optimizer.step()
                losses.append(loss.item())
            train_loss = float(np.mean(losses)) if losses else float("nan")
        else:
            optimizer.zero_grad()
            risk = model(tr_tensors)
            loss = cox_ph_loss(risk, ttr, etr)
            loss.backward()
            optimizer.step()
            train_loss = float(loss.item())

        model.eval()
        with torch.no_grad():
            val_risk = model(va_tensors)
            val_loss = float(cox_ph_loss(val_risk, tva, eva).item())
            val_cindex = concordance_index(
                y.loc[val_ids, "OS_time"].values,
                y.loc[val_ids, "OS_event"].values,
                val_risk.cpu().numpy(),
            )
        history["train_loss"].append(train_loss)
        history["val_loss"].append(val_loss)
        history["val_cindex"].append(val_cindex)

        # Track best by validation C-index (the actual eval metric)
        if val_cindex > best_val_cindex:
            best_val_cindex = val_cindex
            best_state = deepcopy(model.state_dict())

        if epoch == 1 or epoch % 10 == 0:
            log.info(f"epoch {epoch:3d} | train_loss={train_loss:.4f} "
                     f"val_loss={val_loss:.4f} val_C={val_cindex:.4f}")

        # Stop on val loss as a stable signal
        if stopper.step(val_loss):
            log.info(f"Early stopping at epoch {epoch} (best val C={best_val_cindex:.4f})")
            break

    if best_state is not None:
        model.load_state_dict(best_state)

    risk_train = _predict(model, full_train_tensors)
    risk_test = _predict(model, te_tensors)

    y_train = y.loc[train_ids]
    y_test = y.loc[test_ids]

    metrics = summarize_metrics(
        "multibranch", y_train, risk_train, y_test, risk_test,
        extra={"branches": list(modalities.keys()),
               "best_val_cindex": float(best_val_cindex),
               "epochs_trained": len(history["train_loss"])},
    )
    log.info(f"MultiBranch  C-index train={metrics['c_index_train']:.4f}  test={metrics['c_index_test']:.4f}")

    # Persist
    models_dir = Path(cfg["paths"]["models"])
    fig_dir = Path(cfg["paths"]["figures"])
    metrics_dir = Path(cfg["paths"]["metrics"])

    torch.save(model.state_dict(), models_dir / "multibranch.pt")
    save_pickle(
        {"scalers": scalers, "input_dims": input_dims, "config": spec,
         "feature_names": {n: list(scaled[n].columns) for n in scaled}},
        models_dir / "multibranch_meta.pkl",
    )

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].plot(history["train_loss"], label="train")
    axes[0].plot(history["val_loss"], label="val")
    axes[0].set_title("MultiBranch loss")
    axes[0].set_xlabel("epoch")
    axes[0].set_ylabel("Cox loss")
    axes[0].legend()
    axes[0].grid(alpha=0.25)
    axes[1].plot(history["val_cindex"], color="#2ca02c")
    axes[1].set_title("MultiBranch val C-index")
    axes[1].set_xlabel("epoch")
    axes[1].set_ylabel("C-index")
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(fig_dir / "multibranch_training.png", dpi=150)
    plt.close(fig)

    km = plot_km_by_risk(
        y_test["OS_time"].to_numpy(), y_test["OS_event"].to_numpy(),
        risk_test, "MultiBranch (test)", fig_dir / "km_multibranch.png",
    )
    metrics.update({f"km_test_{k}": v for k, v in km.items()})

    with open(metrics_dir / "multibranch.json", "w") as fh:
        json.dump(metrics, fh, indent=2)

    return {
        "model": model, "scalers": scalers, "metrics": metrics, "history": history,
        "risk_train": risk_train, "risk_test": risk_test,
        "y_train": y_train, "y_test": y_test,
    }
