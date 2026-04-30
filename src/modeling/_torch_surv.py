"""Shared PyTorch utilities for the deep survival models.

Centralizes:
    * Cox partial-likelihood loss (Breslow approximation, vectorized),
    * an early-stopping helper,
    * a generic train loop step.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch


def cox_ph_loss(risk: torch.Tensor, time: torch.Tensor, event: torch.Tensor) -> torch.Tensor:
    """Negative Cox partial log-likelihood with the Breslow approximation.

    risk : (N,)  raw model output (log hazard)
    time : (N,)  survival time
    event: (N,)  1 if event observed, 0 if censored

    For each observed event i:
        contribution = risk_i - log( sum_{j: time_j >= time_i} exp(risk_j) )

    After sorting by *descending* time, the risk set j: t_j >= t_i is the
    prefix [0..i], so the denominator becomes a cumulative log-sum-exp
    from the top — implemented stably by `torch.logcumsumexp`.
    """
    risk = risk.view(-1)
    time = time.view(-1)
    event = event.view(-1).float()

    order = torch.argsort(time, descending=True)
    risk = risk[order]
    event = event[order]

    log_cum_hazard = torch.logcumsumexp(risk, dim=0)
    log_lik = (risk - log_cum_hazard) * event
    n_events = event.sum().clamp(min=1.0)
    return -log_lik.sum() / n_events


@dataclass
class EarlyStopper:
    patience: int
    min_delta: float = 1e-4
    best: float = float("inf")
    counter: int = 0
    should_stop: bool = False

    def step(self, metric: float) -> bool:
        """Pass current loss/-cindex (anything where lower is better)."""
        if metric < self.best - self.min_delta:
            self.best = metric
            self.counter = 0
        else:
            self.counter += 1
            if self.counter >= self.patience:
                self.should_stop = True
        return self.should_stop


def to_tensor(arr: np.ndarray, device: torch.device, dtype=torch.float32) -> torch.Tensor:
    return torch.as_tensor(np.asarray(arr), dtype=dtype, device=device)
