"""Figures reporting discrimination, calibration and clinical utility.

The three figures are meant to be read together. A model can look convincing
in any one of them on its own, and the combination is what the reporting
guidance for prediction models asks for (Collins et al., TRIPOD, Annals of
Internal Medicine, 2015).

Every figure carries the provenance of the data it was built from in its
title. When the pipeline has fallen back to simulation this is the only
thing standing between a plausible looking curve and a reader who assumes it
came from a cohort, so it is not decoration.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_curve

from .evaluate import net_benefit, treat_all_net_benefit


def _stamp(figure, provenance: str) -> None:
    """Write the data provenance beneath a figure."""
    figure.text(
        0.5,
        0.005,
        f"Data source: {provenance}",
        ha="center",
        va="bottom",
        fontsize=8,
        style="italic",
    )


def plot_roc_curves(
    y_true: np.ndarray,
    predictions: dict[str, np.ndarray],
    provenance: str,
    path: Path,
    dpi: int,
) -> None:
    """Receiver operating characteristic curves for the compared models."""

    figure, axis = plt.subplots(figsize=(5.5, 5.5))

    for label, probabilities in predictions.items():
        false_positive, true_positive, _ = roc_curve(y_true, probabilities)
        axis.plot(false_positive, true_positive, linewidth=1.6, label=label)

    axis.plot([0, 1], [0, 1], color="grey", linestyle=":", linewidth=1.0)
    axis.set_xlabel("False positive rate")
    axis.set_ylabel("True positive rate")
    axis.set_title("Discrimination of the baseline and extended models")
    axis.legend(loc="lower right", frameon=False, fontsize=9)
    axis.set_aspect("equal")

    _stamp(figure, provenance)
    figure.tight_layout(rect=(0, 0.03, 1, 1))
    figure.savefig(path, dpi=dpi)
    plt.close(figure)


def plot_calibration(
    y_true: np.ndarray,
    predictions: dict[str, np.ndarray],
    provenance: str,
    path: Path,
    dpi: int,
    n_bins: int = 10,
) -> None:
    """Observed event rate against predicted risk, by decile of prediction.

    Equal frequency bins are used rather than equal width bins. Predicted
    risks are concentrated at the low end when the event rate is low, and
    equal width bins would leave the upper bins nearly empty and produce
    points whose scatter is mostly sampling noise.
    """

    figure, axis = plt.subplots(figsize=(5.5, 5.5))

    for label, probabilities in predictions.items():
        quantiles = np.quantile(probabilities, np.linspace(0, 1, n_bins + 1))
        quantiles[-1] += 1e-9
        bin_index = np.digitize(probabilities, quantiles[1:-1])

        observed, expected = [], []
        for index in range(n_bins):
            selected = bin_index == index
            if selected.sum() == 0:
                continue
            observed.append(float(np.mean(y_true[selected])))
            expected.append(float(np.mean(probabilities[selected])))

        axis.plot(expected, observed, marker="o", markersize=4, linewidth=1.4, label=label)

    limit = max(
        max(np.max(values) for values in predictions.values()),
        float(np.mean(y_true)) * 2,
    )
    axis.plot([0, limit], [0, limit], color="grey", linestyle=":", linewidth=1.0)
    axis.set_xlabel("Predicted risk")
    axis.set_ylabel("Observed event rate")
    axis.set_title("Calibration of out of fold predictions")
    axis.legend(loc="upper left", frameon=False, fontsize=9)

    _stamp(figure, provenance)
    figure.tight_layout(rect=(0, 0.03, 1, 1))
    figure.savefig(path, dpi=dpi)
    plt.close(figure)


def plot_decision_curve(
    y_true: np.ndarray,
    predictions: dict[str, np.ndarray],
    thresholds: np.ndarray,
    provenance: str,
    path: Path,
    dpi: int,
) -> None:
    """Net benefit of each model against the two default strategies."""

    figure, axis = plt.subplots(figsize=(6.0, 5.0))

    axis.plot(
        thresholds,
        treat_all_net_benefit(y_true, thresholds),
        color="grey",
        linewidth=1.2,
        label="treat all",
    )
    axis.axhline(0.0, color="black", linewidth=1.0, linestyle="--", label="treat none")

    for label, probabilities in predictions.items():
        axis.plot(
            thresholds,
            net_benefit(y_true, probabilities, thresholds),
            linewidth=1.6,
            label=label,
        )

    axis.set_xlabel("Risk threshold")
    axis.set_ylabel("Net benefit")
    axis.set_title("Decision curve analysis")
    axis.set_ylim(-0.02, float(np.mean(y_true)) * 1.2)
    axis.legend(loc="upper right", frameon=False, fontsize=9)

    _stamp(figure, provenance)
    figure.tight_layout(rect=(0, 0.03, 1, 1))
    figure.savefig(path, dpi=dpi)
    plt.close(figure)
