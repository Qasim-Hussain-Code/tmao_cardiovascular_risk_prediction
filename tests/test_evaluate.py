"""Tests for the performance measures.

The area under the curve and its variance are implemented in this package
rather than taken from a library, so they are checked against an independent
implementation and against cases whose answers are known analytically.
"""

from __future__ import annotations

import numpy as np
import pytest
from sklearn.metrics import roc_auc_score

from tmao_cvd.evaluate import (
    auc_with_ci,
    calibration_intercept,
    calibration_slope,
    category_free_nri,
    delong_test,
    integrated_discrimination_improvement,
    net_benefit,
    treat_all_net_benefit,
)


@pytest.fixture
def sample():
    rng = np.random.default_rng(3)
    y = rng.binomial(1, 0.25, 600)
    signal = rng.normal(y * 0.9, 1.0)
    return y, 1.0 / (1.0 + np.exp(-signal))


def test_delong_auc_agrees_with_scikit_learn(sample):
    y, probabilities = sample
    auc, lower, upper = auc_with_ci(y, probabilities)
    assert auc == pytest.approx(roc_auc_score(y, probabilities), abs=1e-9)
    assert lower < auc < upper


def test_delong_handles_tied_predictions():
    # Midranks matter only when there are ties, so the tied case is the one
    # worth asserting on.
    y = np.array([0, 0, 1, 1, 0, 1])
    tied = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
    auc, _, _ = auc_with_ci(y, tied)
    assert auc == pytest.approx(0.5)


def test_identical_predictions_give_no_difference(sample):
    y, probabilities = sample
    difference, p_value = delong_test(y, probabilities, probabilities)
    assert difference == pytest.approx(0.0, abs=1e-12)
    assert p_value == pytest.approx(1.0)


def test_delong_detects_a_better_model(sample):
    y, weak = sample
    strong = 1.0 / (1.0 + np.exp(-(np.random.default_rng(5).normal(y * 2.5, 1.0))))
    difference, p_value = delong_test(y, weak, strong)
    assert difference > 0
    assert p_value < 0.01


def test_perfectly_calibrated_predictions_recover_slope_one():
    # Outcomes drawn from the predicted probabilities themselves should give
    # an intercept near zero and a slope near one.
    rng = np.random.default_rng(17)
    probabilities = rng.uniform(0.02, 0.85, 40000)
    y = rng.binomial(1, probabilities)
    assert calibration_intercept(y, probabilities) == pytest.approx(0.0, abs=0.06)
    assert calibration_slope(y, probabilities) == pytest.approx(1.0, abs=0.06)


def test_idi_is_zero_when_the_models_agree(sample):
    y, probabilities = sample
    assert integrated_discrimination_improvement(y, probabilities, probabilities) == pytest.approx(0.0)


def test_nri_rewards_movement_in_the_right_direction():
    y = np.array([1, 1, 0, 0])
    baseline = np.array([0.4, 0.4, 0.4, 0.4])
    extended = np.array([0.6, 0.6, 0.2, 0.2])
    events, non_events, total = category_free_nri(y, baseline, extended)
    assert events == pytest.approx(1.0)
    assert non_events == pytest.approx(1.0)
    assert total == pytest.approx(2.0)


def test_net_benefit_of_treating_everyone_matches_the_closed_form():
    y = np.array([1, 0, 0, 1, 0, 0, 0, 0, 1, 0])
    thresholds = np.array([0.05, 0.10, 0.25])
    always_treat = np.ones_like(y, dtype=float)
    assert net_benefit(y, always_treat, thresholds) == pytest.approx(
        treat_all_net_benefit(y, thresholds)
    )


def test_net_benefit_is_zero_when_nobody_is_treated():
    y = np.array([1, 0, 0, 1, 0])
    thresholds = np.array([0.1, 0.3])
    assert net_benefit(y, np.zeros_like(y, dtype=float), thresholds) == pytest.approx(
        np.zeros(2)
    )
