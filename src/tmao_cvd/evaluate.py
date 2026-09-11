"""Performance measures for binary risk prediction models.

A single summary statistic cannot establish that a biomarker is useful. The
measures collected here follow the framework set out by Steyerberg and
colleagues (Epidemiology, 2010), which separates three questions that are
often run together:

1. Discrimination. Does the model rank people who go on to have an event
   above those who do not? Measured by the area under the receiver
   operating characteristic curve, with variance from the method of DeLong
   and colleagues (Biometrics, 1988).
2. Calibration. Do the predicted risks match the observed frequencies? A
   model can discriminate well and still be systematically wrong about
   absolute risk, which is the quantity a clinician acts on. Measured by the
   calibration intercept and slope, and by the Brier score.
3. Clinical utility. Would using the model lead to better decisions than
   treating everyone or nobody? Measured by net benefit across a range of
   decision thresholds (Vickers and Elkin, Medical Decision Making, 2006).

Reclassification measures sit alongside these. The integrated discrimination
improvement and the category free net reclassification improvement are
reported because they are conventional in the biomarker literature, but they
are known to be sensitive to miscalibration and to reward a model that is
merely more extreme in its predictions. They are read here as descriptive
companions to the decision curve, not as the primary evidence.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, logit
from scipy.stats import norm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score

EPSILON = 1e-9


def _clip_probabilities(probabilities: np.ndarray) -> np.ndarray:
    """Keep probabilities away from zero and one before a logit transform."""
    return np.clip(np.asarray(probabilities, dtype=float), EPSILON, 1.0 - EPSILON)


# ---------------------------------------------------------------------------
# Discrimination
# ---------------------------------------------------------------------------


def _midrank(values: np.ndarray) -> np.ndarray:
    """Midranks of ``values``, averaging the ranks within a tie group.

    Ties matter. Two models applied to the same cohort routinely produce
    tied predictions, and ranking them arbitrarily biases the estimated
    covariance between their curves.
    """

    order = np.argsort(values)
    sorted_values = values[order]
    n = len(values)
    ranks = np.zeros(n, dtype=float)

    index = 0
    while index < n:
        stop = index
        while stop < n and sorted_values[stop] == sorted_values[index]:
            stop += 1
        ranks[index:stop] = 0.5 * (index + stop - 1) + 1
        index = stop

    out = np.empty(n, dtype=float)
    out[order] = ranks
    return out


def _fast_delong(sorted_predictions: np.ndarray, n_positive: int):
    """Areas under the curve and their covariance matrix.

    Implements the O(n log n) formulation of DeLong's estimator given by Sun
    and Xu (IEEE Signal Processing Letters, 2014). The original 1988
    presentation is O(n squared) in the number of observations, which is
    avoidable and becomes awkward inside a bootstrap.

    Parameters
    ----------
    sorted_predictions:
        Array of shape ``(n_models, n_observations)`` whose columns are
        ordered so that all positive cases come first.
    n_positive:
        Number of positive cases.
    """

    m = n_positive
    n = sorted_predictions.shape[1] - m
    if m == 0 or n == 0:
        raise ValueError("both outcome classes must be present")

    positives = sorted_predictions[:, :m]
    negatives = sorted_predictions[:, m:]
    n_models = sorted_predictions.shape[0]

    tx = np.empty((n_models, m), dtype=float)
    ty = np.empty((n_models, n), dtype=float)
    tz = np.empty((n_models, m + n), dtype=float)
    for row in range(n_models):
        tx[row] = _midrank(positives[row])
        ty[row] = _midrank(negatives[row])
        tz[row] = _midrank(sorted_predictions[row])

    aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / (2.0 * n)
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m

    sx = np.atleast_2d(np.cov(v01))
    sy = np.atleast_2d(np.cov(v10))
    covariance = sx / m + sy / n
    return aucs, covariance


def _order_by_outcome(y_true: np.ndarray) -> tuple[np.ndarray, int]:
    """Return an ordering that puts positives first, and the positive count."""
    y_true = np.asarray(y_true, dtype=int)
    order = np.argsort(-y_true, kind="mergesort")
    return order, int(y_true.sum())


def auc_with_ci(
    y_true: np.ndarray, predictions: np.ndarray, alpha: float = 0.05
) -> tuple[float, float, float]:
    """Area under the curve with a DeLong confidence interval.

    The interval uses the normal approximation on the area scale, as in the
    original paper. It is adequate in the middle of the range and becomes
    conservative as the area approaches one, where the sampling distribution
    is visibly skewed.
    """

    order, n_positive = _order_by_outcome(y_true)
    stacked = np.asarray(predictions, dtype=float).reshape(1, -1)[:, order]
    aucs, covariance = _fast_delong(stacked, n_positive)

    auc = float(aucs[0])
    standard_error = float(np.sqrt(covariance[0, 0]))
    half_width = norm.ppf(1.0 - alpha / 2.0) * standard_error
    return auc, max(0.0, auc - half_width), min(1.0, auc + half_width)


def delong_test(
    y_true: np.ndarray, predictions_a: np.ndarray, predictions_b: np.ndarray
) -> tuple[float, float]:
    """Compare two correlated curves on the same cohort.

    Returns the difference in area (``b`` minus ``a``) and a two sided p
    value. The test accounts for the correlation between the two sets of
    predictions, which is what distinguishes it from comparing two
    independent confidence intervals by eye.
    """

    order, n_positive = _order_by_outcome(y_true)
    stacked = np.vstack(
        [np.asarray(predictions_a, dtype=float), np.asarray(predictions_b, dtype=float)]
    )[:, order]
    aucs, covariance = _fast_delong(stacked, n_positive)

    contrast = np.array([[-1.0, 1.0]])
    difference = float(aucs[1] - aucs[0])
    variance = float((contrast @ covariance @ contrast.T).item())

    if variance <= 0:
        # The two models produced effectively identical rankings.
        return difference, 1.0

    z_statistic = difference / np.sqrt(variance)
    return difference, float(2.0 * norm.sf(abs(z_statistic)))


def delong_difference_with_ci(
    y_true: np.ndarray,
    predictions_a: np.ndarray,
    predictions_b: np.ndarray,
    alpha: float = 0.05,
) -> tuple[float, float, float, float]:
    """Difference in area with a confidence interval and a two sided p value.

    The interval and the p value come from the same DeLong covariance, so they
    agree exactly: the interval excludes zero if and only if the p value falls
    below ``alpha``. The plan states the question 1 positive branch in terms of
    the interval, so the interval is what gets reported, with the p value
    alongside it rather than in place of it.
    """

    order, n_positive = _order_by_outcome(y_true)
    stacked = np.vstack(
        [np.asarray(predictions_a, dtype=float), np.asarray(predictions_b, dtype=float)]
    )[:, order]
    aucs, covariance = _fast_delong(stacked, n_positive)

    contrast = np.array([[-1.0, 1.0]])
    difference = float(aucs[1] - aucs[0])
    variance = float((contrast @ covariance @ contrast.T).item())

    if variance <= 0:
        return difference, difference, difference, 1.0

    standard_error = float(np.sqrt(variance))
    half_width = norm.ppf(1.0 - alpha / 2.0) * standard_error
    p_value = float(2.0 * norm.sf(abs(difference / standard_error)))
    return difference, difference - half_width, difference + half_width, p_value


# ---------------------------------------------------------------------------
# Calibration
# ---------------------------------------------------------------------------


def calibration_slope(y_true: np.ndarray, probabilities: np.ndarray) -> float:
    """Slope of the observed outcome on the predicted log odds.

    A slope below one indicates predictions that are too extreme, which is
    the usual signature of overfitting. A slope of one does not by itself
    mean the model is well calibrated, since the intercept may still be
    wrong.
    """

    linear_predictor = logit(_clip_probabilities(probabilities)).reshape(-1, 1)
    model = LogisticRegression(penalty=None, solver="lbfgs", max_iter=1000)
    model.fit(linear_predictor, np.asarray(y_true, dtype=int))
    return float(model.coef_[0][0])


def calibration_intercept(y_true: np.ndarray, probabilities: np.ndarray) -> float:
    """Calibration in the large, holding the slope fixed at one.

    Fitted as an intercept only logistic model with the predicted log odds
    as an offset. The score equation is monotonic in the intercept, so a
    bracketed root finder recovers the maximum likelihood estimate directly
    and avoids depending on a general purpose optimiser for a one parameter
    problem.
    """

    y_true = np.asarray(y_true, dtype=float)
    offset = logit(_clip_probabilities(probabilities))

    def score(intercept: float) -> float:
        return float(np.sum(y_true - expit(offset + intercept)))

    return float(brentq(score, -25.0, 25.0, xtol=1e-10))


@dataclass(frozen=True)
class CalibrationSummary:
    """Calibration measures for one set of predictions."""

    intercept: float
    slope: float
    brier: float
    brier_skill: float


def summarise_calibration(
    y_true: np.ndarray, probabilities: np.ndarray
) -> CalibrationSummary:
    """Collect the calibration measures reported for every model.

    The Brier skill score expresses the Brier score relative to the score of
    a model that predicts the observed event rate for everyone. Without that
    reference the Brier score is hard to read, because its scale depends on
    the event rate and low rates make every model look accurate.
    """

    y_true = np.asarray(y_true, dtype=int)
    probabilities = np.asarray(probabilities, dtype=float)

    brier = float(brier_score_loss(y_true, probabilities))
    reference = float(np.mean((y_true - y_true.mean()) ** 2))
    skill = 1.0 - brier / reference if reference > 0 else float("nan")

    return CalibrationSummary(
        intercept=calibration_intercept(y_true, probabilities),
        slope=calibration_slope(y_true, probabilities),
        brier=brier,
        brier_skill=float(skill),
    )


# ---------------------------------------------------------------------------
# Reclassification
# ---------------------------------------------------------------------------


def integrated_discrimination_improvement(
    y_true: np.ndarray, baseline: np.ndarray, extended: np.ndarray
) -> float:
    """Integrated discrimination improvement (Pencina et al., 2008).

    The difference in mean predicted risk between models, computed separately
    among cases and controls and then combined. It rewards a model that
    pushes cases up and controls down.
    """

    y_true = np.asarray(y_true, dtype=int)
    baseline = np.asarray(baseline, dtype=float)
    extended = np.asarray(extended, dtype=float)

    events = y_true == 1
    non_events = ~events
    gain_events = extended[events].mean() - baseline[events].mean()
    gain_non_events = extended[non_events].mean() - baseline[non_events].mean()
    return float(gain_events - gain_non_events)


def category_free_nri(
    y_true: np.ndarray, baseline: np.ndarray, extended: np.ndarray
) -> tuple[float, float, float]:
    """Category free net reclassification improvement (Pencina et al., 2011).

    Returns the event component, the non event component and their sum. The
    components are reported separately on purpose. A total near zero can hide
    a marker that helps materially among cases and harms among controls, and
    the two are not interchangeable when the clinical costs of a false
    negative and a false positive differ.
    """

    y_true = np.asarray(y_true, dtype=int)
    difference = np.asarray(extended, dtype=float) - np.asarray(baseline, dtype=float)

    events = difference[y_true == 1]
    non_events = difference[y_true == 0]

    nri_events = float(np.mean(events > 0) - np.mean(events < 0))
    nri_non_events = float(np.mean(non_events < 0) - np.mean(non_events > 0))
    return nri_events, nri_non_events, nri_events + nri_non_events


def bootstrap_interval(
    statistic,
    y_true: np.ndarray,
    baseline: np.ndarray,
    extended: np.ndarray,
    n_bootstrap: int,
    seed: int,
    alpha: float = 0.05,
) -> tuple[float, float]:
    """Percentile bootstrap interval for a reclassification statistic.

    Resampling is stratified by outcome so that every replicate keeps the
    observed number of events. With an event rate near ten per cent an
    unstratified bootstrap occasionally draws a replicate with too few events
    for the statistic to be defined.
    """

    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true, dtype=int)
    event_index = np.flatnonzero(y_true == 1)
    non_event_index = np.flatnonzero(y_true == 0)

    replicates = np.empty(n_bootstrap, dtype=float)
    for replicate in range(n_bootstrap):
        drawn = np.concatenate(
            [
                rng.choice(event_index, size=event_index.size, replace=True),
                rng.choice(non_event_index, size=non_event_index.size, replace=True),
            ]
        )
        value = statistic(y_true[drawn], baseline[drawn], extended[drawn])
        replicates[replicate] = value[-1] if isinstance(value, tuple) else value

    lower, upper = np.percentile(replicates, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(lower), float(upper)


# ---------------------------------------------------------------------------
# Clinical utility
# ---------------------------------------------------------------------------


def net_benefit(
    y_true: np.ndarray, probabilities: np.ndarray, thresholds: np.ndarray
) -> np.ndarray:
    """Net benefit across decision thresholds.

    At threshold ``t`` the benefit of a true positive is weighted against the
    harm of a false positive by the odds ``t / (1 - t)``. The threshold is
    therefore an explicit statement of how many unnecessary interventions a
    clinician would accept to prevent one event, which is the assumption that
    usually stays implicit when models are compared by area under the curve
    alone.
    """

    y_true = np.asarray(y_true, dtype=int)
    probabilities = np.asarray(probabilities, dtype=float)
    n = y_true.size

    benefits = np.empty(thresholds.size, dtype=float)
    for position, threshold in enumerate(thresholds):
        treated = probabilities >= threshold
        true_positives = int(np.sum(treated & (y_true == 1)))
        false_positives = int(np.sum(treated & (y_true == 0)))
        weight = threshold / (1.0 - threshold)
        benefits[position] = true_positives / n - (false_positives / n) * weight

    return benefits


def treat_all_net_benefit(y_true: np.ndarray, thresholds: np.ndarray) -> np.ndarray:
    """Net benefit of the strategy that treats everyone.

    One of the two reference strategies any model must beat to be worth
    using. The other, treating nobody, has a net benefit of zero everywhere.
    """

    prevalence = float(np.mean(np.asarray(y_true, dtype=int)))
    weight = thresholds / (1.0 - thresholds)
    return prevalence - (1.0 - prevalence) * weight
