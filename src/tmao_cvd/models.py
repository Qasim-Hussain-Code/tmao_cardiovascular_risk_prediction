"""Model specifications and the resampling scheme used to score them.

Penalised regression and tree ensembles are both defensible choices here.
Logistic regression is treated as the primary model because the question is
about the incremental value of one term, and a model whose coefficients can
be read directly makes that question easier to answer honestly. Gradient
boosting is fitted alongside it as a sensitivity analysis, to check that a
null result is not simply an artefact of assuming a linear effect of log
TMAO on the log odds.

All performance figures come from out of fold predictions. Predictions from
the data a model was fitted on are optimistic by an amount that grows with
the number of candidate terms, which is exactly the setting where a new
biomarker is being evaluated.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .config import AnalysisConfig


def make_logistic_model(seed: int) -> Pipeline:
    """Unpenalised logistic regression with imputation and standardisation.

    Standardisation is applied so that coefficients are expressed per
    standard deviation and are comparable across variables measured in
    different units. It is fitted inside the pipeline, and therefore inside
    each training fold, so that no information from the held out fold
    reaches the transformation.

    No penalty is applied. With eleven candidate terms and several hundred
    events the maximum likelihood fit is stable, and an unpenalised fit keeps
    the coefficients interpretable as log odds ratios.
    """

    return Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            (
                "model",
                LogisticRegression(
                    penalty=None, solver="lbfgs", max_iter=2000, random_state=seed
                ),
            ),
        ]
    )


def make_gradient_boosting_model(seed: int) -> Pipeline:
    """Gradient boosted trees, used as a sensitivity analysis.

    The hyperparameters are shrunk relative to the library defaults because
    the cohort is small by the standards of tree ensembles. A deliberately
    constrained learner is the right comparison here: the purpose is to
    detect a non linear signal that logistic regression would miss, not to
    win a prediction contest.
    """

    return Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            (
                "model",
                HistGradientBoostingClassifier(
                    learning_rate=0.05,
                    max_depth=3,
                    max_iter=250,
                    min_samples_leaf=40,
                    l2_regularization=1.0,
                    early_stopping=False,
                    random_state=seed,
                ),
            ),
        ]
    )


@dataclass(frozen=True)
class CrossValidatedPredictions:
    """Out of fold risk predictions from a repeated cross validation.

    Attributes
    ----------
    averaged:
        Mean out of fold prediction per participant across repeats. Used for
        the reported statistics, since averaging removes some of the noise
        introduced by one arbitrary partition of the data.
    per_repeat:
        Array of shape ``(n_repeats, n_participants)``. Retained so that the
        spread of performance across partitions can be reported, which is the
        honest way to convey how much a single split would have told us.
    """

    averaged: np.ndarray
    per_repeat: np.ndarray


def cross_validated_predictions(
    build_model,
    features: pd.DataFrame,
    outcome: np.ndarray,
    config: AnalysisConfig,
) -> CrossValidatedPredictions:
    """Produce out of fold predictions under a repeated stratified scheme.

    The splitter is seeded from the configuration, so every model in a run
    sees identical folds. That is a requirement rather than a convenience:
    the comparison between the baseline and extended models is paired, and
    pairing is only valid if both models were scored on the same partitions.
    """

    x = features.to_numpy(dtype=float)
    y = np.asarray(outcome, dtype=int)

    splitter = RepeatedStratifiedKFold(
        n_splits=config.n_splits,
        n_repeats=config.n_repeats,
        random_state=config.random_seed,
    )

    per_repeat = np.full((config.n_repeats, y.size), np.nan, dtype=float)

    for fold_index, (train_index, test_index) in enumerate(splitter.split(x, y)):
        repeat = fold_index // config.n_splits
        model = build_model(config.random_seed + fold_index)
        model.fit(x[train_index], y[train_index])
        per_repeat[repeat, test_index] = model.predict_proba(x[test_index])[:, 1]

    if np.isnan(per_repeat).any():
        raise RuntimeError("some participants were never held out during resampling")

    return CrossValidatedPredictions(
        averaged=per_repeat.mean(axis=0), per_repeat=per_repeat
    )


def logistic_coefficient_table(
    model: Pipeline, features: pd.DataFrame
) -> pd.DataFrame:
    """Odds ratios per standard deviation, with Wald confidence intervals.

    scikit-learn does not expose standard errors, so the observed information
    is reconstructed as the inverse of ``Z' W Z``, where ``Z`` is the design
    matrix after the pipeline transformations and ``W`` holds the fitted
    variances ``p (1 - p)``. This is the standard asymptotic covariance of an
    unpenalised logistic fit, and it is valid here only because no penalty is
    applied. Reading Wald intervals off a penalised fit would understate the
    uncertainty.
    """

    design = model[:-1].transform(features.to_numpy(dtype=float))
    classifier = model[-1]

    fitted = classifier.predict_proba(design)[:, 1]
    weights = fitted * (1.0 - fitted)

    with_intercept = np.hstack([np.ones((design.shape[0], 1)), design])
    information = (with_intercept.T * weights) @ with_intercept
    covariance = np.linalg.pinv(information)
    standard_errors = np.sqrt(np.diag(covariance))[1:]

    coefficients = classifier.coef_[0]
    lower = coefficients - 1.96 * standard_errors
    upper = coefficients + 1.96 * standard_errors

    return pd.DataFrame(
        {
            "term": list(features.columns),
            "log_odds_per_sd": coefficients,
            "standard_error": standard_errors,
            "odds_ratio_per_sd": np.exp(coefficients),
            "ci_lower": np.exp(lower),
            "ci_upper": np.exp(upper),
        }
    )


#: Models scored in a full run, in the order they are reported.
MODEL_BUILDERS = {
    "logistic regression": make_logistic_model,
    "gradient boosting": make_gradient_boosting_model,
}
