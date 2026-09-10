"""Synthetic cohort generation used to validate the analysis pipeline.

This repository does not distribute participant level data. To keep the
analysis code executable and testable in the absence of a cohort, this
module draws a synthetic population whose generative model is written out
in full below. Every parameter is an assumption chosen to give a plausible
looking dataset. None of them is an empirical estimate, and no number
computed from these data carries biological meaning. The purpose is
software validation, not inference.

Structure of the generative model
---------------------------------
Covariates are drawn in an order that respects the causal ordering usually
assumed in cardiovascular epidemiology. Age and sex are treated as
exogenous. Body mass index follows, then the metabolic and haemodynamic
variables that depend on it, then renal function, and finally the
biomarkers.

Plasma TMAO is generated as a log normal variable whose location depends on
age, renal function and diabetes status. The dependence on estimated
glomerular filtration rate is deliberate rather than decorative. TMAO is
cleared renally, and part of the observed association between TMAO and
cardiovascular outcomes is attributable to renal impairment (Tang et al.,
Circulation Research, 2015). A simulation that ignored this would make the
incremental value question artificially easy and would flatter the code.

The outcome is a binary indicator of a major adverse cardiovascular event
within three years of sampling, drawn from a logistic model on standardised
covariates. The intercept is solved numerically so that the marginal event
rate matches the requested value, which keeps the event rate stable when
the coefficients are edited.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.special import expit

from .config import ID_COLUMN, OUTCOME


@dataclass(frozen=True)
class SimulationParameters:
    """Assumed parameters of the synthetic data generating process.

    The log odds are expressed per standard deviation of the corresponding
    covariate so that their relative magnitudes can be compared directly.
    They are loosely informed by the size of effect usually reported for
    these variables, but they are assumptions and should be read as such.
    """

    n_participants: int = 4000
    seed: int = 20260910
    target_event_rate: float = 0.10

    #: Log odds of the outcome per standard deviation of each covariate.
    outcome_log_odds: dict[str, float] = field(
        default_factory=lambda: {
            "age_years": 0.45,
            "sex_male": 0.30,
            "current_smoker": 0.35,
            "diabetes": 0.40,
            "bmi_kg_m2": 0.10,
            "systolic_bp_mmhg": 0.25,
            "total_cholesterol_mmol_l": 0.20,
            "hdl_cholesterol_mmol_l": -0.22,
            "egfr_ml_min_1_73m2": -0.28,
            "log_hs_crp": 0.18,
            "log_tmao": 0.22,
        }
    )

    #: Residual standard deviation of log transformed TMAO.
    tmao_log_sd: float = 0.55
    #: Median TMAO in micromoles per litre for a reference participant.
    tmao_reference_median: float = 3.5


def _solve_intercept(linear_predictor: np.ndarray, target_rate: float) -> float:
    """Find the intercept that gives the requested marginal event rate.

    Fixing the intercept by hand couples the event rate to the covariate
    coefficients, so any change to the latter silently shifts the former.
    Solving for it keeps the two independent. The objective is monotonic in
    the intercept, so a bracketed root finder is sufficient and exact to
    machine precision.
    """

    def objective(intercept: float) -> float:
        return float(expit(linear_predictor + intercept).mean() - target_rate)

    return float(brentq(objective, -25.0, 25.0, xtol=1e-10))


def _standardise(values: np.ndarray) -> np.ndarray:
    """Centre and scale, guarding against a degenerate constant column."""
    spread = values.std()
    if spread == 0:
        return np.zeros_like(values, dtype=float)
    return (values - values.mean()) / spread


def simulate_cohort(parameters: SimulationParameters | None = None) -> pd.DataFrame:
    """Draw a synthetic cohort.

    Parameters
    ----------
    parameters:
        Generative assumptions. The default set is documented on
        :class:`SimulationParameters`.

    Returns
    -------
    pandas.DataFrame
        One row per participant, using the column names defined in the data
        dictionary so that simulated and real data are interchangeable from
        the point of view of the rest of the pipeline.
    """

    parameters = parameters or SimulationParameters()
    rng = np.random.default_rng(parameters.seed)
    n = parameters.n_participants

    age = np.clip(rng.normal(62.0, 10.0, n), 40.0, 85.0)
    sex_male = rng.binomial(1, 0.55, n)

    bmi = np.clip(rng.normal(28.0, 4.5, n), 17.0, 50.0)

    # Smoking is made mildly more common in younger participants, which is
    # the direction seen in most contemporary cohorts.
    smoking_p = expit(0.4 - 0.035 * (age - 62.0))
    current_smoker = rng.binomial(1, smoking_p)

    diabetes_p = expit(-2.3 + 0.09 * (bmi - 28.0) + 0.02 * (age - 62.0))
    diabetes = rng.binomial(1, diabetes_p)

    systolic_bp = np.clip(
        rng.normal(132.0 + 0.35 * (age - 62.0) + 0.50 * (bmi - 28.0), 16.0, n),
        85.0,
        220.0,
    )
    total_cholesterol = np.clip(rng.normal(5.1, 1.0, n), 2.0, 11.0)
    hdl_cholesterol = np.clip(
        rng.normal(1.35 - 0.15 * sex_male - 0.02 * (bmi - 28.0), 0.33, n),
        0.4,
        3.5,
    )

    # Renal function declines with age and is lower in diabetes.
    egfr = np.clip(
        rng.normal(88.0 - 0.70 * (age - 62.0) - 6.0 * diabetes, 15.0, n),
        15.0,
        130.0,
    )

    hs_crp = np.clip(
        np.exp(rng.normal(0.35 + 0.045 * (bmi - 28.0) + 0.25 * diabetes, 0.85, n)),
        0.1,
        60.0,
    )

    # TMAO rises as renal clearance falls, hence the negative coefficient on
    # eGFR. The intercept is set so that a reference participant sits at the
    # assumed median.
    log_tmao_mean = (
        np.log(parameters.tmao_reference_median)
        + 0.08 * (age - 62.0) / 10.0
        - 0.30 * (egfr - 88.0) / 15.0
        + 0.18 * diabetes
    )
    tmao = np.exp(rng.normal(log_tmao_mean, parameters.tmao_log_sd, n))

    frame = pd.DataFrame(
        {
            ID_COLUMN: [f"SIM{index:05d}" for index in range(1, n + 1)],
            "age_years": age,
            "sex_male": sex_male,
            "current_smoker": current_smoker,
            "diabetes": diabetes,
            "bmi_kg_m2": bmi,
            "systolic_bp_mmhg": systolic_bp,
            "total_cholesterol_mmol_l": total_cholesterol,
            "hdl_cholesterol_mmol_l": hdl_cholesterol,
            "egfr_ml_min_1_73m2": egfr,
            "hs_crp_mg_l": hs_crp,
            "tmao_umol_l": tmao,
        }
    )

    # The outcome model is written on the same transformed scale that the
    # analysis uses, so the simulated effect of TMAO is directly comparable
    # with the coefficient the models are asked to recover.
    design = frame.copy()
    design["log_hs_crp"] = np.log(design["hs_crp_mg_l"])
    design["log_tmao"] = np.log(design["tmao_umol_l"])

    linear_predictor = np.zeros(n, dtype=float)
    for column, log_odds in parameters.outcome_log_odds.items():
        linear_predictor += log_odds * _standardise(design[column].to_numpy(float))

    intercept = _solve_intercept(linear_predictor, parameters.target_event_rate)
    event_probability = expit(linear_predictor + intercept)
    frame[OUTCOME] = rng.binomial(1, event_probability)

    return frame
