"""Feature definitions and the transformations applied before modelling.

Two design decisions are worth stating explicitly, since both affect how the
results should be read.

First, TMAO and high sensitivity C reactive protein are analysed on the
natural log scale. Both are strongly right skewed, and both are conventionally
modelled that way. On the untransformed scale a handful of very high values
dominate the fit of a linear term, which is a property of the measurement
distribution rather than of the biology.

Second, the baseline model is fixed in advance and contains the variables a
clinician would already have. The extended model is the baseline plus log
TMAO and nothing else. Keeping the two models identical apart from the single
term under test is what makes the comparison interpretable as the incremental
value of that term (Hlatky et al., Circulation, 2009).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

#: Established risk factors available before any biomarker is measured.
BASELINE_FEATURES: list[str] = [
    "age_years",
    "sex_male",
    "current_smoker",
    "diabetes",
    "bmi_kg_m2",
    "systolic_bp_mmhg",
    "total_cholesterol_mmol_l",
    "hdl_cholesterol_mmol_l",
    "egfr_ml_min_1_73m2",
    "log_hs_crp",
]

#: The single term whose incremental value the study is designed to estimate.
CANDIDATE_MARKER = "log_tmao"

#: Baseline model plus the candidate marker.
EXTENDED_FEATURES: list[str] = BASELINE_FEATURES + [CANDIDATE_MARKER]


def add_derived_features(frame: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of ``frame`` with the transformed variables added.

    The function does not modify its argument. Log transforms are applied
    to strictly positive quantities only, which the schema already
    guarantees, so no offset is needed and none is added. Adding a small
    constant to avoid a log of zero would change the shape of the lower tail
    for no benefit here.
    """

    out = frame.copy()
    out["log_tmao"] = np.log(out["tmao_umol_l"].astype(float))
    out["log_hs_crp"] = np.log(out["hs_crp_mg_l"].astype(float))
    return out


def describe_by_outcome(frame: pd.DataFrame, outcome: str) -> pd.DataFrame:
    """Summarise the cohort by outcome status.

    Continuous variables are reported as mean and standard deviation, binary
    variables as a count and percentage. The table is descriptive. No tests
    are reported against it, because a significance test on baseline
    characteristics answers a question nobody asked and invites reading a
    p value as though it bore on the study hypothesis.
    """

    rows: list[dict[str, object]] = []
    groups = {
        "no event": frame[frame[outcome] == 0],
        "event": frame[frame[outcome] == 1],
    }

    for column in frame.columns:
        if column == outcome or not pd.api.types.is_numeric_dtype(frame[column]):
            continue

        is_binary = set(frame[column].dropna().unique()).issubset({0, 1})
        row: dict[str, object] = {"variable": column}
        row["summary"] = "n (%)" if is_binary else "mean (SD)"

        for label, group in groups.items():
            values = group[column].dropna().astype(float)
            if is_binary:
                row[label] = f"{int(values.sum())} ({100 * values.mean():.1f})"
            else:
                row[label] = f"{values.mean():.2f} ({values.std():.2f})"
        rows.append(row)

    return pd.DataFrame(rows)
