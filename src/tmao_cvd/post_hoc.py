"""Diagnostics run AFTER the pre-specified questions returned their verdicts.

Nothing in this module was pre-specified. It is separated from
:mod:`tmao_cvd.reanalysis` so that the distinction survives contact with a
reader who is skimming, and every function here is labelled post hoc in the
output it produces.

These diagnostics exist because question 2 returned an out of fold area under
the curve of exactly 1.000 on 750 participants. That is not a performance
figure to be reported and interpreted, it is a signal that something
structural separates the two groups, and leaving it uninvestigated would have
been the more serious failure.

The order in which they were run is preserved below, because each was chosen
in response to what the previous one ruled out.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from .st001420 import OUTCOME, PRECURSORS, TMAO, metabolite_columns

#: Proportion of a panel expected to reach a univariate AUC of 0.6 by chance.
#: Used only as a stated reference point, not as a test.
NULL_PANEL_EXPECTATION = "roughly 5 to 10 per cent"


def missingness_by_outcome(frame: pd.DataFrame) -> pd.DataFrame:
    """Does the pattern of missing values encode the outcome?

    The first hypothesis for a perfect separation, because a metabolite absent
    in one group and present in the other is a perfect predictor that survives
    imputation. Run first so that it could be ruled out before anything more
    elaborate was attempted.
    """

    metabolites = metabolite_columns(frame)
    missing = frame[metabolites].isna()
    outcome = frame[OUTCOME]

    gap = (missing[outcome == 1].mean() - missing[outcome == 0].mean()).abs()
    return pd.DataFrame(
        [{
            "diagnostic": "missingness gap between groups (post hoc)",
            "metabolites": len(metabolites),
            "fully_observed": int((missing.mean() == 0).sum()),
            "absent_for_everyone": int((missing.mean() == 1).sum()),
            "max_group_gap": round(float(gap.max()), 6),
            "metabolites_with_gap_above_0.2": int((gap > 0.2).sum()),
        }]
    )


def univariate_panel_discrimination(frame: pd.DataFrame) -> pd.DataFrame:
    """How much of the panel discriminates the two groups on its own?

    The shape of this distribution is the diagnostic. Real biology produces a
    handful of discriminating metabolites against a background centred on 0.5.
    A panel where most metabolites discriminate is the signature of a
    systematic difference between the two sets of samples rather than of a
    disease process with hundreds of independent markers.
    """

    metabolites = [
        m for m in metabolite_columns(frame) if frame[m].notna().all()
    ]
    outcome = frame[OUTCOME].to_numpy(int)

    # Direction free, since a batch shift can move a metabolite either way.
    auc = pd.Series(
        {
            m: max(
                roc_auc_score(outcome, frame[m]), 1 - roc_auc_score(outcome, frame[m])
            )
            for m in metabolites
        }
    )

    rows = [{
        "diagnostic": "panel-wide univariate discrimination (post hoc)",
        "metabolites_tested": len(metabolites),
        "median_auc": round(float(auc.median()), 4),
        "median_auc_under_null": 0.5,
        "proportion_above_0.6": round(float((auc > 0.6).mean()), 4),
        "null_expectation": NULL_PANEL_EXPECTATION,
        "proportion_above_0.8": round(float((auc > 0.8).mean()), 4),
        "proportion_above_0.9": round(float((auc > 0.9).mean()), 4),
    }]
    return pd.DataFrame(rows), auc


def pathway_position_in_panel(auc: pd.Series) -> pd.DataFrame:
    """Where TMAO and its precursors rank among all measured metabolites.

    This is the diagnostic that bears hardest on question 1. If TMAO is an
    unremarkable member of a globally shifted panel, then its incremental
    contribution over three other members of that panel is not evidence that it
    carries biological information about the outcome.
    """

    rows = []
    for name in [TMAO] + list(PRECURSORS):
        rows.append({
            "metabolite": name,
            "univariate_auc": round(float(auc[name]), 4),
            "percentile_of_panel": round(100 * float((auc < auc[name]).mean()), 1),
            "metabolites_ranking_higher": int((auc > auc[name]).sum()),
        })
    return pd.DataFrame(rows)


def total_intensity_discrimination(frame: pd.DataFrame) -> pd.DataFrame:
    """Does the crude sum of all signal separate the groups?

    A panel-wide intensity difference cannot plausibly be a disease signature.
    It is what differing sample handling, extraction or instrument state
    between two sets of samples looks like.
    """

    metabolites = [m for m in metabolite_columns(frame) if frame[m].notna().all()]
    outcome = frame[OUTCOME].to_numpy(int)
    total = np.log(frame[metabolites]).sum(axis=1)
    auc = roc_auc_score(outcome, total)

    return pd.DataFrame(
        [{
            "diagnostic": "total log intensity across the panel (post hoc)",
            "metabolites_summed": len(metabolites),
            "auc_direction_free": round(float(max(auc, 1 - auc)), 4),
            "auc_under_null": 0.5,
        }]
    )


def run_all(frame: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Run every post hoc diagnostic, in the order they were originally run."""

    panel, auc = univariate_panel_discrimination(frame)
    return {
        "missingness": missingness_by_outcome(frame),
        "panel_discrimination": panel,
        "pathway_position": pathway_position_in_panel(auc),
        "total_intensity": total_intensity_discrimination(frame),
    }
