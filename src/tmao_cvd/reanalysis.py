"""The three pre-specified questions of docs/reanalysis_plan.md.

Each question is implemented as a function returning a
:class:`~tmao_cvd.reporting.QuestionResult` whose verdict is decided by the
thresholds fixed in section 6 of the plan. The thresholds appear here as module
constants so that the code and the plan can be read against each other, and so
that changing one is a visible diff rather than an inline edit buried in a
conditional.

Nothing in this module chooses a threshold, a model or a comparison on the
basis of a result. All of that was fixed before any estimate existed.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import kstest, spearmanr
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import roc_curve
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .config import AnalysisConfig
from .evaluate import (
    auc_with_ci,
    bootstrap_interval,
    category_free_nri,
    delong_test,
    integrated_discrimination_improvement,
    net_benefit,
    summarise_calibration,
)
from .models import cross_validated_predictions, make_logistic_model
from .reporting import (
    PanelContext,
    QuestionResult,
    incremental_value_statement,
    order_structure_statement,
)
from .st001420 import (
    EXTENDED_FEATURES,
    OUTCOME,
    PRECURSORS,
    PRECURSOR_FEATURES,
    SAMPLE_INDEX,
    TMAO,
    metabolite_columns,
)

# --- Thresholds fixed in section 6 of the plan -----------------------------

#: Question 1 net benefit range, bracketing the 28 per cent event rate.
Q1_THRESHOLDS = np.linspace(0.05, 0.60, 56)

#: Question 2, reproduction of the published performance.
Q2_REPRODUCED_FLOOR = 0.85
Q2_CALIBRATION_SLOPE_RANGE = (0.80, 1.25)
#: Question 2, collapse.
Q2_COLLAPSED_FLOOR = 0.75
Q2_COLLAPSED_SLOPE_FLOOR = 0.70

#: Question 3, order structure detection.
Q3_SIGNIFICANT_PROPORTION_LIMIT = 0.10
Q3_KS_P_LIMIT = 0.001
Q3_MINIMUM_BLOCK_OBSERVATIONS = 30



def panel_context(frame: pd.DataFrame) -> PanelContext:
    """Locate the marker under test among every metabolite in the deposit.

    Required by :func:`~tmao_cvd.reporting.incremental_value_statement`, so
    question 1 cannot report its difference in area without it. This is a
    descriptive placement rather than an additional analysis: it re-expresses
    measurements already in the deposit on a common scale.
    """

    from sklearn.metrics import roc_auc_score

    outcome = frame[OUTCOME].to_numpy(int)
    complete = [m for m in metabolite_columns(frame) if frame[m].notna().all()]
    auc = pd.Series(
        {
            m: max(roc_auc_score(outcome, frame[m]), 1 - roc_auc_score(outcome, frame[m]))
            for m in complete
        }
    )
    precursor_auc = auc[list(PRECURSORS)]
    best = precursor_auc.idxmax()

    return PanelContext(
        marker=TMAO,
        univariate_auc=float(auc[TMAO]),
        percentile=100 * float((auc < auc[TMAO]).mean()),
        ranking_higher=int((auc > auc[TMAO]).sum()),
        panel_size=len(complete),
        best_precursor=str(best),
        best_precursor_auc=float(precursor_auc.max()),
        panel_median_auc=float(auc.median()),
    )


def question_one(
    frame: pd.DataFrame, config: AnalysisConfig
) -> tuple[QuestionResult, dict[str, np.ndarray]]:
    """Does TMAO add beyond its own dietary precursors?

    Baseline and extended models differ by log TMAO and nothing else, which is
    what makes the difference readable as the contribution of that term.
    """

    outcome = frame[OUTCOME].to_numpy(int)

    baseline = cross_validated_predictions(
        make_logistic_model, frame[PRECURSOR_FEATURES], outcome, config
    )
    extended = cross_validated_predictions(
        make_logistic_model, frame[EXTENDED_FEATURES], outcome, config
    )

    base_auc, base_lo, base_hi = auc_with_ci(outcome, baseline.averaged)
    ext_auc, ext_lo, ext_hi = auc_with_ci(outcome, extended.averaged)
    difference, p_value = delong_test(outcome, baseline.averaged, extended.averaged)

    base_cal = summarise_calibration(outcome, baseline.averaged)
    ext_cal = summarise_calibration(outcome, extended.averaged)

    idi = integrated_discrimination_improvement(outcome, baseline.averaged, extended.averaged)
    idi_lo, idi_hi = bootstrap_interval(
        integrated_discrimination_improvement, outcome, baseline.averaged,
        extended.averaged, config.n_bootstrap, config.random_seed,
    )
    nri_events, nri_non, nri_total = category_free_nri(
        outcome, baseline.averaged, extended.averaged
    )
    nri_lo, nri_hi = bootstrap_interval(
        category_free_nri, outcome, baseline.averaged, extended.averaged,
        config.n_bootstrap, config.random_seed + 1,
    )

    nb_base = net_benefit(outcome, baseline.averaged, Q1_THRESHOLDS)
    nb_ext = net_benefit(outcome, extended.averaged, Q1_THRESHOLDS)
    # "lies above across a substantial part of the range" is read as a
    # majority of the pre-specified threshold grid, by a margin large enough
    # not to be numerical noise.
    separation = float(np.mean(nb_ext > nb_base + 1e-4))

    # The DeLong interval on the difference is the pre-specified primary test.
    interval_excludes_zero = p_value < 0.05
    curves_separate = separation > 0.5

    if interval_excludes_zero and curves_separate:
        verdict = "POSITIVE"
    elif not interval_excludes_zero:
        verdict = "NEGATIVE"
    else:
        verdict = "DISCORDANT"

    # The difference in area is never formatted into text here. It is passed to
    # the sole constructor in tmao_cvd.reporting, which attaches the panel
    # context unconditionally. See the guard test in tests/test_reporting.py.
    statement = incremental_value_statement(
        verdict, difference, p_value, separation, panel_context(frame)
    )

    table = pd.DataFrame(
        [
            {
                "model": "baseline (precursors only)", "auc": round(base_auc, 4),
                "auc_ci_lower": round(base_lo, 4), "auc_ci_upper": round(base_hi, 4),
                "calibration_slope": round(base_cal.slope, 4),
                "brier_skill": round(base_cal.brier_skill, 4),
            },
            {
                "model": "extended (plus log TMAO)", "auc": round(ext_auc, 4),
                "auc_ci_lower": round(ext_lo, 4), "auc_ci_upper": round(ext_hi, 4),
                "calibration_slope": round(ext_cal.slope, 4),
                "brier_skill": round(ext_cal.brier_skill, 4),
            },
            {
                "model": "difference", "auc": round(difference, 4),
                "auc_ci_lower": np.nan, "auc_ci_upper": np.nan,
                "calibration_slope": np.nan, "brier_skill": np.nan,
            },
        ]
    )
    extra = pd.DataFrame(
        [{
            "delong_p": round(p_value, 5),
            "idi": round(idi, 4), "idi_ci_lower": round(idi_lo, 4),
            "idi_ci_upper": round(idi_hi, 4),
            "nri_events": round(nri_events, 4), "nri_non_events": round(nri_non, 4),
            "nri_total": round(nri_total, 4), "nri_ci_lower": round(nri_lo, 4),
            "nri_ci_upper": round(nri_hi, 4),
            "decision_curve_separation": round(separation, 4),
        }]
    )

    result = QuestionResult(
        question="QUESTION 1 (primary). Does TMAO add beyond its own precursors?",
        verdict=verdict,
        statement=statement,
        table=pd.concat([table, extra], axis=0, ignore_index=True).fillna(""),
    )
    curves = {
        "thresholds": Q1_THRESHOLDS,
        "net_benefit_baseline": nb_base,
        "net_benefit_extended": nb_ext,
        "baseline_predictions": baseline.averaged,
        "extended_predictions": extended.averaged,
    }
    return result, curves


def _panel_model(seed: int) -> Pipeline:
    """L2 penalised logistic regression with the penalty chosen in-fold.

    ``LogisticRegressionCV`` runs its own cross validation over the training
    data it is given. Placed inside the outer resampling loop, that makes the
    penalty selection nested, so no information from a held out fold reaches
    model selection. Imputation and scaling sit inside the same pipeline for
    the same reason.
    """

    return Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            (
                "model",
                LogisticRegressionCV(
                    Cs=np.logspace(-3, 1, 5), cv=5, penalty="l2", solver="lbfgs",
                    scoring="roc_auc", max_iter=2000, n_jobs=-1, random_state=seed,
                ),
            ),
        ]
    )


def question_two(frame: pd.DataFrame, config: AnalysisConfig) -> QuestionResult:
    """What does the deposited panel support under honest out of fold validation?

    A figure of above 89 per cent on accuracy, sensitivity and specificity
    simultaneously has previously been reported for this deposit. What follows
    tests what the deposited artefact supports, and is not a comment on anyone
    who produced or deposited it. Note one asymmetry, stated rather than hidden:
    the operating threshold is chosen by maximising the Youden index on the same
    out of fold predictions that are then scored. That is mildly optimistic and
    biases the comparison towards reproducing the published figures, so it makes
    a collapse verdict harder to reach, not easier.
    """

    outcome = frame[OUTCOME].to_numpy(int)
    panel = frame[metabolite_columns(frame)]

    predictions = cross_validated_predictions(_panel_model, panel, outcome, config)
    probabilities = predictions.averaged

    auc, auc_lo, auc_hi = auc_with_ci(outcome, probabilities)
    calibration = summarise_calibration(outcome, probabilities)

    false_positive, true_positive, cuts = roc_curve(outcome, probabilities)
    youden = int(np.argmax(true_positive - false_positive))
    threshold = float(cuts[youden])

    predicted = probabilities >= threshold
    sensitivity = float(np.mean(predicted[outcome == 1]))
    specificity = float(np.mean(~predicted[outcome == 0]))
    accuracy = float(np.mean(predicted == outcome.astype(bool)))

    trio = (accuracy, sensitivity, specificity)
    slope = calibration.slope

    reproduced = (
        min(trio) >= Q2_REPRODUCED_FLOOR
        and Q2_CALIBRATION_SLOPE_RANGE[0] <= slope <= Q2_CALIBRATION_SLOPE_RANGE[1]
    )
    collapsed = min(trio) < Q2_COLLAPSED_FLOOR or slope < Q2_COLLAPSED_SLOPE_FLOOR

    if reproduced:
        verdict = "REPRODUCED"
        statement = (
            "The reported discrimination is robust to honest out of fold resampling. "
            f"Accuracy {accuracy:.3f}, sensitivity {sensitivity:.3f}, specificity "
            f"{specificity:.3f}, calibration slope {slope:.3f}. Given section 5 of the "
            "plan, this still does not establish a biological basis for the signal, "
            "because acquisition order is confounded with outcome."
        )
    elif collapsed:
        verdict = "COLLAPSED"
        statement = (
            "The published performance does not survive out of fold evaluation on the "
            f"deposited cohort. Accuracy {accuracy:.3f}, sensitivity {sensitivity:.3f}, "
            f"specificity {specificity:.3f}, calibration slope {slope:.3f}, against a "
            "previously reported figure of above 0.890 on all three. This is the "
            "finding. It is reported exactly as prominently as a confirmation would "
            "have been, and it concerns what the deposited artefact supports rather "
            "than anyone who produced or deposited it."
        )
    else:
        verdict = "PARTIAL"
        # The plan's partial branch was written imagining discrimination falling
        # between the floors. The other way to land here is discrimination above
        # the reproduction floor with calibration outside the acceptable band,
        # which is what occurred. The sentence distinguishes the two rather than
        # asserting the one the plan happened to anticipate. Logged as a post hoc
        # correction in section 8 of the plan; the verdict itself is unchanged.
        if min(trio) >= Q2_REPRODUCED_FLOOR:
            reason = (
                "Discrimination is at or above the reproduction floor, but the "
                f"calibration slope of {slope:.3f} lies outside the pre-specified band "
                f"{Q2_CALIBRATION_SLOPE_RANGE}, so the reproduction criterion is not met. "
                "A slope this far from one means the predicted risks are not usable as "
                "probabilities even though the ranking is near perfect."
            )
        else:
            reason = (
                "These fall between the pre-specified reproduction and collapse "
                "thresholds."
            )
        statement = (
            f"Accuracy {accuracy:.3f}, sensitivity {sensitivity:.3f}, specificity "
            f"{specificity:.3f}, calibration slope {slope:.3f}. {reason} Per the plan "
            "the numbers are reported with no verdict attached, and the gap from the "
            "previously reported figure of above 0.890 on all three is stated without "
            "interpretation."
        )

    table = pd.DataFrame(
        [{
            "metric": "accuracy", "out_of_fold": round(accuracy, 4), "published": ">0.890",
        }, {
            "metric": "sensitivity", "out_of_fold": round(sensitivity, 4), "published": ">0.890",
        }, {
            "metric": "specificity", "out_of_fold": round(specificity, 4), "published": ">0.890",
        }, {
            "metric": "auc", "out_of_fold": round(auc, 4),
            "published": f"[{auc_lo:.3f}, {auc_hi:.3f}] CI",
        }, {
            "metric": "calibration_slope", "out_of_fold": round(slope, 4), "published": "not reported",
        }, {
            "metric": "calibration_intercept", "out_of_fold": round(calibration.intercept, 4),
            "published": "not reported",
        }, {
            "metric": "youden_threshold", "out_of_fold": round(threshold, 4), "published": "not reported",
        }]
    )

    return QuestionResult(
        question="QUESTION 2 (replication). What does the deposited panel support?",
        verdict=verdict,
        statement=statement,
        table=table,
    )


def question_three(frame: pd.DataFrame) -> QuestionResult:
    """Is acquisition order structure detectable along the deposited ordering?

    Computed within each outcome block separately so the test is blind to the
    between group difference. Under the null of no order structure the Spearman
    p values are uniform, so both the proportion below 0.05 and a
    Kolmogorov-Smirnov test against the uniform are informative.
    """

    metabolites = metabolite_columns(frame)
    rows: list[dict[str, object]] = []

    for label, mask in (("cases", frame[OUTCOME] == 1), ("controls", frame[OUTCOME] == 0)):
        block = frame.loc[mask]
        index = block[SAMPLE_INDEX].to_numpy(float)

        p_values: list[float] = []
        for name in metabolites:
            values = block[name].to_numpy(float)
            usable = ~np.isnan(values)
            if usable.sum() < Q3_MINIMUM_BLOCK_OBSERVATIONS:
                continue
            if np.unique(values[usable]).size < 3:
                continue
            _, p = spearmanr(index[usable], values[usable])
            if not np.isnan(p):
                p_values.append(float(p))

        p_array = np.array(p_values)
        proportion = float(np.mean(p_array < 0.05))
        ks_p = float(kstest(p_array, "uniform").pvalue)

        rows.append({
            "block": label,
            "n_participants": int(mask.sum()),
            "metabolites_tested": int(p_array.size),
            "proportion_p_below_0.05": round(proportion, 4),
            "null_expectation": 0.05,
            "ks_vs_uniform_p": f"{ks_p:.3g}",
            "exceeds_limit": bool(
                proportion > Q3_SIGNIFICANT_PROPORTION_LIMIT or ks_p < Q3_KS_P_LIMIT
            ),
        })

    table = pd.DataFrame(rows)
    detected = bool(table["exceeds_limit"].any())

    detail = "; ".join(
        f"{r['block']}: {r['proportion_p_below_0.05']:.1%} of "
        f"{r['metabolites_tested']} metabolites at p < 0.05 against a null of 5%, "
        f"KS p = {r['ks_vs_uniform_p']}"
        for r in rows
    )

    return QuestionResult(
        question="QUESTION 3 (diagnostic). Is acquisition order structure detectable?",
        verdict="ORDER STRUCTURE DETECTED" if detected else "NOT DETECTED",
        statement=order_structure_statement(detected, detail),
        table=table,
    )
