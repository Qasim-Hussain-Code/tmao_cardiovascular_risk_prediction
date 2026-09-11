"""Single machine readable summary of every headline figure.

Every number quoted in the README must trace to this file, and the audit
script in scripts/audit_readme.py enforces that by reading the JSON and
checking each figure appears in the README at the same precision. Writing the
summary from the same tables that produce the printed report, rather than
transcribing from a terminal session, is what makes that check meaningful.

Each question carries its verdict exactly as it fired, in the vocabulary fixed
by docs/reanalysis_plan.md before any estimate existed.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from .reporting import SCOPE_STATEMENT, PanelContext, QuestionResult

SUMMARY_NAME = "00_summary.json"


def _row(table: pd.DataFrame, column: str, value: str) -> pd.Series:
    """Fetch exactly one row, failing loudly rather than guessing."""
    matched = table[table[column] == value]
    if len(matched) != 1:
        raise ValueError(f"expected one row where {column} == {value!r}, got {len(matched)}")
    return matched.iloc[0]


def build_summary(
    results: dict[str, QuestionResult],
    panel: PanelContext,
    diagnostics: dict[str, pd.DataFrame],
    metadata,
) -> dict:
    """Assemble every headline figure into one structure."""

    q1, q2, q3 = results["q1"], results["q2"], results["q3"]

    t1 = q1.table
    baseline = _row(t1, "model", "baseline (precursors only)")
    extended = _row(t1, "model", "extended (plus log TMAO)")
    difference = _row(t1, "model", "difference")
    stats = t1[t1["model"] == ""].iloc[0]

    t2 = q2.table
    def metric(name: str) -> float:
        return float(_row(t2, "metric", name)["out_of_fold"])

    t3 = q3.table
    cases = _row(t3, "block", "cases")
    controls = _row(t3, "block", "controls")

    panel_disc = diagnostics["panel_discrimination"].iloc[0]
    missing = diagnostics["missingness"].iloc[0]
    intensity = diagnostics["total_intensity"].iloc[0]

    return {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "scope": SCOPE_STATEMENT,
        "dataset": {
            "source": "a public Metabolomics Workbench deposit",
            "accession": metadata.study_id,
            "licence": metadata.licence,
            "n_participants": metadata.n_participants,
            "n_events": metadata.n_events,
            "n_metabolites": metadata.n_metabolites,
            "measurement_scale": metadata.scale.value,
        },
        "question_1": {
            "question": "Does TMAO add beyond its own dietary precursors?",
            "verdict": q1.verdict,
            "baseline_auc": float(baseline["auc"]),
            "baseline_auc_ci": [float(baseline["auc_ci_lower"]), float(baseline["auc_ci_upper"])],
            "extended_auc": float(extended["auc"]),
            "extended_auc_ci": [float(extended["auc_ci_lower"]), float(extended["auc_ci_upper"])],
            "delta_auc": float(difference["auc"]),
            "delta_auc_ci": [float(difference["auc_ci_lower"]), float(difference["auc_ci_upper"])],
            "delong_p": float(stats["delong_p"]),
            "idi": float(stats["idi"]),
            "idi_ci": [float(stats["idi_ci_lower"]), float(stats["idi_ci_upper"])],
            "nri_total": float(stats["nri_total"]),
            "nri_ci": [float(stats["nri_ci_lower"]), float(stats["nri_ci_upper"])],
            "decision_curve_separation": float(stats["decision_curve_separation"]),
            "panel_context": {
                "marker_univariate_auc": round(panel.univariate_auc, 4),
                "percentile_of_panel": round(panel.percentile, 1),
                "metabolites_ranking_higher": panel.ranking_higher,
                "panel_size": panel.panel_size,
                "best_precursor": panel.best_precursor,
                "best_precursor_auc": round(panel.best_precursor_auc, 4),
                "panel_median_auc": round(panel.panel_median_auc, 4),
            },
        },
        "question_2": {
            "question": "What does the deposited panel support under honest validation?",
            "verdict": q2.verdict,
            "accuracy": metric("accuracy"),
            "sensitivity": metric("sensitivity"),
            "specificity": metric("specificity"),
            "auc": metric("auc"),
            "calibration_slope": metric("calibration_slope"),
            "calibration_intercept": metric("calibration_intercept"),
            "previously_reported_figure": "above 0.890 on accuracy, sensitivity and specificity",
        },
        "question_3": {
            "question": "Is acquisition order structure detectable?",
            "verdict": q3.verdict,
            "null_expectation": 0.05,
            "cases": {
                "n_participants": int(cases["n_participants"]),
                "metabolites_tested": int(cases["metabolites_tested"]),
                "proportion_p_below_0_05": float(cases["proportion_p_below_0.05"]),
                "ks_vs_uniform_p": str(cases["ks_vs_uniform_p"]),
            },
            "controls": {
                "n_participants": int(controls["n_participants"]),
                "metabolites_tested": int(controls["metabolites_tested"]),
                "proportion_p_below_0_05": float(controls["proportion_p_below_0.05"]),
                "ks_vs_uniform_p": str(controls["ks_vs_uniform_p"]),
            },
        },
        "post_hoc": {
            "note": "not pre-specified; prompted by the question 2 result",
            "panel_median_auc": float(panel_disc["median_auc"]),
            "panel_median_auc_under_null": float(panel_disc["median_auc_under_null"]),
            "proportion_above_0_6": float(panel_disc["proportion_above_0.6"]),
            "proportion_above_0_8": float(panel_disc["proportion_above_0.8"]),
            "proportion_above_0_9": float(panel_disc["proportion_above_0.9"]),
            "metabolites_fully_observed": int(missing["fully_observed"]),
            "metabolites_absent_for_everyone": int(missing["absent_for_everyone"]),
            "missingness_max_group_gap": float(missing["max_group_gap"]),
            "total_intensity_auc": float(intensity["auc_direction_free"]),
        },
    }


def write_summary(summary: dict, directory: Path) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / SUMMARY_NAME
    path.write_text(json.dumps(summary, indent=2) + "\n")
    return path
