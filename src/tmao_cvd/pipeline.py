"""End to end analysis: load, fit, evaluate, and write the results.

The order of operations follows the analysis plan in docs/analysis_plan.md.
Nothing in this module chooses between models or thresholds on the basis of
the results it computes. Every specification is fixed before the data are
seen, which is what allows the reported p value for the incremental
contribution of TMAO to be read at face value.
"""

from __future__ import annotations

import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn

from .config import DEFAULT_CONFIG, FIGURES_DIR, OUTCOME, TABLES_DIR, AnalysisConfig
from .data import load_cohort
from .evaluate import (
    auc_with_ci,
    bootstrap_interval,
    category_free_nri,
    delong_test,
    integrated_discrimination_improvement,
    summarise_calibration,
)
from .features import (
    BASELINE_FEATURES,
    CANDIDATE_MARKER,
    EXTENDED_FEATURES,
    add_derived_features,
    describe_by_outcome,
)
from .figures import plot_calibration, plot_decision_curve, plot_roc_curves
from .models import (
    MODEL_BUILDERS,
    cross_validated_predictions,
    logistic_coefficient_table,
    make_logistic_model,
)


def _performance_row(
    model_name: str, feature_set: str, y_true: np.ndarray, predictions, provenance: str
) -> dict[str, object]:
    """Assemble one row of the model performance table."""

    auc, lower, upper = auc_with_ci(y_true, predictions.averaged)
    calibration = summarise_calibration(y_true, predictions.averaged)

    per_repeat_auc = [
        float(auc_with_ci(y_true, repeat)[0]) for repeat in predictions.per_repeat
    ]

    return {
        "model": model_name,
        "feature_set": feature_set,
        "auc": round(auc, 4),
        "auc_ci_lower": round(lower, 4),
        "auc_ci_upper": round(upper, 4),
        "auc_sd_across_repeats": round(float(np.std(per_repeat_auc)), 4),
        "calibration_intercept": round(calibration.intercept, 4),
        "calibration_slope": round(calibration.slope, 4),
        "brier": round(calibration.brier, 4),
        "brier_skill": round(calibration.brier_skill, 4),
        "provenance": provenance,
    }


def run_analysis(
    config: AnalysisConfig | None = None, cohort_path: Path | None = None
) -> dict[str, object]:
    """Run the full analysis and write every table and figure.

    Returns a manifest describing the run. The manifest is also written to
    disk, because a figure without a record of the seed, the software
    versions and the data source that produced it cannot be reproduced and
    should not be trusted.
    """

    config = config or DEFAULT_CONFIG
    config.prepare_output_dirs()

    cohort, provenance = load_cohort(cohort_path)
    cohort = add_derived_features(cohort)
    outcome = cohort[OUTCOME].to_numpy(dtype=int)

    print(f"Cohort: {len(cohort)} participants, {int(outcome.sum())} events")
    print(f"Data source: {provenance}")

    # Table 1. Descriptive characteristics by outcome status.
    reported_columns = [column for column in cohort.columns if column != "participant_id"]
    table_one = describe_by_outcome(cohort[reported_columns], OUTCOME)
    table_one.to_csv(TABLES_DIR / "table_1_cohort_characteristics.csv", index=False)

    # Fit every model under both feature sets on identical folds.
    feature_sets = {"baseline": BASELINE_FEATURES, "extended": EXTENDED_FEATURES}
    predictions: dict[tuple[str, str], object] = {}
    performance_rows: list[dict[str, object]] = []

    for model_name, build_model in MODEL_BUILDERS.items():
        for feature_name, columns in feature_sets.items():
            result = cross_validated_predictions(
                build_model, cohort[columns], outcome, config
            )
            predictions[(model_name, feature_name)] = result
            performance_rows.append(
                _performance_row(model_name, feature_name, outcome, result, provenance)
            )
            print(
                f"  {model_name:<20} {feature_name:<9} "
                f"AUC {performance_rows[-1]['auc']:.3f}"
            )

    performance = pd.DataFrame(performance_rows)
    performance.to_csv(TABLES_DIR / "table_2_model_performance.csv", index=False)

    # Table 3. Incremental value of log TMAO, model by model.
    incremental_rows: list[dict[str, object]] = []
    for model_name in MODEL_BUILDERS:
        baseline = predictions[(model_name, "baseline")].averaged
        extended = predictions[(model_name, "extended")].averaged

        difference, p_value = delong_test(outcome, baseline, extended)
        idi = integrated_discrimination_improvement(outcome, baseline, extended)
        nri_events, nri_non_events, nri_total = category_free_nri(
            outcome, baseline, extended
        )

        idi_lower, idi_upper = bootstrap_interval(
            integrated_discrimination_improvement,
            outcome,
            baseline,
            extended,
            config.n_bootstrap,
            config.random_seed,
        )
        nri_lower, nri_upper = bootstrap_interval(
            category_free_nri,
            outcome,
            baseline,
            extended,
            config.n_bootstrap,
            config.random_seed + 1,
        )

        incremental_rows.append(
            {
                "model": model_name,
                "marker": CANDIDATE_MARKER,
                "delta_auc": round(difference, 4),
                "delong_p_value": round(p_value, 4),
                "idi": round(idi, 4),
                "idi_ci_lower": round(idi_lower, 4),
                "idi_ci_upper": round(idi_upper, 4),
                "nri_events": round(nri_events, 4),
                "nri_non_events": round(nri_non_events, 4),
                "nri_total": round(nri_total, 4),
                "nri_ci_lower": round(nri_lower, 4),
                "nri_ci_upper": round(nri_upper, 4),
                "provenance": provenance,
            }
        )

    incremental = pd.DataFrame(incremental_rows)
    incremental.to_csv(TABLES_DIR / "table_3_incremental_value.csv", index=False)

    # Table 4. Coefficients of the extended logistic model fitted to the
    # whole cohort. Reported for interpretation only. The performance
    # estimates above come from held out predictions, not from this fit.
    final_model = make_logistic_model(config.random_seed)
    final_model.fit(cohort[EXTENDED_FEATURES].to_numpy(dtype=float), outcome)
    coefficients = logistic_coefficient_table(final_model, cohort[EXTENDED_FEATURES])
    coefficients.to_csv(
        TABLES_DIR / "table_4_extended_model_coefficients.csv", index=False
    )

    # Figures, built from the primary model only.
    primary = {
        "baseline model": predictions[("logistic regression", "baseline")].averaged,
        "extended model (with TMAO)": predictions[
            ("logistic regression", "extended")
        ].averaged,
    }
    start, stop, count = config.decision_thresholds
    thresholds = np.linspace(start, stop, count)

    plot_roc_curves(
        outcome, primary, provenance, FIGURES_DIR / "figure_1_roc_curves.png", config.figure_dpi
    )
    plot_calibration(
        outcome,
        primary,
        provenance,
        FIGURES_DIR / "figure_2_calibration.png",
        config.figure_dpi,
        config.calibration_bins,
    )
    plot_decision_curve(
        outcome,
        primary,
        thresholds,
        provenance,
        FIGURES_DIR / "figure_3_decision_curve.png",
        config.figure_dpi,
    )

    manifest = {
        "run_completed_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "data_provenance": provenance,
        "n_participants": int(len(cohort)),
        "n_events": int(outcome.sum()),
        "event_rate": round(float(outcome.mean()), 4),
        "random_seed": config.random_seed,
        "n_splits": config.n_splits,
        "n_repeats": config.n_repeats,
        "n_bootstrap": config.n_bootstrap,
        "baseline_features": BASELINE_FEATURES,
        "candidate_marker": CANDIDATE_MARKER,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "scikit_learn_version": sklearn.__version__,
    }
    (TABLES_DIR / "run_manifest.json").write_text(json.dumps(manifest, indent=2))

    print(f"\nTables written to {TABLES_DIR}")
    print(f"Figures written to {FIGURES_DIR}")
    return manifest
