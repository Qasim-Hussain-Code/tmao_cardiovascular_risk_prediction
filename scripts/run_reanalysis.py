#!/usr/bin/env python3
"""Run the three pre-specified questions of docs/reanalysis_plan.md.

    python scripts/run_reanalysis.py

Downloads ST001420 from Metabolomics Workbench on first run, into data/raw/,
which is not tracked. Every verdict is decided by thresholds fixed in the plan
before any estimate existed, and the claim scope is printed above the results
by the reporting code rather than left to the reader to look up.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from tmao_cvd.config import DEFAULT_CONFIG, TABLES_DIR, AnalysisConfig  # noqa: E402
from tmao_cvd import post_hoc  # noqa: E402
from tmao_cvd.reanalysis import question_one, question_three, question_two  # noqa: E402
from tmao_cvd.reporting import SCOPE_STATEMENT, render_report  # noqa: E402
from tmao_cvd.st001420 import add_log_pathway_features, load_st001420  # noqa: E402


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=DEFAULT_CONFIG.random_seed)
    parser.add_argument("--repeats", type=int, default=DEFAULT_CONFIG.n_repeats)
    parser.add_argument("--bootstrap", type=int, default=DEFAULT_CONFIG.n_bootstrap)
    parser.add_argument(
        "--skip-question-two",
        action="store_true",
        help="skip the 600 metabolite nested cross validation, which is the slow part",
    )
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    config = AnalysisConfig(
        random_seed=arguments.seed,
        n_repeats=arguments.repeats,
        n_bootstrap=arguments.bootstrap,
    )
    config.prepare_output_dirs()

    frame, metadata = load_st001420()
    frame = add_log_pathway_features(frame)
    print(f"Loaded {metadata.provenance}\n")

    results = []

    print("Running question 1 ...")
    q1, curves = question_one(frame, config)
    results.append(q1)

    if not arguments.skip_question_two:
        print("Running question 2 (nested cross validation over 600 metabolites) ...")
        results.append(question_two(frame, config))

    print("Running question 3 ...")
    results.append(question_three(frame))

    report = render_report(results, scope=SCOPE_STATEMENT)
    print("\n" + report)

    # Post hoc, and labelled as such everywhere it appears. Triggered by the
    # question 2 result, not planned in advance. See section 8 of the plan.
    print("\n" + "=" * 78)
    print("POST HOC DIAGNOSTICS (not pre-specified; prompted by the question 2 result)")
    print("=" * 78)
    diagnostics = post_hoc.run_all(frame)
    for name, table in diagnostics.items():
        print(f"\n{name}:")
        print(table.to_string(index=False))
        table.to_csv(TABLES_DIR / f"post_hoc_{name}.csv", index=False)

    (TABLES_DIR / "reanalysis_report.txt").write_text(report)
    for index, result in enumerate(results, start=1):
        result.table.to_csv(TABLES_DIR / f"reanalysis_question_{index}.csv", index=False)

    print(f"\nWritten to {TABLES_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
