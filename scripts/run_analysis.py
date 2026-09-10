#!/usr/bin/env python3
"""Command line entry point for the full analysis.

Run from the repository root:

    python scripts/run_analysis.py

If a cohort export is present at data/raw/cohort.csv it is used. Otherwise
the pipeline simulates a cohort and says so, loudly, in every table, figure
and log line it produces.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from tmao_cvd.config import DEFAULT_CONFIG, AnalysisConfig  # noqa: E402
from tmao_cvd.pipeline import run_analysis  # noqa: E402


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cohort",
        type=Path,
        default=None,
        help="path to a cohort CSV, overriding data/raw/cohort.csv",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_CONFIG.random_seed,
        help="random seed for simulation, resampling and the bootstrap",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=DEFAULT_CONFIG.n_repeats,
        help="number of cross validation repeats",
    )
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=DEFAULT_CONFIG.n_bootstrap,
        help="bootstrap resamples for the reclassification intervals",
    )
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    config = AnalysisConfig(
        random_seed=arguments.seed,
        n_repeats=arguments.repeats,
        n_bootstrap=arguments.bootstrap,
    )
    run_analysis(config=config, cohort_path=arguments.cohort)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
