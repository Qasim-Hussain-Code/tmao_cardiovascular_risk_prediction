"""Configuration shared across the analysis.

Analysis settings are collected here rather than scattered through the
modules that use them. The reason is auditability: a reviewer asking which
seed produced a figure, or how many cross validation repeats stand behind a
confidence interval, should find the answer in one file.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_DIR = PROJECT_ROOT / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
PROCESSED_DATA_DIR = DATA_DIR / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
TABLES_DIR = RESULTS_DIR / "tables"

#: Location the pipeline checks for a real cohort before falling back to
#: simulation. The file is never committed. See data/README.md.
COHORT_FILE = RAW_DATA_DIR / "cohort.csv"

OUTCOME = "mace_3yr"
ID_COLUMN = "participant_id"


@dataclass(frozen=True)
class AnalysisConfig:
    """Settings that govern model fitting and evaluation.

    Attributes
    ----------
    random_seed:
        Seed for every stochastic component, including the simulation, the
        cross validation splits and the bootstrap.
    n_splits:
        Number of folds in each stratified cross validation round.
    n_repeats:
        Number of times the cross validation is repeated with different
        splits. Repetition matters here because a single five fold split of
        a cohort with roughly four hundred events gives noticeably unstable
        estimates of the area under the curve.
    n_bootstrap:
        Bootstrap resamples used for the reclassification statistics, which
        have no convenient closed form variance.
    decision_thresholds:
        Risk thresholds spanned by the decision curve. The range is chosen
        to cover the region where a clinician might plausibly change
        management for primary prevention.
    """

    random_seed: int = 20260910
    n_splits: int = 5
    n_repeats: int = 5
    n_bootstrap: int = 1000
    decision_thresholds: tuple[float, float, int] = (0.01, 0.40, 79)
    calibration_bins: int = 10
    figure_dpi: int = 300
    output_dirs: tuple[Path, ...] = field(
        default_factory=lambda: (FIGURES_DIR, TABLES_DIR)
    )

    def prepare_output_dirs(self) -> None:
        """Create the results directories if a run needs them."""
        for directory in self.output_dirs:
            directory.mkdir(parents=True, exist_ok=True)


DEFAULT_CONFIG = AnalysisConfig()
