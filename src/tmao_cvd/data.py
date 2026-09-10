"""Loading and validation of the analysis dataset.

The pipeline accepts either a real cohort export or a simulated one, and it
holds both to the same schema. Validation runs before anything else because
a silently miscoded column is far more expensive to discover at the point of
interpreting a coefficient than at the point of reading the file.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import COHORT_FILE, ID_COLUMN, OUTCOME
from .simulate import SimulationParameters, simulate_cohort

#: Expected columns, with the plausible range each measurement may take.
#: The bounds are wide on purpose. They are there to catch unit errors and
#: sentinel codes such as 999, not to exclude unusual but genuine values.
SCHEMA: dict[str, tuple[float, float]] = {
    "age_years": (18.0, 110.0),
    "sex_male": (0.0, 1.0),
    "current_smoker": (0.0, 1.0),
    "diabetes": (0.0, 1.0),
    "bmi_kg_m2": (10.0, 70.0),
    "systolic_bp_mmhg": (60.0, 260.0),
    "total_cholesterol_mmol_l": (1.0, 15.0),
    "hdl_cholesterol_mmol_l": (0.2, 5.0),
    "egfr_ml_min_1_73m2": (5.0, 200.0),
    "hs_crp_mg_l": (0.01, 200.0),
    "tmao_umol_l": (0.05, 500.0),
    OUTCOME: (0.0, 1.0),
}


class SchemaError(ValueError):
    """Raised when a dataset does not satisfy the analysis schema."""


def validate_cohort(frame: pd.DataFrame) -> None:
    """Check a cohort against :data:`SCHEMA`.

    All problems are collected before raising rather than failing on the
    first one, so that a user fixing an export sees the full list in one
    pass instead of discovering it one column at a time.
    """

    problems: list[str] = []

    if ID_COLUMN not in frame.columns:
        problems.append(f"missing identifier column '{ID_COLUMN}'")
    elif frame[ID_COLUMN].duplicated().any():
        n_duplicated = int(frame[ID_COLUMN].duplicated().sum())
        problems.append(f"{n_duplicated} duplicated values in '{ID_COLUMN}'")

    for column, (lower, upper) in SCHEMA.items():
        if column not in frame.columns:
            problems.append(f"missing column '{column}'")
            continue

        values = pd.to_numeric(frame[column], errors="coerce")
        if values.isna().all():
            problems.append(f"column '{column}' holds no numeric values")
            continue

        out_of_range = values.dropna()
        out_of_range = out_of_range[(out_of_range < lower) | (out_of_range > upper)]
        if not out_of_range.empty:
            problems.append(
                f"column '{column}' has {len(out_of_range)} values outside "
                f"the plausible range [{lower}, {upper}]"
            )

    if OUTCOME in frame.columns:
        observed = set(pd.to_numeric(frame[OUTCOME], errors="coerce").dropna().unique())
        if not observed.issubset({0.0, 1.0}):
            problems.append(f"outcome '{OUTCOME}' is not coded as 0 and 1")
        elif len(observed) < 2:
            problems.append(f"outcome '{OUTCOME}' has no variation")

    if problems:
        raise SchemaError(
            "the cohort does not satisfy the analysis schema:\n  - "
            + "\n  - ".join(problems)
        )


def load_cohort(path: Path | None = None) -> tuple[pd.DataFrame, str]:
    """Return the analysis cohort and a label describing its provenance.

    A real export at :data:`~tmao_cvd.config.COHORT_FILE` takes precedence.
    If none is present the function falls back to simulation, and the
    returned provenance label says so. Every table and figure the pipeline
    writes carries that label, so a simulated result cannot be mistaken for
    an empirical one further downstream.
    """

    path = path or COHORT_FILE

    if path.exists():
        frame = pd.read_csv(path)
        provenance = f"observed cohort ({path.name})"
    else:
        frame = simulate_cohort(SimulationParameters())
        provenance = "SIMULATED DATA, not an observed cohort"

    validate_cohort(frame)
    return frame, provenance
