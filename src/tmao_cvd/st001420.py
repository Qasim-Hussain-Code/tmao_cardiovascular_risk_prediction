"""Loader for Metabolomics Workbench study ST001420.

The study is described in section 3 of docs/reanalysis_plan.md. In short: 750
patients with stable angina, sampled 48 hours after percutaneous coronary
intervention and followed for nine months, 210 of whom had recurrent angina,
with 600 metabolites measured by targeted LC-MS/MS. Deposited under CC BY 4.0.

Two properties of this deposit are enforced here rather than documented and
hoped for.

First, the values are relative peak areas, not concentrations. The loader
declares that scale and verifies it against the data, so that a calibrated
file substituted later fails at load time instead of silently changing the
meaning of every downstream number.

Second, sample identifiers are perfectly confounded with outcome: S1 to S210
are exactly the cases and S211 to S750 exactly the controls, and no batch, run
order or acquisition date was deposited. The loader therefore exposes
``sample_index`` as a column rather than discarding it, because the diagnostic
in section 6 of the plan needs it. It is emphatically not a predictor, and
:data:`METABOLITE_COLUMNS` excludes it.
"""

from __future__ import annotations

import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .config import RAW_DATA_DIR
from .scales import MeasurementScale, assert_peak_area_scale

STUDY_ID = "ST001420"
MWTAB_URL = (
    f"https://www.metabolomicsworkbench.org/rest/study/study_id/{STUDY_ID}/mwtab/txt"
)
LOCAL_MWTAB = RAW_DATA_DIR / f"{STUDY_ID}.mwtab.txt"

#: Everything in this deposit is instrument response. See scales.py.
SCALE = MeasurementScale.RELATIVE_PEAK_AREA

OUTCOME = "recurrent_angina_9mo"
SAMPLE_INDEX = "sample_index"

TMAO = "Trimethylamine N-oxide"
#: The dietary and microbial precursors of TMAO, measured in the same run.
#: These form the baseline model for question 1 of the reanalysis plan.
PRECURSORS = ["Choline", "Betaine", "Carnitine"]
PATHWAY = [TMAO] + PRECURSORS

#: Expected shape of the deposit, checked on every load. These are not
#: defensive guesses; they are the numbers stated in the study record, and a
#: mismatch means the deposit has been revised and the plan needs rereading.
EXPECTED_SAMPLES = 750
EXPECTED_METABOLITES = 600
EXPECTED_EVENTS = 210


class DepositError(RuntimeError):
    """Raised when the deposit does not match its documented structure."""


@dataclass(frozen=True)
class DepositMetadata:
    """Provenance travelling with the loaded matrix."""

    study_id: str
    scale: MeasurementScale
    n_participants: int
    n_events: int
    n_metabolites: int
    source: str
    licence: str = "CC BY 4.0"

    @property
    def provenance(self) -> str:
        """The label stamped onto every table and figure built from this."""
        return (
            f"{self.study_id} (observed cohort, {self.n_participants} participants, "
            f"{self.n_events} events; values are {self.scale.value})"
        )


def fetch_mwtab(destination: Path | None = None, force: bool = False) -> Path:
    """Download the mwTab record unless it is already present.

    The file is written under data/raw/, which is not tracked, so the data are
    never committed while remaining one command away from any clone.
    """

    destination = destination or LOCAL_MWTAB
    destination.parent.mkdir(parents=True, exist_ok=True)

    if destination.exists() and not force and destination.stat().st_size > 0:
        return destination

    with urllib.request.urlopen(MWTAB_URL, timeout=180) as response:
        payload = response.read()

    if len(payload) < 100_000:
        raise DepositError(
            f"download from {MWTAB_URL} returned only {len(payload)} bytes, which is "
            "far short of the expected record. Refusing to write a truncated file."
        )

    destination.write_bytes(payload)
    return destination


def _read_data_block(text: str) -> tuple[list[str], list[str], list[list[str]]]:
    """Extract the sample names, factor labels and metabolite rows."""

    lines = text.split("\n")
    try:
        start = next(
            i for i, line in enumerate(lines) if line.startswith("MS_METABOLITE_DATA_START")
        )
        end = next(
            i for i, line in enumerate(lines) if line.startswith("MS_METABOLITE_DATA_END")
        )
    except StopIteration as error:
        raise DepositError("no MS_METABOLITE_DATA block found in the record") from error

    block = lines[start + 1 : end]
    samples = block[0].split("\t")[1:]
    factors = block[1].split("\t")[1:]
    rows = [line.split("\t") for line in block[2:] if line.strip()]
    return samples, factors, rows


def load_st001420(path: Path | None = None) -> tuple[pd.DataFrame, DepositMetadata]:
    """Return the ST001420 matrix and its provenance.

    The returned frame is indexed by sample identifier and carries the 600
    metabolite columns, the binary outcome, and the sample index needed by the
    acquisition order diagnostic.
    """

    path = fetch_mwtab(path)
    samples, factors, rows = _read_data_block(
        path.read_text(encoding="utf-8", errors="replace")
    )

    frame = pd.DataFrame(
        {row[0]: pd.to_numeric(pd.Series(row[1:]), errors="coerce").to_numpy() for row in rows},
        index=pd.Index(samples, name="sample_id"),
    )

    distinct = set(factors)
    if distinct != {"Disease:Angina", "Disease:Angina-free"}:
        raise DepositError(f"unexpected outcome labels in the deposit: {sorted(distinct)}")

    frame[OUTCOME] = [1 if f == "Disease:Angina" else 0 for f in factors]
    frame[SAMPLE_INDEX] = [int(s.lstrip("S")) for s in samples]

    metabolites = [c for c in frame.columns if c not in (OUTCOME, SAMPLE_INDEX)]
    _validate(frame, metabolites)

    metadata = DepositMetadata(
        study_id=STUDY_ID,
        scale=SCALE,
        n_participants=int(len(frame)),
        n_events=int(frame[OUTCOME].sum()),
        n_metabolites=len(metabolites),
        source=MWTAB_URL,
    )
    return frame, metadata


def _validate(frame: pd.DataFrame, metabolites: list[str]) -> None:
    """Check the deposit against its documented structure and declared scale."""

    problems: list[str] = []

    if len(frame) != EXPECTED_SAMPLES:
        problems.append(f"expected {EXPECTED_SAMPLES} participants, found {len(frame)}")
    if len(metabolites) != EXPECTED_METABOLITES:
        problems.append(
            f"expected {EXPECTED_METABOLITES} metabolites, found {len(metabolites)}"
        )

    events = int(frame[OUTCOME].sum())
    if events != EXPECTED_EVENTS:
        problems.append(f"expected {EXPECTED_EVENTS} events, found {events}")

    for name in PATHWAY:
        if name not in frame.columns:
            problems.append(f"pathway metabolite '{name}' is absent")
        elif frame[name].isna().any():
            n_missing = int(frame[name].isna().sum())
            problems.append(f"pathway metabolite '{name}' has {n_missing} missing values")

    if problems:
        raise DepositError(
            "the deposit does not match its documented structure:\n  - "
            + "\n  - ".join(problems)
        )

    # Verifies the declared scale against the values themselves.
    assert_peak_area_scale(frame, PATHWAY, SCALE)


def metabolite_columns(frame: pd.DataFrame) -> list[str]:
    """The 600 metabolite columns, excluding the outcome and the sample index.

    Kept as a function rather than a constant so that ``sample_index`` cannot
    reach a model by being left in a hand written column list.
    """

    return [c for c in frame.columns if c not in (OUTCOME, SAMPLE_INDEX)]


def add_log_pathway_features(frame: pd.DataFrame) -> pd.DataFrame:
    """Add natural log transformed pathway variables.

    Peak areas span orders of magnitude and are right skewed, so the log scale
    is the sensible one for a linear term. The transform is monotone, so it
    changes nothing about the arbitrary units: a logged peak area remains
    uncomparable with any external concentration.
    """

    out = frame.copy()
    out["log_tmao"] = np.log(out[TMAO].astype(float))
    for precursor in PRECURSORS:
        out[f"log_{precursor.lower()}"] = np.log(out[precursor].astype(float))
    return out


#: Baseline model for question 1: the precursors alone.
PRECURSOR_FEATURES = [f"log_{p.lower()}" for p in PRECURSORS]
#: Extended model for question 1: precursors plus the product under test.
EXTENDED_FEATURES = PRECURSOR_FEATURES + ["log_tmao"]
