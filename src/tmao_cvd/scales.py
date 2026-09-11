"""Measurement scale as a typed, enforced property of the data.

Metabolite panels arrive on two very different scales. Targeted assays with
calibration curves report physical concentrations, usually micromoles per
litre. Most deposited metabolomics reports relative peak area, an arbitrary
instrument response that is comparable within a run and meaningless outside
it.

The distinction matters because the literature is full of TMAO cut points
expressed in micromoles per litre, and applying one of those to peak area
data produces a number rather than an error. Nothing about the arithmetic
objects. The result is simply wrong, and wrong in a way that survives review
because it looks like every other threshold comparison.

This module makes the scale explicit and makes the illegal operation raise.
The reasoning behind that choice is recorded in section 4 of
docs/reanalysis_plan.md: a comment relies on the next reader having read it,
and an exception does not.
"""

from __future__ import annotations

from enum import Enum

import pandas as pd


class MeasurementScale(Enum):
    """The scale on which a set of metabolite values is expressed."""

    #: Arbitrary instrument response. Comparable within a run, not outside it.
    RELATIVE_PEAK_AREA = "relative peak area (arbitrary instrument units)"
    #: A physical concentration, comparable against external reference values.
    MICROMOLES_PER_LITRE = "micromoles per litre"

    @property
    def is_concentration(self) -> bool:
        return self is MeasurementScale.MICROMOLES_PER_LITRE


class MeasurementScaleError(RuntimeError):
    """Raised when an operation is attempted on the wrong measurement scale."""


#: No plausible concentration in micromoles per litre reaches this magnitude,
#: so a median above it is positive evidence that values are instrument
#: response. Used to catch a calibrated file being substituted for a peak area
#: one without the declared scale being updated to match.
PEAK_AREA_PLAUSIBILITY_FLOOR = 1.0e4


def require_concentration_scale(scale: MeasurementScale, operation: str) -> None:
    """Raise unless ``scale`` is a physical concentration.

    Call this at the top of anything that compares a measurement against an
    externally derived value: a clinical cut point, a tertile boundary from
    another cohort, a published reference range.
    """

    if scale.is_concentration:
        return

    raise MeasurementScaleError(
        f"{operation} requires values in {MeasurementScale.MICROMOLES_PER_LITRE.value}, "
        f"but these values are on the {scale.value} scale. Peak areas are arbitrary "
        "instrument units and carry no physical meaning outside the run that produced "
        "them, so no external threshold, cut point or reference range may be applied "
        "to them. See section 4 of docs/reanalysis_plan.md."
    )


def apply_concentration_threshold(
    values: pd.Series, threshold_umol_per_litre: float, scale: MeasurementScale
) -> pd.Series:
    """Dichotomise ``values`` at an externally derived concentration threshold.

    The only way to apply a published cut point inside this codebase, and it
    refuses to run on peak area data. There is no bypass argument on purpose.
    """

    require_concentration_scale(
        scale, f"applying the concentration threshold {threshold_umol_per_litre}"
    )
    return values >= threshold_umol_per_litre


def assert_peak_area_scale(
    frame: pd.DataFrame, columns: list[str], scale: MeasurementScale
) -> None:
    """Check that data declared as peak area actually looks like peak area.

    The declaration in the loader is an assertion about the file, and a file
    can be replaced. This verifies the claim against the values themselves, so
    that swapping in a concentration calibrated export fails at load time
    rather than silently changing what every downstream number means.
    """

    if scale is not MeasurementScale.RELATIVE_PEAK_AREA:
        return

    for column in columns:
        median = float(pd.to_numeric(frame[column], errors="coerce").median())
        if median < PEAK_AREA_PLAUSIBILITY_FLOOR:
            raise MeasurementScaleError(
                f"column '{column}' is declared as {scale.value} but its median is "
                f"{median:.4g}, below the plausibility floor of "
                f"{PEAK_AREA_PLAUSIBILITY_FLOOR:.0e}. Either the declared scale is "
                "wrong or the file has been replaced with calibrated concentrations. "
                "Resolve this before any result from it is used."
            )
