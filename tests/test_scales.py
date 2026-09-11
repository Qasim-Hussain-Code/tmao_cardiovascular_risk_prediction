"""Tests for the measurement scale guard.

The guard exists because applying a published concentration cut point to peak
area data produces a plausible number rather than an error. These tests check
that the illegal operation actually raises, and that there is no way around it,
since a guard with a bypass argument is a comment with extra steps.
"""

from __future__ import annotations

import inspect

import pandas as pd
import pytest

from tmao_cvd.scales import (
    PEAK_AREA_PLAUSIBILITY_FLOOR,
    MeasurementScale,
    MeasurementScaleError,
    apply_concentration_threshold,
    assert_peak_area_scale,
    require_concentration_scale,
)

PEAK_AREA = MeasurementScale.RELATIVE_PEAK_AREA
CONCENTRATION = MeasurementScale.MICROMOLES_PER_LITRE


def test_concentration_scale_is_permitted():
    require_concentration_scale(CONCENTRATION, "a threshold comparison")


def test_peak_area_scale_is_refused():
    with pytest.raises(MeasurementScaleError, match="micromoles per litre"):
        require_concentration_scale(PEAK_AREA, "a threshold comparison")


def test_published_cut_point_cannot_be_applied_to_peak_areas():
    # 6.18 micromoles per litre is a cut point that appears in the TMAO
    # literature. On peak area data it must not silently produce an answer.
    values = pd.Series([7.7e7, 9.4e7, 1.2e8])
    with pytest.raises(MeasurementScaleError):
        apply_concentration_threshold(values, 6.18, PEAK_AREA)


def test_cut_point_works_on_genuine_concentrations():
    values = pd.Series([3.1, 6.5, 9.0])
    flagged = apply_concentration_threshold(values, 6.18, CONCENTRATION)
    assert flagged.tolist() == [False, True, True]


def test_guard_has_no_bypass_argument():
    # If someone adds force= or skip_check= later, this fails and makes them
    # justify it in review rather than in a commit nobody reads.
    parameters = set(inspect.signature(apply_concentration_threshold).parameters)
    assert parameters == {"values", "threshold_umol_per_litre", "scale"}


def test_declared_peak_area_scale_is_verified_against_the_values():
    # A calibrated file substituted for a peak area one must fail at load,
    # not quietly change what every downstream number means.
    calibrated = pd.DataFrame({"Trimethylamine N-oxide": [3.2, 4.1, 5.5]})
    with pytest.raises(MeasurementScaleError, match="plausibility floor"):
        assert_peak_area_scale(calibrated, ["Trimethylamine N-oxide"], PEAK_AREA)


def test_genuine_peak_areas_pass_verification():
    peaks = pd.DataFrame({"Trimethylamine N-oxide": [7.7e7, 9.4e7, 8.0e7]})
    assert_peak_area_scale(peaks, ["Trimethylamine N-oxide"], PEAK_AREA)


def test_plausibility_floor_sits_above_any_real_concentration():
    # TMAO in micromoles per litre reaches perhaps a few hundred in extreme
    # renal failure. The floor must sit far above that to be a safe test.
    assert PEAK_AREA_PLAUSIBILITY_FLOOR > 1000
