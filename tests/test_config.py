"""Tests that the recorded settings match the settings that actually run.

config/analysis_config.yaml exists so that a reader can see how a run was
configured without reading Python. That convenience creates a hazard: two
records of the same numbers can drift apart, and the written one is the one
a reader will trust.

These tests make drift impossible to miss. The dataclass remains the single
source of truth for what executes; the YAML file is held to it.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from tmao_cvd.config import DEFAULT_CONFIG, OUTCOME
from tmao_cvd.features import CANDIDATE_MARKER

CONFIG_FILE = Path(__file__).resolve().parents[1] / "config" / "analysis_config.yaml"


@pytest.fixture(scope="module")
def recorded():
    return yaml.safe_load(CONFIG_FILE.read_text())


def test_seed_matches(recorded):
    assert recorded["random_seed"] == DEFAULT_CONFIG.random_seed


@pytest.mark.parametrize("setting", ["n_splits", "n_repeats", "n_bootstrap"])
def test_resampling_settings_match(recorded, setting):
    assert recorded["resampling"][setting] == getattr(DEFAULT_CONFIG, setting)


def test_decision_curve_settings_match(recorded):
    start, stop, count = DEFAULT_CONFIG.decision_thresholds
    assert recorded["decision_curve"]["threshold_min"] == pytest.approx(start)
    assert recorded["decision_curve"]["threshold_max"] == pytest.approx(stop)
    assert recorded["decision_curve"]["n_thresholds"] == count


@pytest.mark.parametrize("setting", ["calibration_bins", "figure_dpi"])
def test_reporting_settings_match(recorded, setting):
    assert recorded["reporting"][setting] == getattr(DEFAULT_CONFIG, setting)


def test_outcome_and_marker_names_match(recorded):
    # A mismatch here would be read by a reviewer as the study having tested
    # a different variable than the one it reports.
    assert recorded["outcome"] == OUTCOME
    assert recorded["candidate_marker"] == CANDIDATE_MARKER
