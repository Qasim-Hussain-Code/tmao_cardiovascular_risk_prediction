"""Tests for the derived variables and the schema check."""

from __future__ import annotations

import numpy as np
import pytest

from tmao_cvd.data import SchemaError, validate_cohort
from tmao_cvd.features import BASELINE_FEATURES, EXTENDED_FEATURES, add_derived_features
from tmao_cvd.simulate import SimulationParameters, simulate_cohort


def test_extended_model_adds_exactly_one_term():
    # The whole comparison rests on the two feature sets differing by the
    # candidate marker and nothing else.
    assert set(EXTENDED_FEATURES) - set(BASELINE_FEATURES) == {"log_tmao"}
    assert len(EXTENDED_FEATURES) == len(BASELINE_FEATURES) + 1


def test_derived_features_do_not_modify_the_input():
    cohort = simulate_cohort(SimulationParameters(n_participants=100, seed=2))
    before = cohort.columns.tolist()
    derived = add_derived_features(cohort)
    assert cohort.columns.tolist() == before
    assert "log_tmao" in derived.columns
    assert np.allclose(np.exp(derived["log_tmao"]), derived["tmao_umol_l"])


def test_validation_rejects_an_implausible_measurement():
    cohort = simulate_cohort(SimulationParameters(n_participants=100, seed=3))
    cohort.loc[0, "egfr_ml_min_1_73m2"] = 999.0
    with pytest.raises(SchemaError, match="egfr"):
        validate_cohort(cohort)


def test_validation_reports_every_problem_at_once():
    cohort = simulate_cohort(SimulationParameters(n_participants=100, seed=4))
    cohort = cohort.drop(columns=["hs_crp_mg_l", "diabetes"])
    with pytest.raises(SchemaError) as error:
        validate_cohort(cohort)
    assert "hs_crp_mg_l" in str(error.value)
    assert "diabetes" in str(error.value)


def test_validation_rejects_duplicated_identifiers():
    cohort = simulate_cohort(SimulationParameters(n_participants=50, seed=5))
    cohort.loc[1, "participant_id"] = cohort.loc[0, "participant_id"]
    with pytest.raises(SchemaError, match="duplicated"):
        validate_cohort(cohort)
