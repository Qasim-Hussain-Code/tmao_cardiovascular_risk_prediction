"""Tests for the synthetic cohort generator.

These tests check that the simulator produces data the rest of the pipeline
can consume and that its stated properties hold. They do not check that the
simulated cohort resembles any real population, because it is not intended
to.
"""

from __future__ import annotations

import numpy as np

from tmao_cvd.data import validate_cohort
from tmao_cvd.simulate import SimulationParameters, simulate_cohort


def test_simulated_cohort_satisfies_the_schema():
    cohort = simulate_cohort(SimulationParameters(n_participants=800, seed=1))
    validate_cohort(cohort)
    assert len(cohort) == 800


def test_event_rate_matches_the_requested_target():
    # The intercept is solved rather than fixed, so the realised rate should
    # sit close to the target up to binomial sampling error.
    parameters = SimulationParameters(n_participants=20000, seed=7, target_event_rate=0.10)
    cohort = simulate_cohort(parameters)
    assert abs(cohort["mace_3yr"].mean() - 0.10) < 0.01


def test_simulation_is_reproducible_from_its_seed():
    first = simulate_cohort(SimulationParameters(n_participants=500, seed=42))
    second = simulate_cohort(SimulationParameters(n_participants=500, seed=42))
    assert first.equals(second)


def test_tmao_rises_as_renal_function_falls():
    # The generative model makes TMAO depend on eGFR because the metabolite
    # is renally cleared. If that dependence disappeared the simulation would
    # no longer pose the confounding problem the analysis is built to handle.
    cohort = simulate_cohort(SimulationParameters(n_participants=5000, seed=11))
    correlation = np.corrcoef(cohort["egfr_ml_min_1_73m2"], np.log(cohort["tmao_umol_l"]))[0, 1]
    assert correlation < -0.2
