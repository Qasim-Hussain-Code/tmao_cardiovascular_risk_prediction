"""End to end check that the modelling code recovers a known effect.

The individual units are tested elsewhere. This test covers the join between
them, which is where mistakes in an analysis pipeline usually live: a
transformation applied twice, a feature set built from the wrong columns, a
scaler fitted on the wrong axis. None of those would be caught by a test of
any single function.

The simulator plants a known log odds ratio per standard deviation for each
covariate. If simulation, feature construction, scaling and estimation are
mutually consistent, an unpenalised fit to a large synthetic cohort must
return those values up to sampling error. If any stage disagrees with the
others, the estimates drift and this test fails.
"""

from __future__ import annotations

import pytest

from tmao_cvd.config import OUTCOME
from tmao_cvd.features import EXTENDED_FEATURES, add_derived_features
from tmao_cvd.models import logistic_coefficient_table, make_logistic_model
from tmao_cvd.simulate import SimulationParameters, simulate_cohort


@pytest.fixture(scope="module")
def fitted():
    # Large enough that the standard errors are near 0.02, so a tolerance of
    # 0.05 leaves room for sampling error without hiding a real defect.
    parameters = SimulationParameters(n_participants=50000, seed=99)
    cohort = add_derived_features(simulate_cohort(parameters))
    outcome = cohort[OUTCOME].to_numpy(int)

    model = make_logistic_model(0)
    model.fit(cohort[EXTENDED_FEATURES].to_numpy(float), outcome)
    table = logistic_coefficient_table(model, cohort[EXTENDED_FEATURES])
    return parameters, table.set_index("term")


@pytest.mark.parametrize("term", EXTENDED_FEATURES)
def test_planted_log_odds_are_recovered(fitted, term):
    parameters, table = fitted
    planted = parameters.outcome_log_odds[term]
    assert table.loc[term, "log_odds_per_sd"] == pytest.approx(planted, abs=0.05)


def test_marker_odds_ratio_is_reported_on_the_odds_scale(fitted):
    _, table = fitted
    row = table.loc["log_tmao"]
    assert row["ci_lower"] < row["odds_ratio_per_sd"] < row["ci_upper"]
    assert row["ci_lower"] > 1.0
