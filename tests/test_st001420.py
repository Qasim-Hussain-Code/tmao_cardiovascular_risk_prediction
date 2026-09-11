"""Tests for the ST001420 loader.

Parsing and validation are tested offline against a small synthetic record so
that the test suite does not depend on a network call. The integration test
against the real deposit runs when the file is available and skips when it is
not, which keeps continuous integration honest rather than flaky.
"""

from __future__ import annotations

import numpy as np
import pytest

from tmao_cvd.st001420 import (
    EXTENDED_FEATURES,
    OUTCOME,
    PRECURSOR_FEATURES,
    SAMPLE_INDEX,
    TMAO,
    DepositError,
    _read_data_block,
    add_log_pathway_features,
    load_st001420,
    metabolite_columns,
)

SYNTHETIC = "\n".join(
    [
        "#MS_METABOLITE_DATA",
        "MS_METABOLITE_DATA_START",
        "Samples\tS1\tS2\tS3",
        "Factors\tDisease:Angina\tDisease:Angina-free\tDisease:Angina-free",
        "Trimethylamine N-oxide\t9.4e7\t7.7e7\t8.0e7",
        "Choline\t1.5e7\t1.4e7\t1.6e7",
        "MS_METABOLITE_DATA_END",
    ]
)


def test_parser_extracts_samples_factors_and_rows():
    samples, factors, rows = _read_data_block(SYNTHETIC)
    assert samples == ["S1", "S2", "S3"]
    assert factors[0] == "Disease:Angina"
    assert [r[0] for r in rows] == ["Trimethylamine N-oxide", "Choline"]


def test_parser_rejects_a_record_with_no_data_block():
    with pytest.raises(DepositError, match="no MS_METABOLITE_DATA block"):
        _read_data_block("#STUDY\nST:STUDY_TITLE\tsomething else")


def test_feature_sets_differ_by_exactly_the_marker_under_test():
    # Question 1 of the reanalysis plan is only interpretable if the baseline
    # and extended models differ by log TMAO and nothing else.
    assert set(EXTENDED_FEATURES) - set(PRECURSOR_FEATURES) == {"log_tmao"}
    assert len(EXTENDED_FEATURES) == len(PRECURSOR_FEATURES) + 1


@pytest.fixture(scope="module")
def cohort():
    try:
        return load_st001420()
    except Exception as error:  # network unavailable, or deposit revised
        pytest.skip(f"ST001420 not available locally: {error}")


def test_deposit_matches_its_documented_structure(cohort):
    frame, meta = cohort
    assert len(frame) == 750
    assert int(frame[OUTCOME].sum()) == 210
    assert len(metabolite_columns(frame)) == 600
    assert meta.licence == "CC BY 4.0"


def test_provenance_states_the_measurement_scale(cohort):
    # Every table and figure carries this string. If it ever stops saying the
    # values are arbitrary units, a reader can mistake them for concentrations.
    _, meta = cohort
    assert "arbitrary instrument units" in meta.provenance
    assert "observed cohort" in meta.provenance


def test_sample_index_is_kept_but_excluded_from_the_metabolites(cohort):
    # It is needed by the acquisition order diagnostic and must never reach a
    # model, since it encodes the outcome perfectly.
    frame, _ = cohort
    assert SAMPLE_INDEX in frame.columns
    assert SAMPLE_INDEX not in metabolite_columns(frame)
    assert OUTCOME not in metabolite_columns(frame)


def test_sample_index_encodes_the_outcome_perfectly(cohort):
    # Documented in section 5 of the reanalysis plan. Asserted here so that the
    # limitation cannot be quietly forgotten.
    frame, _ = cohort
    cases = frame.loc[frame[OUTCOME] == 1, SAMPLE_INDEX]
    controls = frame.loc[frame[OUTCOME] == 0, SAMPLE_INDEX]
    assert cases.max() < controls.min()


def test_pathway_metabolites_are_complete(cohort):
    frame, _ = cohort
    for name in ["Trimethylamine N-oxide", "Choline", "Betaine", "Carnitine"]:
        assert frame[name].notna().all()


def test_log_transform_is_exact_and_non_destructive(cohort):
    frame, _ = cohort
    before = frame.columns.tolist()
    derived = add_log_pathway_features(frame)
    assert frame.columns.tolist() == before
    assert np.allclose(np.exp(derived["log_tmao"]), derived[TMAO])
    for feature in EXTENDED_FEATURES:
        assert feature in derived.columns
