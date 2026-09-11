"""Tests for the reporting rules the plan requires to be enforced.

Both rules tested here exist because a document cannot enforce itself. Section
7 requires the claim scope to travel with every result, and question 3 of
section 6 requires a null diagnostic never to be emitted without the section 5
limitation attached.
"""

from __future__ import annotations

import pandas as pd
import pytest

from tmao_cvd.reporting import (
    REQUIRED_SCOPE_PHRASES,
    PanelContext,
    incremental_value_statement,
    SCOPE_STATEMENT,
    SECTION_5_CAVEAT,
    QuestionResult,
    ReportingError,
    order_structure_statement,
    render_report,
    verify_scope,
)


def test_the_shipped_scope_statement_is_complete():
    verify_scope(SCOPE_STATEMENT)


@pytest.mark.parametrize("phrase", REQUIRED_SCOPE_PHRASES)
def test_dropping_any_single_scope_element_raises(phrase):
    # Trimming the scope must fail loudly. Every element is load bearing.
    trimmed = SCOPE_STATEMENT.lower().replace(phrase, "")
    with pytest.raises(ReportingError, match="missing required elements"):
        verify_scope(trimmed)


def test_a_report_cannot_be_rendered_without_a_valid_scope():
    result = QuestionResult("Q", "VERDICT", "statement", pd.DataFrame([{"a": 1}]))
    with pytest.raises(ReportingError):
        render_report([result], scope="results from a cohort")


def test_rendered_report_leads_with_the_scope():
    result = QuestionResult("Q1", "NEGATIVE", "nothing found", pd.DataFrame([{"a": 1}]))
    text = render_report([result])
    scope_position = text.index("SCOPE OF EVERY RESULT BELOW")
    verdict_position = text.index("VERDICT:")
    assert scope_position < verdict_position


def test_null_order_structure_result_carries_the_section_5_caveat():
    # The central rule: a clean null must never stand alone, because it reads
    # as exoneration and the confound was not in fact ruled out.
    statement = order_structure_statement(False, "cases: 4.8% at p < 0.05")
    assert SECTION_5_CAVEAT in statement
    assert "does not rule the confound out" in statement


def test_caveat_and_result_appear_in_the_same_sentence():
    # Section 6 requires same sentence, equal prominence, never a footnote.
    statement = order_structure_statement(False, "detail here")
    sentences = [s for s in statement.split(". ") if s.strip()]
    carrying = [s for s in sentences if "No acquisition order structure" in s]
    assert len(carrying) == 1
    assert "does not rule the confound out" in carrying[0]


def test_detected_branch_states_the_consequence_for_other_results():
    statement = order_structure_statement(True, "detail here")
    assert "confounded" in statement
    assert "every other result" in statement.lower()


def test_there_is_only_one_constructor_for_the_q3_statement():
    # If a second code path starts building this sentence, the caveat coupling
    # can be bypassed. Keep construction funnelled through one function.
    import inspect

    from tmao_cvd import reanalysis

    source = inspect.getsource(reanalysis)
    assert "No acquisition order structure" not in source
    assert "order_structure_statement(" in source


# --- Question 1 context guard, mirroring the question 3 caveat guard ---------

CONTEXT = PanelContext(
    marker="Trimethylamine N-oxide",
    univariate_auc=0.7404,
    percentile=64.2,
    ranking_higher=144,
    panel_size=405,
    best_precursor="Betaine",
    best_precursor_auc=0.8427,
    panel_median_auc=0.6852,
)


@pytest.mark.parametrize("verdict", ["POSITIVE", "NEGATIVE", "DISCORDANT"])
def test_every_q1_verdict_carries_the_panel_context(verdict):
    # The difference in area is persuasive on its own beyond what the evidence
    # supports, so no verdict may emit it bare.
    statement = incremental_value_statement(verdict, 0.0256, 0.00547, 0.9821, CONTEXT)
    assert "percentile" in statement
    assert "144" in statement
    assert "Betaine" in statement
    assert "not evidence about" in statement


def test_the_difference_and_its_context_are_inseparable():
    statement = incremental_value_statement("POSITIVE", 0.0256, 0.00547, 0.9821, CONTEXT)
    assert "+0.0256" in statement
    # The context must follow the figure in the same emitted string, not be a
    # separate paragraph a reader may not reach.
    assert statement.index("+0.0256") < statement.index("percentile")


def test_panel_context_has_no_optional_fields():
    # A default would let a caller construct an empty context and satisfy the
    # signature while defeating the guard.
    import dataclasses

    for field in dataclasses.fields(PanelContext):
        assert field.default is dataclasses.MISSING
        assert field.default_factory is dataclasses.MISSING


def test_no_second_code_path_formats_the_q1_difference():
    # The structural guarantee: the analysis module must not build any string
    # containing the difference in area. Only the reporting constructor may.
    import inspect

    from tmao_cvd import reanalysis

    source = inspect.getsource(reanalysis)
    assert "{difference:" not in source
    assert "difference in out of fold area" not in source
    assert "incremental_value_statement(" in source


def test_rendered_report_never_shows_the_difference_without_context():
    statement = incremental_value_statement("POSITIVE", 0.0256, 0.00547, 0.9821, CONTEXT)
    result = QuestionResult("QUESTION 1", "POSITIVE", statement, pd.DataFrame([{"a": 1}]))
    text = render_report([result])
    assert "+0.0256" in text
    assert "percentile" in text
