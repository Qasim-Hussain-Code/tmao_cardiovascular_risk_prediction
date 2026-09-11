"""Reporting rules that the plan requires to be enforced, not remembered.

Two rules from docs/reanalysis_plan.md are implemented here rather than left
to whoever writes the manuscript.

Section 7 requires that the claim scope travels with every reported result. A
collapsed question 2 verdict beside an unscoped headline is the thing that gets
screenshotted out of a repository and circulated without its qualifications, so
the renderer takes the scope as an argument and verifies that no element of it
has been dropped.

Question 3 of section 6 requires that a null diagnostic result is never emitted
without the section 5 limitation attached, in the same sentence. A null there
is the result most likely to be misread as exoneration. The null verdict string
is therefore constructed with the caveat embedded, and there is no code path
that produces one without the other.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

#: The claim scope from section 7, printed at the head of every report.
SCOPE_STATEMENT = (
    "SCOPE OF EVERY RESULT BELOW: these findings come from a single cohort of "
    "post-PCI secondary prevention patients with stable angina, of Asian ethnicity, "
    "uniformly on dual antiplatelet therapy, with a symptom-driven endpoint "
    "(recurrent angina at nine months). Metabolite values are relative peak areas in "
    "arbitrary instrument units, not concentrations, so no published cut point "
    "applies. There is no external validation, because the two other cohorts in the "
    "source publication were never deposited. Acquisition order is perfectly "
    "confounded with outcome and no batch metadata exists. Nothing here transfers to "
    "primary prevention, to hard cardiovascular endpoints, or to other populations."
)

#: Every element section 7 requires. Checked against the scope text so that
#: quietly trimming the statement fails loudly instead of silently.
REQUIRED_SCOPE_PHRASES = (
    "secondary prevention",
    "asian",
    "dual antiplatelet",
    "symptom-driven",
    "peak areas",
    "single cohort",
    "no external validation",
    "acquisition order",
)

#: The section 5 limitation, in the words that must accompany a null question 3.
SECTION_5_CAVEAT = (
    "this does not rule the confound out, because it was tested under the "
    "unverifiable assumption that the deposited order is the acquisition order, and "
    "the deposit contains no batch identifier, run order or acquisition date that "
    "could establish whether it is"
)


class ReportingError(RuntimeError):
    """Raised when a reporting rule from the plan would be violated."""


@dataclass(frozen=True)
class QuestionResult:
    """One pre-specified question, its verdict, and the table behind it."""

    question: str
    verdict: str
    statement: str
    table: pd.DataFrame


def verify_scope(scope: str) -> None:
    """Raise unless every element of the section 7 scope is present."""

    lowered = scope.lower()
    missing = [phrase for phrase in REQUIRED_SCOPE_PHRASES if phrase not in lowered]
    if missing:
        raise ReportingError(
            "the scope statement is missing required elements "
            f"{missing}. Section 7 of docs/reanalysis_plan.md requires the full claim "
            "scope to accompany every reported result. Restore the missing elements "
            "rather than removing this check."
        )


def order_structure_statement(detected: bool, detail: str) -> str:
    """Build the question 3 verdict sentence.

    The null branch embeds the section 5 caveat in the same sentence as the
    result. This is the only constructor for that sentence in the codebase, so
    a bare null cannot be emitted by any path.
    """

    if detected:
        return (
            f"Acquisition order structure IS detectable along the deposited ordering "
            f"({detail}). Any between group difference in this dataset is therefore "
            "confounded with it to an unknown degree, and every other result in this "
            "analysis is reported under that caveat."
        )

    return (
        f"No acquisition order structure was detected along the deposited ordering "
        f"({detail}), but {SECTION_5_CAVEAT}."
    )


#: Emitted against every other question when question 3 detects order
#: structure. Section 6 of the plan requires it ("every other result in this
#: analysis is then reported under that caveat"), so the renderer attaches it
#: rather than relying on whoever writes the manuscript to remember.
ORDER_STRUCTURE_CAVEAT = (
    "!! READ WITH QUESTION 3: acquisition order structure was detected in this "
    "dataset. This result is confounded with it to an unknown degree and cannot be "
    "interpreted as biological."
)

ORDER_STRUCTURE_VERDICT = "ORDER STRUCTURE DETECTED"


@dataclass(frozen=True)
class PanelContext:
    """Where the marker under test sits among all metabolites in the panel.

    Every field is required. The class exists so that the question 1 statement
    cannot be constructed without this context, rather than relying on the
    caller to remember to append it.
    """

    marker: str
    univariate_auc: float
    percentile: float
    ranking_higher: int
    panel_size: int
    best_precursor: str
    best_precursor_auc: float
    panel_median_auc: float


#: The context that must accompany the question 1 difference in area wherever
#: it appears. An incremental contribution measured inside a panel that is
#: globally displaced is not evidence about the marker, and the bare figure
#: invites exactly that reading.
Q1_CONTEXT_TEMPLATE = (
    "This figure must be read with the panel context: {marker} has a univariate "
    "area under the curve of {auc:.4f}, placing it at the {percentile:.1f}th "
    "percentile of the {panel_size} metabolites in this deposit, with {higher} "
    "ranking above it, among them {best_precursor}, one of its own precursors, at "
    "{best_precursor_auc:.4f}. The panel median is {median:.4f} against a null of "
    "0.5, so an incremental contribution measured inside it is not evidence about "
    "the biology of the marker."
)


def incremental_value_statement(
    verdict: str,
    delta_auc: float,
    delong_p: float,
    separation: float,
    context: PanelContext,
) -> str:
    """Build the question 1 verdict sentence, with panel context attached.

    The only constructor for any sentence containing the question 1 difference
    in area. ``context`` is required and unconditional, so there is no code
    path that emits the difference without it. This mirrors
    :func:`order_structure_statement`, and exists for the same reason: the
    number is more persuasive on its own than the evidence supports.
    """

    if verdict == "POSITIVE":
        head = (
            f"In this cohort the marker carries information about the outcome beyond "
            f"its dietary precursors. The difference in out of fold area under the "
            f"curve is {delta_auc:+.4f} (DeLong p = {delong_p:.3g}), and the extended "
            f"model's decision curve lies above the baseline model's across "
            f"{separation:.0%} of the pre-specified threshold range."
        )
    elif verdict == "NEGATIVE":
        head = (
            f"The marker adds nothing beyond its precursors here. The difference in "
            f"out of fold area under the curve is {delta_auc:+.4f} with DeLong p = "
            f"{delong_p:.3g}, which does not exclude zero."
        )
    else:
        head = (
            f"A statistically detectable but clinically negligible contribution. The "
            f"difference in area is {delta_auc:+.4f} (DeLong p = {delong_p:.3g}), but "
            f"the decision curves separate across only {separation:.0%} of the "
            f"threshold range. Per the plan the decision curve governs the wording."
        )

    return head + " " + Q1_CONTEXT_TEMPLATE.format(
        marker=context.marker,
        auc=context.univariate_auc,
        percentile=context.percentile,
        panel_size=context.panel_size,
        higher=context.ranking_higher,
        best_precursor=context.best_precursor,
        best_precursor_auc=context.best_precursor_auc,
        median=context.panel_median_auc,
    )


def render_report(
    results: list[QuestionResult], scope: str = SCOPE_STATEMENT
) -> str:
    """Render the verdicts beneath the claim scope.

    The scope is verified before anything else is written, so a report cannot
    be produced without it. If question 3 detected acquisition order structure,
    every other result carries the caveat section 6 requires.
    """

    verify_scope(scope)
    order_structure = any(r.verdict == ORDER_STRUCTURE_VERDICT for r in results)

    lines: list[str] = ["=" * 78, "REANALYSIS OF ST001420", "=" * 78, ""]
    for sentence in _wrap(scope, 78):
        lines.append(sentence)
    lines += ["", "=" * 78, ""]

    for result in results:
        lines.append(f"{result.question}")
        lines.append("-" * 78)
        lines.append(f"VERDICT: {result.verdict}")
        if order_structure and result.verdict != ORDER_STRUCTURE_VERDICT:
            lines += _wrap(ORDER_STRUCTURE_CAVEAT, 78)
        lines += _wrap(result.statement, 78)
        lines.append("")
        lines.append(result.table.to_string(index=False))
        lines.append("")

    return "\n".join(lines)


def _wrap(text: str, width: int) -> list[str]:
    """Simple greedy wrap, kept local to avoid a dependency for six lines."""

    words = text.split()
    lines: list[str] = []
    current: list[str] = []
    for word in words:
        if sum(len(w) + 1 for w in current) + len(word) > width:
            lines.append(" ".join(current))
            current = [word]
        else:
            current.append(word)
    if current:
        lines.append(" ".join(current))
    return lines
