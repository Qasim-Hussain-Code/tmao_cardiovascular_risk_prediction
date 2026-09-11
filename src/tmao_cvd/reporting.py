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
