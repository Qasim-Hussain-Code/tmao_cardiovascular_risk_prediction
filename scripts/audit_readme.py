#!/usr/bin/env python3
"""Verify every figure in the README against results/metrics.

    python scripts/audit_readme.py

Reads every JSON file under results/metrics, extracts each literal figure, and
confirms it appears in the README at the same decimal precision. The point is
to make transcription errors impossible to commit: a number that drifts
between the pipeline and the prose fails the audit rather than sitting in the
document until someone happens to recompute it.

Matching is bounded on both sides, so 0.5 does not pass by matching inside
0.55, and 750 does not pass by matching inside 1750.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
METRICS_DIR = PROJECT_ROOT / "results" / "metrics"
README = PROJECT_ROOT / "README.md"

#: Keys whose values are prose or provenance rather than reportable figures.
#: Each is excluded deliberately, and the reason matters: requiring a run
#: timestamp or the full scope paragraph to appear verbatim in the README would
#: make the audit fail for reasons unconnected to numerical accuracy.
SKIP_KEYS = {
    "generated_utc",   # changes on every run
    "scope",           # prose, reproduced in the README in its own words
    "note",            # prose label
    "question",        # prose
    "source",          # prose
    "licence",         # quoted in the README as text, not as a figure
    "measurement_scale",
    "previously_reported_figure",
    "best_precursor",
    "accession",
}


def extract_figures(node, path: str = "") -> list[tuple[str, str]]:
    """Walk a JSON structure and return every reportable figure with its path."""

    found: list[tuple[str, str]] = []

    if isinstance(node, dict):
        for key, value in node.items():
            if key in SKIP_KEYS:
                continue
            found += extract_figures(value, f"{path}.{key}" if path else key)
    elif isinstance(node, list):
        for index, value in enumerate(node):
            found += extract_figures(value, f"{path}[{index}]")
    elif isinstance(node, bool):
        pass
    elif isinstance(node, (int, float)):
        found.append((path, repr(node)))
    elif isinstance(node, str) and re.fullmatch(r"[\d.]+e[+-]\d+", node):
        # Scientific notation is stored as a string to preserve its printed form.
        found.append((path, node))

    return found


def appears_in(text: str, figure: str) -> bool:
    """Is ``figure`` present, not embedded in a longer number?"""

    pattern = r"(?<![\d.])" + re.escape(figure) + r"(?![\d])"
    return re.search(pattern, text) is not None


def main() -> int:
    if not METRICS_DIR.exists():
        print(f"FAIL: {METRICS_DIR} does not exist")
        return 1

    readme = README.read_text()
    files = sorted(METRICS_DIR.glob("*.json"))
    if not files:
        print(f"FAIL: no JSON files under {METRICS_DIR}")
        return 1

    passed: list[str] = []
    failed: list[tuple[str, str, str]] = []

    for path in files:
        summary = json.loads(path.read_text())
        for location, figure in extract_figures(summary):
            if appears_in(readme, figure):
                passed.append(f"{path.name}:{location}")
            else:
                failed.append((path.name, location, figure))

    print(f"Audited {len(files)} file(s) under results/metrics against README.md")
    print(f"  figures checked : {len(passed) + len(failed)}")
    print(f"  passed          : {len(passed)}")
    print(f"  failed          : {len(failed)}")

    if failed:
        print("\nFigures absent from the README, or present at a different precision:")
        for name, location, figure in failed:
            print(f"  {figure:<14} {name}:{location}")
        print(
            "\nFix the README so each figure appears exactly as recorded. Do not "
            "adjust the recorded figures to match the prose."
        )
        return 1

    print("\nEvery recorded figure appears in the README at the recorded precision.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
