"""Test configuration.

The package lives under src/ and the tests are run without installing it,
so the source directory is placed on the path here rather than in each test
module.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
