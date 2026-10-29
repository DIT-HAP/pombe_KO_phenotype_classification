"""pytest configuration — add project src/ and scripts/ to sys.path."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
for folder in (ROOT / "scripts", ROOT / "src"):
    if str(folder) not in sys.path:
        sys.path.insert(0, str(folder))
