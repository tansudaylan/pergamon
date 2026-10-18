import importlib.util
from pathlib import Path

import numpy as np
import pytest

EXAMPLE_SCRIPT = (
    Path(__file__).parents[1] / "examples" / "compact_object_signatures.py"
)
SPEC = importlib.util.spec_from_file_location("compact_object_signatures", EXAMPLE_SCRIPT)
example = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(example)


def test_compact_object_signatures_write_bokeh_html(tmp_path):
    pytest.importorskip("bokeh")
    html_path = tmp_path / "compact_object_signatures.html"

    signatures = example.run_example(
        tmp_path / "compact_object_signatures.png",
        bokeh_output_path=html_path,
    )

    assert set(signatures) == {"beaming", "ellipsoidal", "self_lensing"}
    html = html_path.read_text(encoding="utf-8")
    assert all(title in html for title in ("Doppler beaming", "Ellipsoidal variation", "Self-lensing"))