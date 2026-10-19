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


@pytest.mark.parametrize('typeplotback', ['white', 'dark'])
def test_compact_object_signatures_write_bokeh_html(tmp_path, monkeypatch, typeplotback):
    pytest.importorskip("bokeh")
    from bokeh.models import Axis, Plot
    import tdpy.plotting as plotting

    html_path = tmp_path / "compact_object_signatures.html"
    save = plotting.save_bokeh_html

    def capture_and_save(figure, *args, **kwargs):
        result = save(figure, *args, **kwargs)
        background, foreground = plotting.plot_background_colors(typeplotback)
        assert all(plot.background_fill_color == background for plot in figure.select({'type': Plot}))
        assert all(axis.major_label_text_color == foreground for axis in figure.select({'type': Axis}))
        return result

    monkeypatch.setattr(plotting, 'save_bokeh_html', capture_and_save)

    signatures = example.run_example(
        tmp_path / "compact_object_signatures.png",
        bokeh_output_path=html_path,
        typeplotback=typeplotback,
    )

    assert set(signatures) == {"beaming", "ellipsoidal", "self_lensing"}
    html = html_path.read_text(encoding="utf-8")
    assert all(title in html for title in ("Doppler beaming", "Ellipsoidal variation", "Self-lensing"))