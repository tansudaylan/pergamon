import importlib.util
from pathlib import Path

import matplotlib.image as mpimg
import numpy as np
import pytest


EXAMPLE_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "examples"
    / "synthetic_vetting_performance.py"
)
SPECIFICATION = importlib.util.spec_from_file_location(
    "pergamon_synthetic_vetting_example", EXAMPLE_SCRIPT
)
example = importlib.util.module_from_spec(SPECIFICATION)
SPECIFICATION.loader.exec_module(example)


def make_result() -> dict:
    return {
        "typesyst": "PlanetarySystem",
        "listnameclastruetype": ["PlanetarySystem", "StellarBinary"],
        "indxclastruetypetarg": np.array([0, 0, 0, 1]),
        "listlablclasdisp": ["Low threshold", "High threshold"],
        "boolpositarg": [
            np.array([True, True, False, False]),
            np.array([True, False, False, False]),
        ],
    }


def test_synthetic_pipeline_routes_all_outputs_to_temporary_path(tmp_path, monkeypatch):
    captured_arguments = {}

    def fake_init(**arguments):
        captured_arguments.update(arguments)
        return make_result()

    monkeypatch.setattr(example.vetting, "run_candidate_vetting", fake_init)

    result = example.run_synthetic_pipeline(tmp_path)

    assert result["typesyst"] == "PlanetarySystem"
    assert captured_arguments["pathbase"] == str(tmp_path)
    assert captured_arguments["dicttroiinpt"]["pathbase"] == str(tmp_path)


@pytest.mark.parametrize('typeplotback', ['white', 'dark'])
@pytest.mark.parametrize('inherit', [False, True])
def test_synthetic_vetting_example_writes_nonblank_png(tmp_path, capsys, typeplotback, inherit):
    from tdpy.plotting import plot_style

    output_path = tmp_path / "synthetic_vetting_performance.png"

    options = {} if inherit else {'typeplotback': typeplotback}
    outer_theme = typeplotback if inherit else ('dark' if typeplotback == 'white' else 'white')
    with plot_style(outer_theme):
        original = dict(example.plt.rcParams)
        summary = example.run_example(output_path, result=make_result(), **options)
        assert dict(example.plt.rcParams) == original

    image = mpimg.imread(output_path)
    assert image.shape[0] > 100
    assert image.shape[1] > 100
    assert image[..., :3].min() < 0.8
    np.testing.assert_array_equal(image[0, 0], (0, 0, 0, 1) if typeplotback == 'dark' else (1, 1, 1, 1))
    np.testing.assert_allclose(summary["precision"], np.array([1.0, 1.0]))
    np.testing.assert_allclose(summary["recall"], np.array([2.0 / 3.0, 1.0 / 3.0]))
    assert f"Writing to {output_path}..." in capsys.readouterr().out


@pytest.mark.parametrize('typeplotback', ['white', 'dark'])
def test_synthetic_vetting_example_writes_bokeh_html(tmp_path, monkeypatch, typeplotback):
    pytest.importorskip("bokeh")
    from bokeh.models import Axis, Plot
    import tdpy.plotting as plotting

    save = plotting.save_bokeh_html

    def capture_and_save(figure, *args, **kwargs):
        result = save(figure, *args, **kwargs)
        background, foreground = plotting.plot_background_colors(typeplotback)
        assert all(plot.background_fill_color == background for plot in figure.select({'type': Plot}))
        assert all(axis.major_label_text_color == foreground for axis in figure.select({'type': Axis}))
        return result

    monkeypatch.setattr(plotting, 'save_bokeh_html', capture_and_save)
    output_path = tmp_path / "synthetic_vetting_performance.html"

    example.run_example(
        tmp_path / "synthetic_vetting_performance.png",
        result=make_result(),
        bokeh_output_path=output_path,
        typeplotback=typeplotback,
    )

    html = output_path.read_text(encoding="utf-8")
    assert "Input population" in html
    assert "Low-threshold confusion matrix" in html
    assert "Vetting performance" in html