import importlib.util
from pathlib import Path

import matplotlib.image as mpimg
import numpy as np
import pytest


EXAMPLE_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "examples"
    / "simulated_stellar_populations.py"
)
SPECIFICATION = importlib.util.spec_from_file_location(
    "pergamon_simulated_stellar_populations", EXAMPLE_SCRIPT
)
example = importlib.util.module_from_spec(SPECIFICATION)
SPECIFICATION.loader.exec_module(example)


@pytest.mark.parametrize('typeplotback', ['white', 'dark'])
@pytest.mark.parametrize('inherit', [False, True])
def test_simulated_stellar_populations_example_runs_pipeline(tmp_path, capsys, typeplotback, inherit):
    from tdpy.plotting import plot_style

    output_path = tmp_path / "simulated_stellar_populations.png"

    options = {} if inherit else {'typeplotback': typeplotback}
    outer_theme = typeplotback if inherit else ('dark' if typeplotback == 'white' else 'white')
    with plot_style(outer_theme):
        original = dict(example.plt.rcParams)
        result = example.run_example(output_path, **options)
        assert dict(example.plt.rcParams) == original

    image = mpimg.imread(output_path)
    assert image.shape[0] > 100
    assert image.shape[1] > 100
    assert image[..., :3].min() < 0.8
    np.testing.assert_array_equal(image[0, 0], (0, 0, 0, 1) if typeplotback == 'dark' else (1, 1, 1, 1))
    np.testing.assert_allclose(
        result["reference"]["massstar"][0],
        np.linspace(0.75, 1.15, 40),
    )
    output = capsys.readouterr().out
    assert "Reading from" in output
    assert f"Writing to {output_path}..." in output


@pytest.mark.parametrize('typeplotback', ['white', 'dark'])
def test_simulated_stellar_populations_bokeh_variant_writes_html(tmp_path, monkeypatch, typeplotback):
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
    output_path = example.run_bokeh_example(tmp_path / "stellar_populations.html", typeplotback=typeplotback)

    assert output_path.is_file()
    html = output_path.read_text(encoding="utf-8")
    assert "Reference" in html
    assert "Radius-enhanced" in html
