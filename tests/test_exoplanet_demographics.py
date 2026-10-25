import json
from pathlib import Path

import matplotlib.image as mpimg
import numpy as np
import pandas as pd
import pytest

from pergamon import exoplanet_demographics as demographics


def test_demographics_notebook_uses_public_pipeline_and_documents_limits():
    repository = Path(__file__).parents[1]
    notebook_path = repository / "examples" / "exoplanet_demographics.ipynb"
    print(f"Reading from {notebook_path}...")
    notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
    source = "\n".join(
        "".join(cell["source"])
        for cell in notebook["cells"]
        if cell["cell_type"] in {"code", "markdown"}
    )

    assert "load_confirmed_exoplanets(limit=catalog_limit)" in source
    assert "plot_exoplanet_demographics(planets, visuals_directory" in source
    assert "not occurrence rates" in source
    assert "/Users/" not in json.dumps(notebook)

    readme_path = repository / "README.md"
    print(f"Reading from {readme_path}...")
    readme = readme_path.read_text(encoding="utf-8")
    assert "examples/exoplanet_demographics.ipynb" in readme


def test_derived_quantities_exclude_reported_limits():
    planets = pd.DataFrame({
        "pl_masse": [1.0, 8.0, np.nan], "pl_masselim": [0, 1, 0],
        "pl_rade": [1.0, 2.0, 1.5], "pl_radelim": [0, 0, 0],
    })

    result = demographics.add_demographic_quantities(planets)

    assert result.loc[0, "bulk_density_g_cm3"] == demographics.EARTH_DENSITY_G_CM3
    assert result.loc[0, "escape_velocity_km_s"] == demographics.EARTH_ESCAPE_VELOCITY_KM_S
    assert np.isnan(result.loc[1, "bulk_density_g_cm3"])
    assert np.isnan(result.loc[2, "escape_velocity_km_s"])


def test_repeated_demographic_plots_preserve_earlier_outputs(tmp_path, monkeypatch):
    planets = pd.DataFrame({"pl_masse": [1.0], "pl_rade": [1.0],
                            "pl_orbper": [365.25], "pl_insol": [1.0]})
    first = demographics.plot_exoplanet_demographics(planets, tmp_path)
    original = {path: path.read_bytes() for path in first.values()}
    second = demographics.plot_exoplanet_demographics(planets, tmp_path)
    assert set(first.values()).isdisjoint(second.values())
    assert all(path.read_bytes() == content for path, content in original.items())
    tag = "20261025_203000"
    collision = tmp_path / f"exoplanet_cosmic_shoreline_{tag}.png"
    collision.write_bytes(b"preserved earlier output")
    monkeypatch.setattr(demographics.plt, "subplots", lambda *args, **kwargs: pytest.fail(
        "An existing output must be rejected before creating figures"))
    before = {path: path.read_bytes() for path in tmp_path.iterdir()}
    with pytest.raises(FileExistsError):
        demographics.plot_exoplanet_demographics(planets, tmp_path, time_tag=tag)
    assert {path: path.read_bytes() for path in tmp_path.iterdir()} == before


@pytest.mark.parametrize("time_tag", ["../outside", "", "20261025", 123])
def test_demographic_time_tag_rejected(tmp_path, time_tag):
    with pytest.raises(ValueError, match="time_tag"):
        demographics.plot_exoplanet_demographics(pd.DataFrame(), tmp_path, time_tag=time_tag)


def test_demographics_cli_forwards_timestamp_and_saved_catalog(tmp_path, monkeypatch):
    import sys

    catalog = tmp_path / "planets.csv"
    planets = pd.DataFrame({"pl_masse": [1.0], "pl_rade": [1.0]})
    planets.to_csv(catalog, index=False)
    calls = []
    monkeypatch.setattr(sys, "argv", ["demographics", "--catalog", str(catalog),
                                     "--output-directory", str(tmp_path),
                                     "--time-tag", "20261025_203000", "--typefileplot", "pdf"])
    monkeypatch.setattr(demographics, "plot_exoplanet_demographics",
                        lambda *args, **kwargs: calls.append((args, kwargs)))
    monkeypatch.setattr(demographics, "load_confirmed_exoplanets",
                        lambda: pytest.fail("Saved-catalog replay must stay offline"))
    assert demographics.main() == 0
    pd.testing.assert_frame_equal(calls[0][0][0], planets)
    assert calls[0][0][1:] == (tmp_path, "pdf")
    assert calls[0][1] == {"time_tag": "20261025_203000", "typeplotback": "white"}


@pytest.mark.parametrize('typeplotback', ['white', 'dark'])
@pytest.mark.parametrize('inherit', [False, True])
def test_catalog_generates_all_three_demographic_figures(tmp_path, monkeypatch, typeplotback, inherit):
    from matplotlib.colors import to_rgba
    from tdpy.plotting import plot_background_colors, plot_style

    planets = pd.DataFrame({
        "pl_name": ["Earth twin", "Hot sub-Neptune", "Mass limit"],
        "pl_masse": [1.0, 8.0, 3.0], "pl_masselim": [0, 0, 1],
        "pl_rade": [1.0, 2.0, 1.4], "pl_radelim": [0, 0, 0],
        "pl_orbper": [365.25, 10.0, 3.0], "pl_orbperlim": [0, 0, 0],
        "pl_insol": [1.0, 100.0, 20.0],
    })

    original_planets = planets.copy(deep=True)
    background, foreground = plot_background_colors(typeplotback)
    save = demographics.plt.Figure.savefig

    def capture_and_save(figure, path, **kwargs):
        assert figure.get_facecolor() == to_rgba(background)
        assert all(axis.xaxis.label.get_color() == foreground for axis in figure.axes)
        axis = figure.axes[0]
        if path.stem.startswith("exoplanet_mass_radius_density_"):
            np.testing.assert_array_equal(axis.collections[0].get_offsets(), [[1, 1], [8, 2]])
            np.testing.assert_allclose(axis.collections[0].get_array(),
                                       [demographics.EARTH_DENSITY_G_CM3] * 2)
        elif path.stem.startswith("exoplanet_cosmic_shoreline_"):
            np.testing.assert_allclose(axis.collections[0].get_offsets(),
                                       [[1, demographics.EARTH_ESCAPE_VELOCITY_KM_S],
                                        [100, 2 * demographics.EARTH_ESCAPE_VELOCITY_KM_S]])
            guide = axis.lines[0]
            np.testing.assert_allclose(guide.get_ydata(),
                                       demographics.EARTH_ESCAPE_VELOCITY_KM_S * guide.get_xdata() ** 0.25)
            figure.canvas.draw()
            legend_box = axis.get_legend().get_window_extent(figure.canvas.get_renderer())
            assert legend_box.y0 > axis.get_window_extent().y1
        if path.stem.startswith(("exoplanet_mass_radius_density_", "exoplanet_cosmic_shoreline_")):
            offsets = axis.collections[0].get_offsets()
            assert np.all((offsets[:, 0] > axis.get_xlim()[0]) & (offsets[:, 0] < axis.get_xlim()[1]))
            assert np.all((offsets[:, 1] > axis.get_ylim()[0]) & (offsets[:, 1] < axis.get_ylim()[1]))
        return save(figure, path, **kwargs)

    monkeypatch.setattr(demographics.plt.Figure, 'savefig', capture_and_save)
    options = {} if inherit else {'typeplotback': typeplotback}
    outer_theme = typeplotback if inherit else ('dark' if typeplotback == 'white' else 'white')
    with plot_style(outer_theme):
        original_settings = dict(demographics.plt.rcParams)
        paths = demographics.plot_exoplanet_demographics(planets, tmp_path, **options)
        assert dict(demographics.plt.rcParams) == original_settings

    assert set(paths) == {
        "exoplanet_mass_radius_density", "exoplanet_cosmic_shoreline", "exoplanet_period_radius",
    }
    for path in paths.values():
        image = mpimg.imread(path)
        assert image.shape[0] > 100 and image.shape[1] > 100
        assert image[..., :3].min() < 0.8
        np.testing.assert_array_equal(image[0, 0], to_rgba(background))
    pd.testing.assert_frame_equal(planets, original_planets)