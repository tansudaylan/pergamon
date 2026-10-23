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