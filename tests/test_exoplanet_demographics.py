import matplotlib.image as mpimg
import numpy as np
import pandas as pd

from pergamon import exoplanet_demographics as demographics


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


def test_catalog_generates_all_three_demographic_figures(tmp_path):
    planets = pd.DataFrame({
        "pl_name": ["Earth twin", "Hot sub-Neptune", "Mass limit"],
        "pl_masse": [1.0, 8.0, 3.0], "pl_masselim": [0, 0, 1],
        "pl_rade": [1.0, 2.0, 1.4], "pl_radelim": [0, 0, 0],
        "pl_orbper": [365.25, 10.0, 3.0], "pl_orbperlim": [0, 0, 0],
        "pl_insol": [1.0, 100.0, 20.0],
    })

    paths = demographics.plot_exoplanet_demographics(planets, tmp_path)

    assert set(paths) == {
        "exoplanet_mass_radius_density", "exoplanet_cosmic_shoreline", "exoplanet_period_radius",
    }
    for path in paths.values():
        image = mpimg.imread(path)
        assert image.shape[0] > 100 and image.shape[1] > 100
        assert image[..., :3].min() < 0.8