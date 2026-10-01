"""Synthetic catalog rows exercise enrichment without pretending to be observations."""

from datetime import datetime, timezone
import sys
from urllib.parse import parse_qs, urlparse

import numpy as np
import pandas as pd
import pytest

from pergamon.exoplanet_catalog import (
    add_measured_host_features,
    add_measured_orbit_orientation,
    attach_atmosphere_observations,
    build_planet_pair_catalog,
    build_system_catalog,
    enrich_exoplanets,
    load_confirmed_exoplanets,
    load_mast_spectra,
)


def example_planets():
    return pd.DataFrame({
        "pl_name": ["Test c", "Test b", "Test d", "Other b"],
        "hostname": ["Test", "Test", "Test", "Other"],
        "ra": [266.4051, 266.4051, 266.4051, np.nan],
        "dec": [-28.936175, -28.936175, -28.936175, np.nan],
        "pl_orbper": [2.0, 1.0, 4.0, 10.0],  # [day]
        "pl_masse": [2.0, 1.0, np.nan, np.nan],  # [Earth mass]
        "pl_rade": [1.2, 1.0, 1.5, 1.0],  # [Earth radius]
        "st_mass": [1.0, 1.0, 1.0, 1.0],  # [solar mass]
        "st_rad": [1.0, 1.0, 1.0, 1.0],  # [solar radius]
        "st_rotp": [20.0, 20.0, 20.0, np.nan],  # [day]
        "pl_orbeccen": [0.1, 0.0, np.nan, np.nan],
        "sy_pnum": [3, 3, 3, 1],
    })


def test_archive_query_uses_default_solutions_and_limits(monkeypatch):
    captured = []
    monkeypatch.setattr(pd, "read_csv", lambda url: captured.append(url) or pd.DataFrame())
    load_confirmed_exoplanets(limit=2)
    query = parse_qs(urlparse(captured[0]).query)["query"][0]
    assert "top 2" in query and "from ps where default_flag=1" in query
    assert "pl_masse" in query and "pl_bmassprov" in query
    with pytest.raises(ValueError, match="positive integer"):
        load_confirmed_exoplanets(limit=0)


def test_orbital_features_require_true_mass_and_preserve_missing_values():
    catalog = enrich_exoplanets(example_planets())
    first = catalog.set_index("pl_name").loc["Test c"]
    assert abs(first["galactic_longitude_deg"]) < 1.0  # [deg]
    assert first["hill_radius_au"] > 0
    assert first["safronov_number"] > 0
    assert first["eccentricity_amd_lower_bound_kg_m2_s"] > 0
    assert np.isnan(catalog.set_index("pl_name").loc["Test d", "hill_radius_au"])
    assert np.isnan(catalog.set_index("pl_name").loc["Other b", "galactic_latitude_deg"])
    assert bool(catalog.set_index("pl_name").loc["Test c", "has_inner_companion"])
    assert bool(catalog.set_index("pl_name").loc["Test c", "has_outer_companion"])
    assert not bool(catalog.set_index("pl_name").loc["Test b", "has_inner_companion"])


def test_reported_mass_limit_cannot_generate_derived_hill_radii():
    planets = example_planets().iloc[:2].copy()
    planets["pl_masselim"] = [1, 0]
    assert np.isnan(enrich_exoplanets(planets).loc[planets.index[0], "hill_radius_au"])
    assert np.isnan(build_planet_pair_catalog(planets).loc[0, "mutual_hill_separation"])


def test_eccentricity_bound_is_not_a_measured_angular_momentum_deficit():
    planets = example_planets().iloc[:1].copy()
    planets["pl_orbeccenlim"] = [1]
    enriched = enrich_exoplanets(planets)
    assert np.isnan(enriched["eccentricity_amd_lower_bound_kg_m2_s"].iloc[0])
    oriented = add_measured_orbit_orientation(enriched, pd.DataFrame({
        "pl_name": ["Test c"], "inclination_invariable_plane_deg": [0.0],
        "source_reference": ["orbital fit"],
    }))
    assert np.isnan(oriented["angular_momentum_deficit_kg_m2_s"].iloc[0])


def test_partial_system_does_not_assert_absent_companions():
    catalog = enrich_exoplanets(example_planets().iloc[:2])
    assert pd.isna(catalog.set_index("pl_name").loc["Test b", "has_inner_companion"])
    assert bool(catalog.set_index("pl_name").loc["Test b", "has_outer_companion"])
    assert not bool(build_system_catalog(catalog).loc[0, "system_complete"])


def test_pairs_preserve_all_members_and_measure_only_supported_hill_separations():
    catalog = example_planets()
    systems = build_system_catalog(catalog)
    assert systems.set_index("hostname").loc["Test", "downloaded_planet_count"] == 3
    pairs = build_planet_pair_catalog(catalog)
    test_pairs = pairs[pairs["hostname"] == "Test"]
    assert len(test_pairs) == 3
    assert len(build_planet_pair_catalog(catalog, neighbors_only=True)) == 2
    assert test_pairs.set_index(["inner_planet", "outer_planet"]).loc[("Test b", "Test c"),
                                                                      "mutual_hill_separation"] > 0
    assert np.isnan(test_pairs.set_index(["inner_planet", "outer_planet"]).loc[("Test c", "Test d"),
                                                                               "mutual_hill_separation"])
    unordered = catalog.iloc[:2].copy()
    unordered.loc[unordered.index[0], "pl_orbper"] = np.nan
    pair = build_planet_pair_catalog(unordered).iloc[0]
    assert pd.isna(pair["inner_planet"]) and pd.isna(pair["outer_planet"])
    assert {pair["planet_a"], pair["planet_b"]} == {"Test b", "Test c"}


def test_measured_features_require_provenance_and_valid_inputs():
    planets = enrich_exoplanets(example_planets())
    hosts = pd.DataFrame({"hostname": ["Test"], "source_reference": ["measured source"],
                          "host_xuv_luminosity_erg_s": [1e28],
                          "host_convective_turnover_days": [10.0],
                          "host_alfven_radius_stellar_radii": [10.0],
                          "host_flare_rate_per_day": [0.2]})
    joined = add_measured_host_features(planets, hosts)
    first = joined.set_index("pl_name").loc["Test c"]
    assert first["host_rossby_number"] == 2.0
    assert first["xuv_flux_erg_s_cm2"] > 0
    assert first["host_alfven_radius_au"] > 0
    assert first["host_source_reference"] == "measured source"
    oriented = add_measured_orbit_orientation(joined, pd.DataFrame({
        "pl_name": ["Test c"], "inclination_invariable_plane_deg": [0.0],
        "source_reference": ["orbital fit"],
    }))
    first = oriented.set_index("pl_name").loc["Test c"]
    assert first["orientation_source_reference"] == "orbital fit"
    assert first["angular_momentum_deficit_kg_m2_s"] == pytest.approx(
        first["eccentricity_amd_lower_bound_kg_m2_s"]
    )
    assert pd.isna(oriented.set_index("pl_name").loc["Test b", "angular_momentum_deficit_kg_m2_s"])
    with pytest.raises(ValueError, match="unique"):
        add_measured_host_features(planets, pd.concat([hosts, hosts], ignore_index=True))


def test_atmosphere_observations_keep_multiple_and_unmatched_sources(monkeypatch):
    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def read(self):
            return b'{"filenames":["spectrum1.txt","spectrum2.txt"]}'

    import pergamon.exoplanet_catalog as module
    monkeypatch.setattr(module, "urlopen", lambda url, timeout: Response())
    spectra = load_mast_spectra("Test b")
    assert len(spectra) == 2 and spectra["source_url"].str.contains("/file/").all()
    extra = pd.DataFrame({"pl_name": ["Unknown b"], "source": ["EAOT"],
                          "observation_id": ["eaot-1"], "source_url": ["https://example.org/observation"]})
    enriched, observations = attach_atmosphere_observations(
        example_planets(), pd.concat([spectra, extra], ignore_index=True)
    )
    assert enriched.set_index("pl_name").loc["Test b", "atmosphere_observation_count"] == 2
    assert len(observations) == 3 and not observations.loc[2, "matched_archive_planet"]


def test_snapshot_command_retains_previous_outputs(monkeypatch, tmp_path):
    import pergamon.exoplanet_catalog as module

    class FixedTime:
        @staticmethod
        def now(zone):
            return datetime(2026, 9, 30, tzinfo=timezone.utc)

    monkeypatch.setattr(module, "datetime", FixedTime)
    monkeypatch.setattr(module, "load_confirmed_exoplanets", lambda limit: example_planets())
    monkeypatch.setattr(sys, "argv", ["exoplanet_catalog", "--output-directory", str(tmp_path)])
    module.main()
    snapshots = sorted(tmp_path.glob("exoplanet_*.csv"))
    assert len(snapshots) == 4
    planets = tmp_path / "exoplanet_planets_20260930_000000.csv"
    print(f"Reading from {planets}...")
    assert len(pd.read_csv(planets)) == len(example_planets())
    with pytest.raises(FileExistsError, match="already exists"):
        module.main()
    assert sorted(tmp_path.glob("exoplanet_*.csv")) == snapshots