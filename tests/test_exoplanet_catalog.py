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
    monkeypatch.setattr(pd, "read_csv", lambda url, **kwargs: captured.append(url) or pd.DataFrame())
    load_confirmed_exoplanets(limit=2)
    query = parse_qs(urlparse(captured[0]).query)["query"][0]
    assert "top 2" in query and "from ps where default_flag=1" in query
    assert "pl_masse" in query and "pl_bmassprov" in query
    assert "st_radlim" in query and "st_rotplim" in query
    with pytest.raises(ValueError, match="positive integer"):
        load_confirmed_exoplanets(limit=0)


def test_archive_reader_preserves_large_gaia_identifiers_with_missing_rows(monkeypatch):
    from io import StringIO

    identifier = '5853498713190525696'
    original_read_csv = pd.read_csv

    def read_csv(url, **kwargs):
        assert kwargs['dtype'] == {'gaia_dr3_id': 'string'}
        return original_read_csv(StringIO(f'pl_name,gaia_dr3_id\nTest b,{identifier}\nOther b,\n'), **kwargs)

    monkeypatch.setattr(pd, 'read_csv', read_csv)
    result = load_confirmed_exoplanets()
    assert result['gaia_dr3_id'].iloc[0] == identifier
    assert pd.isna(result['gaia_dr3_id'].iloc[1])


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


@pytest.mark.parametrize('eccentricity, inclination', [(1e-10, 0.0), (0.0, 1e-8), (1e-10, 1e-8)])
def test_tiny_orbital_deficits_remain_positive(eccentricity, inclination):
    from astropy import constants, units

    planets = example_planets().iloc[:1].copy()
    planets['pl_orbeccen'] = eccentricity
    enriched = enrich_exoplanets(planets)
    oriented = add_measured_orbit_orientation(enriched, pd.DataFrame({
        'pl_name': ['Test c'], 'inclination_invariable_plane_deg': [inclination],
        'source_reference': ['synthetic invariant-plane fit'],
    }))
    semimajor_axis = enriched['kepler_semimajor_axis_au'].iloc[0] * units.au
    momentum = (2 * units.M_earth * np.sqrt(constants.G * units.M_sun * semimajor_axis)).to_value(
        units.kg * units.m**2 / units.s)
    small_angle_fraction = (eccentricity**2 + np.deg2rad(inclination)**2) / 2
    actual = oriented['angular_momentum_deficit_kg_m2_s'].iloc[0]
    assert actual > 0
    assert actual == pytest.approx(momentum * small_angle_fraction, rel=1e-12)
    assert enriched['eccentricity_amd_lower_bound_kg_m2_s'].iloc[0] == pytest.approx(
        momentum * eccentricity**2 / 2, rel=1e-12, abs=0)


@pytest.mark.parametrize('eccentricity, inclination', [(0.0, 0.0), (0.1, 20.0), (0.8, 90.0), (0.9, 180.0)])
def test_orbital_deficit_fraction_preserves_ordinary_formula(eccentricity, inclination):
    from pergamon.exoplanet_catalog import _angular_momentum_deficit_fraction

    actual = _angular_momentum_deficit_fraction(np.array([eccentricity]), inclination)[0]
    expected = 1 - np.sqrt(1 - eccentricity**2) * np.cos(np.deg2rad(inclination))
    assert actual == pytest.approx(expected, rel=1e-13, abs=1e-15)


def test_partial_system_does_not_assert_absent_companions():
    catalog = enrich_exoplanets(example_planets().iloc[:2])
    assert pd.isna(catalog.set_index("pl_name").loc["Test b", "has_inner_companion"])
    assert bool(catalog.set_index("pl_name").loc["Test b", "has_outer_companion"])
    assert not bool(build_system_catalog(catalog).loc[0, "system_complete"])


@pytest.mark.parametrize('invalid_period', [np.inf, -np.inf, np.nan])
def test_unknown_periods_do_not_create_companion_order_or_complete_neighbors(invalid_period):
    planets = example_planets().iloc[:3].copy()
    planets.loc[planets['pl_name'] == 'Test d', 'pl_orbper'] = invalid_period
    before = planets.copy(deep=True)

    enriched = enrich_exoplanets(planets).set_index('pl_name')
    pairs = build_planet_pair_catalog(planets)

    assert pd.isna(enriched.loc['Test b', 'has_inner_companion'])
    assert bool(enriched.loc['Test b', 'has_outer_companion'])
    assert bool(enriched.loc['Test c', 'has_inner_companion'])
    assert pd.isna(enriched.loc['Test c', 'has_outer_companion'])
    assert pd.isna(enriched.loc['Test d', 'has_inner_companion'])
    assert pd.isna(enriched.loc['Test d', 'has_outer_companion'])
    assert len(pairs) == 3
    assert pairs['are_neighbors'].isna().all()
    assert build_planet_pair_catalog(planets, neighbors_only=True).empty
    unknown_pairs = pairs[(pairs['planet_a'] == 'Test d') | (pairs['planet_b'] == 'Test d')]
    assert len(unknown_pairs) == 2
    assert unknown_pairs[['inner_planet', 'outer_planet', 'period_ratio', 'mutual_hill_separation']].isna().all().all()
    finite_pair = pairs[(pairs['inner_planet'] == 'Test b') & (pairs['outer_planet'] == 'Test c')].iloc[0]
    assert finite_pair['period_ratio'] == 2.0
    assert finite_pair['mutual_hill_separation'] > 0.0
    pd.testing.assert_frame_equal(planets, before)
    pd.testing.assert_series_equal(enriched['pl_orbper'].reset_index(drop=True), before['pl_orbper'].reset_index(drop=True))


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


@pytest.mark.parametrize('reference', ['', ' ', '\t\n'])
@pytest.mark.parametrize('measurement_kind', ['host', 'orientation'])
def test_measured_features_reject_blank_provenance(reference, measurement_kind):
    planets = enrich_exoplanets(example_planets())
    if measurement_kind == 'host':
        measurements = pd.DataFrame({
            'hostname': ['Test'], 'source_reference': [reference],
            'host_convective_turnover_days': [10.0],  # [day]
        })
        attach = add_measured_host_features
    else:
        measurements = pd.DataFrame({
            'pl_name': ['Test c'], 'source_reference': [reference],
            'inclination_invariable_plane_deg': [0.0],  # [deg]
        })
        attach = add_measured_orbit_orientation
    original_planets = planets.copy(deep=True)
    original_measurements = measurements.copy(deep=True)
    with pytest.raises(ValueError, match='source reference|provenance'):
        attach(planets, measurements)
    pd.testing.assert_frame_equal(planets, original_planets)
    pd.testing.assert_frame_equal(measurements, original_measurements)


def test_measured_features_require_provenance_and_valid_inputs():
    planets = enrich_exoplanets(example_planets())
    hosts = pd.DataFrame({"hostname": ["Test"], "source_reference": [" measured source "],
                          "host_xuv_luminosity_erg_s": [1e28],
                          "host_convective_turnover_days": [10.0],
                          "host_alfven_radius_stellar_radii": [10.0],
                          "host_flare_rate_per_day": [0.2]})
    joined = add_measured_host_features(planets, hosts)
    first = joined.set_index("pl_name").loc["Test c"]
    assert first["host_rossby_number"] == 2.0
    assert first["xuv_flux_erg_s_cm2"] > 0
    assert first["host_alfven_radius_au"] > 0
    assert first["host_source_reference"] == " measured source "
    oriented = add_measured_orbit_orientation(joined, pd.DataFrame({
        "pl_name": ["Test c"], "inclination_invariable_plane_deg": [0.0],
        "source_reference": [" orbital fit "],
    }))
    first = oriented.set_index("pl_name").loc["Test c"]
    assert first["orientation_source_reference"] == " orbital fit "
    assert first["angular_momentum_deficit_kg_m2_s"] == pytest.approx(
        first["eccentricity_amd_lower_bound_kg_m2_s"]
    )
    assert pd.isna(oriented.set_index("pl_name").loc["Test b", "angular_momentum_deficit_kg_m2_s"])
    with pytest.raises(ValueError, match="unique"):
        add_measured_host_features(planets, pd.concat([hosts, hosts], ignore_index=True))


@pytest.mark.parametrize('column, value, derived', [
    ('st_rotplim', 1, 'host_rossby_number'), ('st_radlim', -1, 'host_alfven_radius_au'),
    ('st_rotp', np.inf, 'host_rossby_number'), ('st_rad', np.inf, 'host_alfven_radius_au'),
])
def test_host_feature_derivations_exclude_limits_and_nonfinite_stellar_values(column, value, derived):
    planets = enrich_exoplanets(example_planets().iloc[:1].copy())
    planets[column] = value
    before = planets.copy(deep=True)
    measurements = pd.DataFrame({'hostname': ['Test'], 'source_reference': ['synthetic measurements'],
                                 'host_convective_turnover_days': [10.0],
                                 'host_alfven_radius_stellar_radii': [10.0]})
    result = add_measured_host_features(planets, measurements)
    assert np.isnan(result[derived].iloc[0])
    pd.testing.assert_frame_equal(planets, before)


@pytest.mark.parametrize('field', ['pl_name', 'source', 'observation_id', 'source_url'])
@pytest.mark.parametrize('blank', [' ', '\t\n'])
def test_atmosphere_observations_reject_whitespace_provenance(field, blank):
    planets = example_planets()
    observations = pd.DataFrame({
        'pl_name': ['Test b'], 'source': ['synthetic supplied source'],
        'observation_id': ['synthetic-spectrum'],
        'source_url': ['https://example.org/synthetic-spectrum'],
    })
    observations[field] = blank
    original_planets = planets.copy(deep=True)
    original_observations = observations.copy(deep=True)
    with pytest.raises(ValueError, match='cannot be empty'):
        attach_atmosphere_observations(planets, observations)
    pd.testing.assert_frame_equal(planets, original_planets)
    pd.testing.assert_frame_equal(observations, original_observations)


@pytest.mark.parametrize('second_planet', ['Test c', 'Test b'])
def test_atmosphere_observations_keep_multiple_and_unmatched_sources(monkeypatch, second_planet):
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
    second = load_mast_spectra(second_planet)
    assert spectra['observation_id'].tolist() == second['observation_id'].tolist()
    extra = pd.DataFrame({"pl_name": ["Unknown b"], "source": ["EAOT"],
                          "observation_id": ["eaot-1"], "source_url": ["https://example.org/observation"]})
    planets = example_planets()
    supplied = pd.concat([spectra, second, extra], ignore_index=True)
    planet_snapshot, observation_snapshot = planets.copy(deep=True), supplied.copy(deep=True)
    if second_planet == 'Test b':
        with pytest.raises(ValueError, match='unique'):
            attach_atmosphere_observations(planets, supplied)
    else:
        assert not spectra['source_url'].equals(second['source_url'])
        enriched, observations = attach_atmosphere_observations(planets, supplied)
        counts = enriched.set_index('pl_name')['atmosphere_observation_count']
        assert counts['Test b'] == counts['Test c'] == 2
        assert len(observations) == 5 and not observations.loc[4, 'matched_archive_planet']
        assert observations.loc[:3, 'matched_archive_planet'].all()
        assert observations['source_url'].tolist() == supplied['source_url'].tolist()
    pd.testing.assert_frame_equal(planets, planet_snapshot)
    pd.testing.assert_frame_equal(supplied, observation_snapshot)


def test_mast_transmission_reader_preserves_measurements_and_source(monkeypatch):
    from hashlib import sha256
    from io import BytesIO
    import pergamon.exoplanet_catalog as module

    text = ("# Synthetic parser fixture, not observations\n"
            "# Wavelength (microns)  Delta Wavelength (microns)  (Rp/Rs)^2   (Rp/Rs)^2 +/-uncertainty\n"
            "1.5 0.05 -0.001 0.002\n1.2 0.04 0.02 0.003\n")
    payload = text.encode("utf-8")
    calls = []

    def read_source(url, timeout):
        calls.append((url, timeout))
        return BytesIO(payload)

    monkeypatch.setattr(module, "urlopen", read_source)
    result = module.load_mast_transmission_spectrum("Test b", "published spectrum.txt")
    np.testing.assert_array_equal(result.to_numpy(), [[1.5, 0.05, -0.001, 0.002], [1.2, 0.04, 0.02, 0.003]])
    assert list(result) == ["wavelength_um", "delta_wavelength_um", "transit_depth", "transit_depth_error"]
    assert calls == [("https://exo.mast.stsci.edu/api/v0.1/spectra/Test%20b/file/published%20spectrum.txt", 30)]
    assert result.attrs["source_text"] == text
    assert result.attrs["source_sha256"] == sha256(payload).hexdigest()
    assert result.attrs["source_header"].startswith("# Synthetic parser fixture")
    assert result.attrs["pl_name"] == "Test b"
    assert result.attrs["observation_id"] == "published spectrum.txt"
    assert result.attrs["source_url"] == calls[0][0]
    assert result.attrs["retrieved_at_utc"]


@pytest.mark.parametrize("verbosity", [None, "0", "1"])
def test_mast_transmission_reader_narrates_only_when_enabled(monkeypatch, capsys, verbosity):
    from io import BytesIO
    import pergamon.exoplanet_catalog as module

    if verbosity is None:
        monkeypatch.delenv("TDPY_VERBOSITY", raising=False)
    else:
        monkeypatch.setenv("TDPY_VERBOSITY", verbosity)
    text = "# Wavelength (microns) Delta Wavelength (microns) (Rp/Rs)^2 (Rp/Rs)^2 +/-uncertainty\n1 0.1 0.02 0.001"
    monkeypatch.setattr(module, "urlopen", lambda *args, **kwargs: BytesIO(text.encode()))
    result = module.load_mast_transmission_spectrum("Test b", "spectrum.txt")
    assert capsys.readouterr().out == (f"Reading from {result.attrs['source_url']}...\n" if verbosity == "1" else "")


@pytest.mark.parametrize("row", [
    "1 0.1 0.02", "1 0.1 0.02 0.001 5", "nan 0.1 0.02 0.001",
    "-1 0.1 0.02 0.001", "1 -0.1 0.02 0.001", "1 0.1 inf 0.001",
    "1 0.1 0.02 0", "1 0.1 0.02 -0.001", "1 0.1 invalid 0.001", "",
])
def test_mast_transmission_reader_rejects_malformed_values(monkeypatch, row):
    from io import BytesIO
    import pergamon.exoplanet_catalog as module

    header = "# Wavelength (microns) Delta Wavelength (microns) (Rp/Rs)^2 (Rp/Rs)^2 +/-uncertainty\n"
    monkeypatch.setattr(module, "urlopen", lambda *args, **kwargs: BytesIO((header + row).encode()))
    with pytest.raises(ValueError):
        module.load_mast_transmission_spectrum("Test b", "spectrum.txt")


@pytest.mark.parametrize("header", ["", "# Wavelength (nm) Flux Error", "# Wavelength (microns) Flux Error"])
def test_mast_transmission_reader_rejects_unestablished_units_and_columns(monkeypatch, header):
    from io import BytesIO
    import pergamon.exoplanet_catalog as module

    monkeypatch.setattr(module, "urlopen", lambda *args, **kwargs: BytesIO((header + "\n1 0.1 0.02 0.001").encode()))
    with pytest.raises(ValueError, match="header or units"):
        module.load_mast_transmission_spectrum("Test b", "spectrum.txt")


@pytest.mark.parametrize("field", ["planet_name", "observation_id"])
@pytest.mark.parametrize("value", [None, "", " \t", 12])
def test_mast_transmission_reader_validates_identifiers_before_download(monkeypatch, field, value):
    import pergamon.exoplanet_catalog as module

    def unexpected_request(*args, **kwargs):
        pytest.fail("Invalid identifiers reached the network")

    monkeypatch.setattr(module, "urlopen", unexpected_request)
    arguments = {"planet_name": "Test b", "observation_id": "spectrum.txt"}
    arguments[field] = value
    with pytest.raises(ValueError, match="nonempty strings"):
        module.load_mast_transmission_spectrum(**arguments)


@pytest.mark.parametrize("typefileplot", ["png", "pdf"])
@pytest.mark.parametrize("typeplotback", ["white", "dark"])
def test_transmission_plot_preserves_measurements_and_separates_source_legend(
    monkeypatch, tmp_path, typefileplot, typeplotback,
):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.colors import to_rgba
    from PIL import Image
    import pergamon.exoplanet_catalog as module

    spectrum = pd.DataFrame({
        "wavelength_um": [1.5, 0.5, 2.2], "transit_depth": [0.02, -0.001, 0.021],
        "transit_depth_error": [0.002, 0.001, 0.003],
    })
    spectrum.attrs["pl_name"] = "Synthetic layout fixture"
    original = spectrum.copy(deep=True)
    figures = []
    save = Figure.savefig

    def capture(figure, output, **kwargs):
        figures.append(figure)
        save(figure, output, **kwargs)

    monkeypatch.setattr(Figure, "savefig", capture)
    output = module.plot_mast_transmission_spectrum(
        spectrum, tmp_path / "spectrum", source_label="Explicit synthetic parser/layout fixture",
        typefileplot=typefileplot, typeplotback=typeplotback,
    )
    pd.testing.assert_frame_equal(spectrum, original)
    figure = figures[0]
    canvas = FigureCanvasAgg(figure)
    canvas.draw()
    upper, lower = figure.axes
    assert upper.get_shared_x_axes().joined(upper, lower)
    assert upper.get_xscale() == lower.get_xscale() == "log"
    np.testing.assert_array_equal(upper.lines[0].get_xdata(), spectrum["wavelength_um"])
    np.testing.assert_array_equal(upper.lines[0].get_ydata(), spectrum["transit_depth"] * 1e6)
    np.testing.assert_array_equal(lower.lines[0].get_ydata(), spectrum["transit_depth_error"] * 1e6)
    assert upper.lines[0].get_linestyle() == "None"
    foreground = "white" if typeplotback == "dark" else "black"
    np.testing.assert_allclose(to_rgba(upper.lines[0].get_color()), to_rgba(foreground))
    for axis in figure.axes:
        assert not any(line.get_visible() for line in axis.get_xgridlines() + axis.get_ygridlines())
    segments = upper.collections[0].get_segments()
    for segment, wavelength, depth, error in zip(segments, spectrum["wavelength_um"],
                                               spectrum["transit_depth"], spectrum["transit_depth_error"]):
        np.testing.assert_allclose(segment, [[wavelength, (depth - error) * 1e6],
                                            [wavelength, (depth + error) * 1e6]])
    legend = upper.get_legend()
    assert legend.get_frame().get_alpha() == 1
    assert "synthetic" in legend.get_texts()[0].get_text()
    legend_box = legend.get_window_extent(canvas.get_renderer())
    assert not legend_box.overlaps(upper.get_window_extent())
    assert not legend_box.overlaps(upper._left_title.get_window_extent(canvas.get_renderer()))
    assert upper.get_ylabel().endswith("[ppm]")
    assert lower.get_xlabel().endswith("[micron]")
    if typefileplot == "png":
        with Image.open(output) as image:
            assert image.info["dpi"] == pytest.approx((300, 300), abs=0.01)
    else:
        assert output.read_bytes().startswith(b"%PDF-")


@pytest.mark.parametrize("column,value", [
    ("wavelength_um", 0), ("wavelength_um", np.nan),
    ("transit_depth", np.inf), ("transit_depth", 1 + 2j),
    ("transit_depth_error", 0), ("transit_depth_error", -0.01),
])
def test_transmission_plot_rejects_invalid_data_before_writing(tmp_path, column, value):
    import pergamon.exoplanet_catalog as module

    spectrum = pd.DataFrame({"wavelength_um": [1.0], "transit_depth": [0.02], "transit_depth_error": [0.001]})
    spectrum[column] = [value]
    with pytest.raises(ValueError, match="finite real"):
        module.plot_mast_transmission_spectrum(spectrum, tmp_path / "unused/figure", source_label="Synthetic fixture")
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("source_label", [None, "", " \t"])
def test_transmission_plot_requires_source_label(tmp_path, source_label):
    import pergamon.exoplanet_catalog as module

    with pytest.raises(ValueError, match="source_label"):
        module.plot_mast_transmission_spectrum(pd.DataFrame(), tmp_path / "unused/figure", source_label=source_label)
    assert list(tmp_path.iterdir()) == []


def test_transmission_plot_preserves_existing_output_and_rejects_unsupported_format(tmp_path):
    import pergamon.exoplanet_catalog as module

    spectrum = pd.DataFrame({"wavelength_um": [1.0], "transit_depth": [0.02], "transit_depth_error": [0.001]})
    output = tmp_path / "previous.png"
    output.write_bytes(b"preserved prior artifact")
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        module.plot_mast_transmission_spectrum(spectrum, output, source_label="Synthetic fixture")
    with pytest.raises(ValueError, match="typefileplot"):
        module.plot_mast_transmission_spectrum(spectrum, tmp_path / "unused/figure", source_label="Synthetic fixture", typefileplot="svg")
    assert output.read_bytes() == b"preserved prior artifact"
    assert sorted(tmp_path.iterdir()) == [output]


def test_transmission_readme_links_source_cited_products():
    import re
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    readme = (root / "README.md").read_text()
    links = re.findall(r"\]\((examples/visuals/wasp39b_wakeford2018_transmission_uncertainty_[^)]+)\)", readme)
    assert len(links) == 2
    assert {Path(link).suffix for link in links} == {".png", ".pdf"}
    assert all((root / link).is_file() for link in links)
    assert "Wakeford et al. (2018)" in readme
    assert "load_mast_transmission_spectrum" in readme
    assert "plot_mast_transmission_spectrum" in readme


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