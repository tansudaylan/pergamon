from urllib.parse import parse_qs, urlparse

import numpy as np
import pandas as pd
import pytest

from pergamon import gaia_kinematics as gaia


def _astrometry(identifiers):
    frame = pd.DataFrame({name: np.ones(len(identifiers)) for name in gaia.GAIA_COLUMNS})
    frame['source_id'] = pd.Series(identifiers, dtype='string')
    frame['ref_epoch'] = 2016.0
    frame['gaia_source_url'] = gaia.GAIA_TAP_URL
    frame['gaia_retrieved_at_utc'] = '2026-10-23T00:00:00+00:00'
    return frame


def test_gaia_fetch_preserves_large_ids_deduplicates_and_records_provenance(monkeypatch):
    identifiers = ['5853498713190525696', '2260370629320898304']
    queries = []

    def read_csv(url, dtype):
        assert dtype == {'source_id': 'string'}
        query = parse_qs(urlparse(url).query)['QUERY'][0]
        queries.append(query)
        selected = [identifier for identifier in identifiers if identifier in query]
        return _astrometry(selected)

    monkeypatch.setattr(gaia.pd, 'read_csv', read_csv)
    result = gaia.load_gaia_dr3_astrometry([
        f'Gaia DR3 {identifiers[0]}', identifiers[0], identifiers[1], None,
    ], batch_size=1)
    assert len(queries) == 2
    assert set(result['source_id']) == set(identifiers)
    assert result['gaia_source_url'].str.startswith(gaia.GAIA_TAP_URL).all()
    assert result['gaia_retrieved_at_utc'].str.endswith('+00:00').all()


@pytest.mark.parametrize('identifier', [1.0, 5853498713190525696.0, True, 'Gaia DR2 123', '1);select x', '-2', '0'])
def test_gaia_fetch_rejects_inexact_or_invalid_ids_before_network(monkeypatch, identifier):
    monkeypatch.setattr(gaia.pd, 'read_csv', lambda *args, **kwargs: pytest.fail('Network request for invalid ID'))
    with pytest.raises(ValueError):
        gaia.load_gaia_dr3_astrometry([identifier])


def test_gaia_fetch_rejects_unrequested_response_ids(monkeypatch):
    monkeypatch.setattr(gaia.pd, 'read_csv', lambda *args, **kwargs: _astrometry(['999']))
    with pytest.raises(ValueError, match='unrequested'):
        gaia.load_gaia_dr3_astrometry(['123'])


def test_gaia_empty_request_keeps_schema_without_network(monkeypatch):
    monkeypatch.setattr(gaia.pd, 'read_csv', lambda *args, **kwargs: pytest.fail('Unexpected network request'))
    result = gaia.load_gaia_dr3_astrometry([None, pd.NA, np.nan, ''])
    assert result.empty and set(gaia.GAIA_COLUMNS).issubset(result)


def test_phase_space_preserves_planet_rows_indices_and_velocity_units():
    identifier = '5853498713190525696'
    planets = pd.DataFrame({'pl_name': ['a b', 'a c', 'unknown b'],
                            'gaia_dr3_id': [identifier, 'Gaia DR3 ' + identifier, None]}, index=[9, 9, 3])
    measured = _astrometry([identifier])
    measured.loc[0, ['ra', 'dec', 'parallax', 'parallax_error', 'pmra', 'pmdec', 'radial_velocity']] = [0, 0, 100, 1, 100, 0, 20]
    original = planets.copy(deep=True)
    result = gaia.attach_gaia_phase_space(planets, measured)
    pd.testing.assert_frame_equal(planets, original)
    assert result.index.tolist() == [9, 9, 3]
    assert result['gaia_6d_status'].tolist() == ['complete', 'complete', 'unmatched']
    assert result['pl_name'].tolist() == planets['pl_name'].tolist()
    position = result.iloc[0][['gaia_barycentric_galactic_' + axis + '_pc' for axis in ('x', 'y', 'z')]].to_numpy(dtype=float)
    velocity = result.iloc[0][['gaia_barycentric_galactic_v' + axis + '_km_s' for axis in ('x', 'y', 'z')]].to_numpy(dtype=float)
    assert np.linalg.norm(position) == pytest.approx(10)
    assert np.linalg.norm(velocity) == pytest.approx(np.hypot(20, 4.74047046), rel=1e-7)
    rotation = np.array([
        [-0.0548755604, -0.8734370902, -0.4838350155],
        [0.4941094279, -0.4448296300, 0.7469822445],
        [-0.8676661490, -0.1980763734, 0.4559837762],
    ])
    np.testing.assert_allclose(position, rotation @ [10, 0, 0], atol=1e-6)
    np.testing.assert_allclose(velocity, rotation @ [20, 4.74047046, 0], atol=2e-6)
    assert np.isnan(result['gaia_barycentric_galactic_vx_km_s'].iloc[-1])


@pytest.mark.parametrize('field,value,status', [
    ('radial_velocity', np.nan, 'incomplete_astrometry'),
    ('parallax', -1, 'invalid_parallax'), ('parallax_error', 0.3, 'parallax_precision_cut'),
    ('dec', 100, 'invalid_position'), ('pmra', np.inf, 'incomplete_astrometry'),
])
def test_incomplete_or_unsuitable_astrometry_never_fabricates_six_dimensions(field, value, status):
    measured = _astrometry(['123'])
    measured.loc[0, field] = value
    result = gaia.attach_gaia_phase_space(pd.DataFrame({'gaia_dr3_id': ['123']}), measured)
    assert result['gaia_6d_status'].iloc[0] == status
    assert not result['gaia_6d_available'].iloc[0]
    assert np.isnan(result['gaia_barycentric_galactic_vx_km_s'].iloc[0])


def test_phase_space_rejects_duplicate_measurements_and_missing_provenance():
    planets = pd.DataFrame({'gaia_dr3_id': ['123']})
    with pytest.raises(ValueError, match='unique'):
        gaia.attach_gaia_phase_space(planets, _astrometry(['123', '123']))
    measured = _astrometry(['123'])
    measured.loc[0, 'gaia_source_url'] = ''
    with pytest.raises(ValueError, match='provenance'):
        gaia.attach_gaia_phase_space(planets, measured)


@pytest.mark.parametrize('name', ['gaia_ra', 'gaia_6d_status', 'gaia_barycentric_galactic_vx_km_s'])
def test_phase_space_never_overwrites_existing_measurements_or_results(name):
    planets = pd.DataFrame({'gaia_dr3_id': ['123'], name: ['existing']})
    with pytest.raises(ValueError, match='already contains'):
        gaia.attach_gaia_phase_space(planets, _astrometry(['123']))
    assert planets[name].iloc[0] == 'existing'


@pytest.mark.parametrize('batch_size', [0, 101, 1.0, True])
def test_gaia_batch_size_is_bounded_before_network(monkeypatch, batch_size):
    monkeypatch.setattr(gaia.pd, 'read_csv', lambda *args, **kwargs: pytest.fail('Unexpected network request'))
    with pytest.raises(ValueError, match='batch_size'):
        gaia.load_gaia_dr3_astrometry(['123'], batch_size=batch_size)


def test_parallax_precision_boundary_and_stricter_user_threshold():
    planets = pd.DataFrame({'gaia_dr3_id': ['123', '456']})
    measured = _astrometry(['123', '456'])
    measured['parallax_error'] = [0.2, np.nextafter(0.2, 1)]
    result = gaia.attach_gaia_phase_space(planets, measured)
    assert result['gaia_6d_status'].tolist() == ['complete', 'parallax_precision_cut']
    stricter = gaia.attach_gaia_phase_space(planets, measured, max_fractional_parallax_error=0.1)
    assert not stricter['gaia_6d_available'].any()


@pytest.mark.parametrize('threshold', [0, -0.1, 0.3, True, np.nan, np.inf])
def test_rejects_invalid_parallax_threshold(threshold):
    with pytest.raises(ValueError, match='threshold'):
        gaia.attach_gaia_phase_space(pd.DataFrame({'gaia_dr3_id': ['123']}),
                                    _astrometry(['123']), max_fractional_parallax_error=threshold)


def _rank_planets():
    planets = pd.DataFrame({
        'pl_name': ['First b', 'First c', 'Second b', 'Third b'],
        'gaia_dr3_id': ['123', '123', '456', '789'],
        'pl_orbper': [1.0, 3.0, 4.0, 6.0],  # [day]
        'pl_orbperlim': [0, 0, 0, 0],
    })
    measured = _astrometry(['123', '456', '789'])
    measured['parallax_error'] = 0.01  # [mas]
    measured['parallax'] = [100.0, 50.0, 25.0]  # [mas]
    return gaia.attach_gaia_phase_space(planets, measured)


def test_rank_correlation_counts_hosts_not_planets_and_retains_provenance():
    planets = _rank_planets()
    original = planets.copy(deep=True)
    result = gaia.host_phase_space_rank_correlation(planets, 'pl_orbper', 'gaia_inverse_parallax_distance_pc')
    pd.testing.assert_frame_equal(planets, original)
    assert result['n_hosts'] == 3 and result['n_used_planets'] == 4
    assert result['spearman_rho'] == pytest.approx(1)
    assert result['status'] == 'descriptive_only'
    assert result['hosts']['planet_value'].tolist() == [2, 4, 6]
    assert result['hosts']['measured_planet_count'].tolist() == [2, 1, 1]
    assert result['hosts']['gaia_source_url'].notna().all()
    single_planets = planets.drop(index=0)
    single = gaia.host_phase_space_rank_correlation(single_planets, 'pl_orbper', 'gaia_inverse_parallax_distance_pc')
    assert single['n_hosts'] == result['n_hosts']


def test_rank_correlation_excludes_limits_and_keeps_small_sample_result_undefined():
    planets = _rank_planets()
    planets.loc[3, 'pl_orbperlim'] = 1
    result = gaia.host_phase_space_rank_correlation(planets, 'pl_orbper', 'gaia_inverse_parallax_distance_pc')
    assert result['n_hosts'] == 2 and result['n_used_planets'] == 3
    assert np.isnan(result['spearman_rho'])
    assert result['status'] == 'insufficient_hosts'


def test_rank_correlation_rejects_conflicting_host_measurements_and_duplicate_planets():
    planets = _rank_planets()
    planets.loc[1, 'gaia_inverse_parallax_distance_pc'] = 20
    with pytest.raises(ValueError, match='conflicting'):
        gaia.host_phase_space_rank_correlation(planets, 'pl_orbper', 'gaia_inverse_parallax_distance_pc')
    planets = _rank_planets()
    planets.loc[1, 'pl_name'] = 'First b'
    with pytest.raises(ValueError, match='unique'):
        gaia.host_phase_space_rank_correlation(planets, 'pl_orbper', 'gaia_inverse_parallax_distance_pc')


@pytest.mark.parametrize('value', [np.nan, np.inf, 0, -1])
def test_rank_comparison_excludes_invalid_planet_values(value):
    planets = _rank_planets()
    planets.loc[3, 'pl_orbper'] = value
    result = gaia.host_phase_space_rank_correlation(planets, 'pl_orbper', 'gaia_inverse_parallax_distance_pc')
    assert result['n_hosts'] == 2 and result['n_used_planets'] == 3


def test_rank_comparison_keeps_constant_empty_and_unavailable_cases_explicit():
    planets = _rank_planets()
    planets['pl_orbper'] = 1.0  # [day]
    result = gaia.host_phase_space_rank_correlation(planets, 'pl_orbper', 'gaia_inverse_parallax_distance_pc')
    assert result['status'] == 'constant_values' and np.isnan(result['spearman_rho'])
    planets.loc[3, 'gaia_6d_available'] = False
    result = gaia.host_phase_space_rank_correlation(planets, 'pl_orbper', 'gaia_inverse_parallax_distance_pc')
    assert result['n_hosts'] == 2
    result = gaia.host_phase_space_rank_correlation(planets.iloc[:0], 'pl_orbper', 'gaia_inverse_parallax_distance_pc')
    assert result['n_hosts'] == 0 and result['status'] == 'insufficient_hosts'


@pytest.mark.parametrize('planet_column,host_column', [
    ('pl_bmasse', 'gaia_inverse_parallax_distance_pc'), ('pl_orbper', 'pl_orbper'),
])
def test_rank_comparison_rejects_unsupported_or_circular_comparisons(planet_column, host_column):
    with pytest.raises(ValueError, match='Choose'):
        gaia.host_phase_space_rank_correlation(_rank_planets(), planet_column, host_column)


@pytest.mark.parametrize('theme,extension', [('white', 'png'), ('dark', 'png'), ('white', 'pdf')])
def test_rank_plot_formats_themes_bounds_and_timestamp_protection(tmp_path, monkeypatch, theme, extension):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.colors import to_rgb
    from PIL import Image
    from tdpy import plotting

    comparison = gaia.host_phase_space_rank_correlation(_rank_planets(), 'pl_orbper', 'gaia_inverse_parallax_distance_pc')
    original_save = plotting.save_plot

    def save(figure, path, **kwargs):
        FigureCanvasAgg(figure)
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        axis = figure.axes[0]
        for text in [axis.title, axis.xaxis.label, axis.yaxis.label, *figure.texts]:
            bounds = text.get_window_extent(renderer)
            assert bounds.x0 >= 0 and bounds.y0 >= 0
            assert bounds.x1 <= figure.bbox.width and bounds.y1 <= figure.bbox.height
        expected = 'black' if theme == 'dark' else 'white'
        np.testing.assert_allclose(figure.get_facecolor()[:3], to_rgb(expected))
        original_save(figure, path, **kwargs)

    monkeypatch.setattr(plotting, 'save_plot', save)
    paths = gaia.plot_host_phase_space_rank_correlation(
        comparison, tmp_path, typefileplot=extension, typeplotback=theme, time_tag='20261024_010000',
    )
    if extension == 'png':
        image = np.asarray(Image.open(paths['static']).convert('RGB'))
        assert image.shape[0] > 1000 and image.shape[1] > 1000
        assert np.any(image != image[0, 0])
    else:
        assert paths['static'].read_bytes().startswith(b'%PDF')
    before = paths['static'].read_bytes()
    with pytest.raises(FileExistsError):
        gaia.plot_host_phase_space_rank_correlation(
            comparison, tmp_path, typefileplot=extension, typeplotback=theme, time_tag='20261024_010000',
        )
    assert paths['static'].read_bytes() == before


def test_bokeh_plot_uses_the_exact_same_host_values(tmp_path, monkeypatch):
    embed = pytest.importorskip('bokeh.embed')
    original_html = embed.file_html
    comparison = gaia.host_phase_space_rank_correlation(_rank_planets(), 'pl_orbper', 'gaia_inverse_parallax_distance_pc')

    def html(plot, *args, **kwargs):
        source = plot.renderers[0].data_source.data
        np.testing.assert_array_equal(source['host_value'], comparison['hosts']['host_value'])
        np.testing.assert_array_equal(source['planet_value'], comparison['hosts']['planet_value'])
        assert list(source['gaia_source_id']) == ['123', '456', '789']
        assert not any(grid.visible for grid in plot.grid)
        return original_html(plot, *args, **kwargs)

    monkeypatch.setattr(embed, 'file_html', html)
    paths = gaia.plot_host_phase_space_rank_correlation(comparison, tmp_path, write_bokeh=True)
    text = paths['bokeh'].read_text()
    assert 'Barycentric distance [pc]' in text and 'Host median orbital period [day]' in text


def test_empty_comparison_does_not_produce_a_blank_figure(tmp_path):
    comparison = gaia.host_phase_space_rank_correlation(_rank_planets().iloc[:0], 'pl_orbper', 'gaia_inverse_parallax_distance_pc')
    with pytest.raises(ValueError, match='No eligible hosts'):
        gaia.plot_host_phase_space_rank_correlation(comparison, tmp_path)
    assert list(tmp_path.iterdir()) == []


def test_signed_velocity_ranks_are_not_filtered_like_positive_planet_properties():
    planets = _rank_planets()
    planets['gaia_barycentric_galactic_vz_km_s'] = [-30, -30, -10, 20]  # [km/s]
    result = gaia.host_phase_space_rank_correlation(planets, 'pl_orbper', 'gaia_barycentric_galactic_vz_km_s')
    assert result['n_hosts'] == 3 and result['spearman_rho'] == pytest.approx(1)


def test_rank_comparison_rejects_string_availability_flags():
    planets = _rank_planets()
    planets['gaia_6d_available'] = 'False'
    with pytest.raises(ValueError, match='Boolean'):
        gaia.host_phase_space_rank_correlation(planets, 'pl_orbper', 'gaia_inverse_parallax_distance_pc')


def test_notebook_snapshot_replay_never_uses_network_and_preserves_ids(tmp_path, monkeypatch):
    import json
    from pathlib import Path
    import IPython.display
    from pergamon import exoplanet_catalog

    path = Path(__file__).parents[1] / 'examples' / 'gaia_host_correlations.ipynb'
    notebook = json.loads(path.read_text())
    source = next(''.join(cell['source']) for cell in notebook['cells'] if cell['cell_type'] == 'code')
    examples = tmp_path / 'examples'
    examples.mkdir()
    (examples / path.name).write_text('{}')
    data = tmp_path / 'data'
    data.mkdir()
    tag = '20261024_010000'
    identifier = '5853498713190525696'
    pd.DataFrame({'pl_name': ['Test b', 'Unknown b'], 'gaia_dr3_id': [identifier, None]}).to_csv(
        data / f'gaia_archive_planets_{tag}.csv', index=False,
    )
    measured = _astrometry([identifier])
    measured['parallax_error'] = 0.01  # [mas]
    measured.to_csv(data / f'gaia_dr3_astrometry_{tag}.csv', index=False)
    monkeypatch.chdir(examples)
    monkeypatch.setenv('PERGAMON_GAIA_SNAPSHOT_TAG', tag)
    monkeypatch.setattr(IPython.display, 'display', lambda *args: None)
    def no_network(*args, **kwargs):
        pytest.fail('Snapshot replay attempted a network request')
    monkeypatch.setattr(exoplanet_catalog, 'load_confirmed_exoplanets', no_network)
    monkeypatch.setattr(gaia, 'load_gaia_dr3_astrometry', no_network)
    namespace = {}
    exec(compile(source, str(path), 'exec'), namespace)
    assert namespace['planets']['gaia_dr3_id'].iloc[0] == identifier
    assert pd.isna(namespace['planets']['gaia_dr3_id'].iloc[1])
    assert namespace['joined']['gaia_6d_status'].tolist() == ['complete', 'unmatched']
    assert namespace['astrometry']['gaia_retrieved_at_utc'].iloc[0] == measured['gaia_retrieved_at_utc'].iloc[0]