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