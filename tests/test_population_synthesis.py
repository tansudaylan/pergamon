import numpy as np
import pytest

from pergamon.population_synthesis import retr_dictpoplstarcomp


@pytest.mark.parametrize('radius', [
    0.0, -0.8, np.nan, np.inf, -np.inf, True,
    np.array(True), 0.8 + 0.1j, np.ma.array(0.8, mask=True), np.array([0.8, 0.9]),
])
def test_synthetic_stellar_radius_rejects_invalid_values_before_sampling(monkeypatch, radius):
    from importlib import import_module

    model = import_module('pergamon.population_synthesis')
    original_data = np.ma.getdata(radius).copy()
    original_mask = np.ma.getmaskarray(radius).copy()
    monkeypatch.setattr(model.np.random, 'rand',
                        lambda *args: pytest.fail('Invalid radius must fail before sampling'))

    with pytest.raises(ValueError, match='radistar'):
        retr_dictpoplstarcomp('PlanetarySystem', numbsyst=3, typestar='drawkrou', radistar=radius)

    np.testing.assert_array_equal(np.ma.getdata(radius), original_data)
    np.testing.assert_array_equal(np.ma.getmaskarray(radius), original_mask)


@pytest.mark.parametrize('radius', [0.8, np.float32(0.8), np.array(0.8), np.ma.array(0.8, mask=False)])
def test_synthetic_stellar_radius_preserves_supported_real_scalars(radius):
    random_state = np.random.get_state()
    try:
        np.random.seed(42)
        result = retr_dictpoplstarcomp(
            'PlanetarySystem', numbsyst=3, typestar='drawkrou', radistar=radius,
            typesamporbtcomp='peripowr', minmpericomp=3.0, maxmpericomp=10.0,
            booltrancomp=True,
        )
    finally:
        np.random.set_state(random_state)

    stars = result['dictpopl']['star']['star_SyntheticPopulation_All']
    np.testing.assert_array_equal(stars['radistar'][0], np.full(3, float(radius)))
    np.testing.assert_allclose(stars['densstar'][0],
                               1.4 * stars['massstar'][0] / stars['radistar'][0]**3)


@pytest.mark.parametrize('stellar_radius', [None, 0.8])
def test_mass_drawn_population_preserves_mass_radius_density_relation(stellar_radius):
    options = {} if stellar_radius is None else {'radistar': stellar_radius}
    state = np.random.get_state()
    try:
        np.random.seed(42)
        result = retr_dictpoplstarcomp(
            'PlanetarySystem', numbsyst=3, typestar='drawkrou',
            typesamporbtcomp='peripowr', minmpericomp=3.0,
            maxmpericomp=10.0, booltrancomp=True, **options,
        )
    finally:
        np.random.set_state(state)

    stars = result['dictpopl']['star']['star_SyntheticPopulation_All']
    masses = stars['massstar'][0]
    radii = stars['radistar'][0]
    densities = stars['densstar'][0]
    np.testing.assert_allclose(densities, 1.4 * masses / radii**3)
    if stellar_radius is None:
        np.testing.assert_allclose(densities, 1.4 * masses**(-0.7))
    else:
        np.testing.assert_array_equal(radii, np.full(3, stellar_radius))
    assert np.isfinite(radii).all() and (radii > 0).all()
    assert stars['radistar'][1] == r'$R_{\odot}$'
    assert stars['densstar'][1] == 'g cm$^{-3}$'


def test_seeded_planet_population_preserves_verified_periods_and_units(monkeypatch):
    from unittest.mock import Mock
    import tdpy

    shared_values = Mock(wraps=tdpy.population_feature_values)
    monkeypatch.setattr(tdpy, 'population_feature_values', shared_values)
    state = np.random.get_state()
    try:
        np.random.seed(42)
        result = retr_dictpoplstarcomp(
            'PlanetarySystem', numbsyst=3, typesamporbtcomp='peripowr',
            minmpericomp=3.0, maxmpericomp=10.0, booltrancomp=True,
        )
    finally:
        np.random.set_state(state)

    stars = result['dictpopl']['star']['star_SyntheticPopulation_All']
    planets = result['dictpopl']['comp']['compstar_SyntheticPopulation_All']
    np.testing.assert_allclose(planets['pericomp'][0], [3.32463646, 4.68723713, 8.93441805])
    np.testing.assert_array_equal(stars['massstar'][0], np.ones(3))
    assert stars['radistar'][1] == r'$R_{\odot}$'
    assert planets['booltran'][0].all()
    assert result['dictnumbsamp']['comp']['compstar_SyntheticPopulation_All'] == 3
    assert shared_values.call_count > 0
    star_alias = result['dictpopl']['star']['starSynthetictotl']
    planet_alias = result['dictpopl']['comp']['compstarSynthetictotl']
    assert star_alias['massstar'] is stars['massstar'][0]
    assert planet_alias['pericomp'] is planets['pericomp'][0]
    assert len(stars['massstar']) == len(planets['pericomp']) == 2