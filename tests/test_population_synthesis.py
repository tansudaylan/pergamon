import numpy as np

from pergamon.population_synthesis import retr_dictpoplstarcomp


def test_seeded_planet_population_preserves_verified_periods_and_units():
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