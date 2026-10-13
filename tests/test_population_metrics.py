import numpy as np
from pergamon import emission_spectroscopy as metrics


def test_missing_planet_metrics_remain_nan_without_runtime_warning():
    values = {
        'masscomp': [1.0, 1.0],
        'radicomp': [2.0, 2.0],
        'tmptcomp': [800.0, np.nan],
        'radistar': [1.0, 1.0],
        'magtsystJbnd': [10.0, 10.0],
        'magtsystKbnd': [10.0, 10.0],
        'tmptstar': [5700.0, 5700.0],
    }
    population = {name: [np.asarray(data)] for name, data in values.items()}

    with np.errstate(all='raise'):
        metrics.calc_tsmmesmm(population)

    for name in ('tsmm', 'stdvtsmm', 'esmm', 'stdvesmm'):
        assert np.isfinite(population[name][0][0])
        assert np.isnan(population[name][0][1])
