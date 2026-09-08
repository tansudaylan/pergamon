import numpy as np

import pergamon


def test_init_accepts_flattened_array_population():
    dictpopl = {
        'pop': {
            'radistar': np.array([1.0, 1.2]),
            'massstar': np.array([1.0, 1.1]),
        }
    }

    result = pergamon.init(
        typeanls='defa',
        dictpopl=dictpopl,
        boolmakeplot=False,
        booldiag=False,
    )

    assert 'pop' in result
    np.testing.assert_allclose(result['pop']['radistar'][0], np.array([1.0, 1.2]))
    np.testing.assert_allclose(result['pop']['massstar'][0], np.array([1.0, 1.1]))
