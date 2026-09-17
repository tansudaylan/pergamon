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


def test_init_falls_back_when_env_data_path_is_missing(monkeypatch):
    monkeypatch.delenv('PERGAMON_DATA_PATH', raising=False)

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


def test_retr_subp_accepts_python_list_indices():
    dictpopl = {'pop': {'radistar': [np.array([1.0, 2.0, 3.0]), '']}}
    dictnumbsamp = {}
    dictindxsamp = {}

    pergamon.retr_subp(dictpopl, 'pop', 'small', [0, 2], dictnumbsamp, dictindxsamp)

    np.testing.assert_array_equal(dictpopl['small']['radistar'][0], np.array([1.0, 3.0]))
    np.testing.assert_array_equal(dictindxsamp['pop']['small'], np.array([0, 2]))
    assert dictnumbsamp['small'] == 2
