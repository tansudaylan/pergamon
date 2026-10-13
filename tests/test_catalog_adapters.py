import json

import numpy as np
import pytest

from pergamon import catalogs


def test_tic_catalog_first_read_and_cached_read_have_identical_value_label_contract(monkeypatch, tmp_path):
    monkeypatch.setenv('PERGAMON_PATH', str(tmp_path))
    monkeypatch.delenv('EPHESOS_PATH', raising=False)
    rows = [{'ID': 12, 'ra': 10.0, 'dec': -20.0, 'Tmag': 5.0, 'rad': 0.8, 'mass': 0.9}]
    calls = []

    def query(request):
        calls.append(request)
        return [], json.dumps({'data': rows})

    monkeypatch.setattr(catalogs, 'quer_mast', query)
    first = catalogs.retr_dictpopltic8('TIC_m060', typeverb=0)
    cached = catalogs.retr_dictpopltic8('TIC_m060', typeverb=0)

    assert len(calls) == 1
    assert first.keys() == cached.keys()
    for name in first:
        np.testing.assert_array_equal(first[name][0], cached[name][0])
        assert first[name][1] == cached[name][1]
    assert first['radistar'][1] == r'$R_{\odot}$'


def test_host_catalog_groups_planets_and_roundtrips_companion_arrays(monkeypatch, tmp_path):
    monkeypatch.setenv('PERGAMON_PATH', str(tmp_path))
    monkeypatch.delenv('EPHESOS_PATH', raising=False)
    values = {
        'namestar': ['host', 'host', 'other'], 'radistar': [1.0, 1.0, 0.5],
        'massstar': [1.0, 1.0, 0.5], 'booltran': [True, False, True],
        'epocmtracomp': [1.0, 2.0, 3.0], 'pericomp': [3.0, 10.0, 4.0],
        'duratrantotl': [1.0, np.nan, 2.0], 'radicomp': [1.0, 2.0, 1.5],
        'masscomp': [1.0, 8.0, 3.0],
    }
    monkeypatch.setattr(catalogs, 'retr_dictexar', lambda: {name: [np.asarray(data), ''] for name, data in values.items()})

    first = catalogs.retr_dicthostplan('exar', typeverb=0)
    cached = catalogs.retr_dicthostplan('exar', typeverb=0)

    np.testing.assert_array_equal(first['numbplanstar'], [2, 1])
    np.testing.assert_array_equal(first['numbplantranstar'], [1, 1])
    for actual, expected in zip(cached['pericomp'], first['pericomp']):
        np.testing.assert_array_equal(actual, expected)


def test_toi_matching_and_faint_star_selection_use_value_arrays(monkeypatch):
    catalog = {'TICID': [np.array([12, 34]), ''], 'TOIID': [np.array([8.01, 9.01]), ''],
               'strgcomm': [np.array(['found in faint-star QLP search', 'other']), ''],
               'nametoii': [np.array(['TOI-8.01', 'TOI-9.01']), '']}
    monkeypatch.setattr(catalogs, 'retr_dicttoii', lambda **kwargs: catalog)

    assert catalogs.retr_toiitici(12, typeverb=0) == 8
    np.testing.assert_array_equal(catalogs.retr_toiifstr(), ['8.01'])


def test_missing_local_target_returns_none_before_derived_calculations(monkeypatch):
    import pandas as pd

    monkeypatch.setattr(catalogs.pd, 'read_csv', lambda *args, **kwargs: pd.DataFrame({'hostname': ['known']}))
    assert catalogs.retr_dictexar('unknown', typeverb=0) is None