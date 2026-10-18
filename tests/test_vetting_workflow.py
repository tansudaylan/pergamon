from types import SimpleNamespace

import numpy as np
import pytest

from pergamon.vetting_workflow import synthetic_target_truth, validated_classification


def test_synthetic_truth_uses_each_class_local_target_including_negatives():
    def population(first, second):
        return {
            'listnamefeatbody': ['massstar'], 'listnamefeatlimbonly': ['masscomp'],
            'dictpopl': {
                'star': {'star_SyntheticPopulation_All': {'massstar': [np.array(first), '']}},
                'comp': {'compstar_SyntheticPopulation_All': {'masscomp': [np.array(second), '']}},
            },
        }

    state = SimpleNamespace(
        listnameclastruetype=['PlanetarySystem', 'StellarBinary'],
        dictindxtarg={'PlanetarySystem': np.array([0, 1]), 'StellarBinary': np.array([2, 3])},
        dicttroy={'true': {'PlanetarySystem': population([1, 2], [10, 20]),
                           'StellarBinary': population([3, 4], [30, 40])}},
        namepoplstartotl='star_SyntheticPopulation_All',
        namepoplcomptotl='compstar_SyntheticPopulation_All',
    )
    assert [synthetic_target_truth(state, index)['masscomp'].tolist() for index in range(4)] == [[10], [20], [30], [40]]
    assert synthetic_target_truth(state, 3)['typemodl'] == 'StellarBinary'


def test_miletos_classification_requires_real_complete_results():
    result = {'boolcalclspe': False, 'boolsrchboxsperi': False,
              'boolsrchoutlperi': True, 'boolposianls': [True, False],
              'dictoutlperi': {'minmfrddtimeoutlsort': [0.07], 'boolposi': True}}
    assert validated_classification(result, 2).tolist() == [True, False]
    for invalid in ({}, {**result, 'dictoutlperi': {}},
                    {**result, 'dictoutlperi': {'boolposi': np.nan}},
                    {**result, 'boolsrchoutlperi': False, 'boolposianls': [True]}):
        with pytest.raises(ValueError, match='classification'):
            validated_classification(invalid, 2)


@pytest.mark.parametrize('lspe, boxs, expected', [
    (False, True, [False, True]),
    (True, False, [True, False]),
    (True, True, [True, False, False, True]),
])
def test_detector_classification_uses_only_enabled_miletos_slots(lspe, boxs, expected):
    result = {'boolcalclspe': lspe, 'boolsrchboxsperi': boxs,
              'boolsrchoutlperi': False,
              'boolposianls': np.array([False, True, False, True])}
    assert validated_classification(result, len(expected)).tolist() == expected
    with pytest.raises(ValueError, match='classification'):
        validated_classification({**result, 'boolposianls': [False]}, len(expected))


def test_balanced_worker_analyzes_each_synthetic_class(monkeypatch):
    import pergamon.vetting_workflow as main

    def population(stellar_mass):
        return {'listnamefeatbody': ['massstar'], 'listnamefeatlimbonly': ['masscomp'],
                'dictpopl': {'star': {'star_SyntheticPopulation_All': {'massstar': [np.array([stellar_mass]), '']}},
                             'comp': {'compstar_SyntheticPopulation_All': {'masscomp': [np.array([5.0]), '']}}}}

    state = SimpleNamespace(
        listindxtarg=[np.array([0, 1])], boolsimusome=True, indxtypeclastrue=range(1),
        boolreletarg=[[False, False]], dictindxtarg={'rele': [np.array([0])],
                                                       'PlanetarySystem': np.array([0]),
                                                       'StellarBinary': np.array([1])},
        listnameclastruetype=['PlanetarySystem', 'StellarBinary'],
        dicttroy={'true': {'PlanetarySystem': population(1.0),
                           'StellarBinary': population(2.0)}},
        namepoplstartotl='star_SyntheticPopulation_All',
        namepoplcomptotl='compstar_SyntheticPopulation_All',
        typepopl='SyntheticPopulation', labltarg=['planet', 'binary'],
        strgtarg=['planet', 'binary'], numbtarget=2, numbtarg=2,
        numbtargplot=0, boolplotmile=False, numbtarganimtrue=0,
        dictmileinptglob={}, indxdatatser=[], listlablinst=[],
    )
    calls = []

    def analyze(**inputs):
        calls.append(inputs['dicttrue'])
        positive = inputs['dicttrue']['typemodl'] == 'PlanetarySystem'
        return {'boolcalclspe': False, 'boolsrchboxsperi': False,
                'boolsrchoutlperi': True,
                'dictoutlperi': {'boolposi': positive,
                                 'minmfrddtimeoutlsort': [0.05 if positive else 0.5]}}

    monkeypatch.setattr(main, '_get_miletos', lambda: SimpleNamespace(init=analyze))
    main.mile_work(state, 0)
    assert [item['typemodl'] for item in calls] == ['PlanetarySystem', 'StellarBinary']
    assert state.booltargproc.tolist() == [True, True]
    assert state.boolpositarg[0].tolist() == [True, False]
    assert state.boolpositarg[1].tolist() == [False, True]