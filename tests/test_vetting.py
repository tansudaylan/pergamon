from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import pergamon
from pergamon import vetting


def test_summarize_vetting_result_preserves_classifier_metrics():
    result = {
        'typesyst': 'PlanetarySystem',
        'listnameclastruetype': ['PlanetarySystem', 'StellarBinary'],
        'indxclastruetypetarg': np.array([0, 0, 1, 1]),
        'listlablclasdisp': ['Selective', 'Inclusive'],
        'boolpositarg': [np.array([True, False, True, False]), np.array([True, True, True, True])],
    }

    summary = vetting.summarize_vetting_result(result)

    np.testing.assert_array_equal(summary['confusion_matrices'], [[[1, 1], [1, 1]], [[0, 2], [0, 2]]])
    np.testing.assert_allclose(summary['precision'], [0.5, 0.5])
    np.testing.assert_allclose(summary['recall'], [0.5, 1.0])


@pytest.mark.parametrize('indices', [
    [0.5, 1], [np.nan, 1], [-1, 1], [0, 2], [[0, 1]], ['0', '1'],
    np.ma.array([0, 1], mask=[True, False]),
])
def test_vetting_summary_rejects_invalid_truth_indices(indices):
    with pytest.raises(ValueError, match='Truth class indices'):
        vetting.summarize_vetting_result({
            'typesyst': 'PlanetarySystem', 'listnameclastruetype': ['PlanetarySystem', 'StellarBinary'],
            'indxclastruetypetarg': indices, 'listlablclasdisp': ['Classifier'],
            'boolpositarg': [[True, False]],
        })


@pytest.mark.parametrize('predictions', [
    [np.nan, 0], [np.inf, 0], [0.5, 0], [-1, 0], [2, 0], ['False', 'True'],
    np.ma.array([True, False], mask=[False, True]),
])
def test_vetting_summary_rejects_invalid_predictions(predictions):
    with pytest.raises(ValueError, match='Prediction arrays'):
        vetting.summarize_vetting_result({
            'typesyst': 'PlanetarySystem', 'listnameclastruetype': ['PlanetarySystem', 'StellarBinary'],
            'indxclastruetypetarg': [0, 1], 'listlablclasdisp': ['Classifier'],
            'boolpositarg': [predictions],
        })


@pytest.mark.parametrize('predictions', [
    np.array([True, False, True, False]), np.array([1, 0, 1, 0]),
    np.array([1.0, 0.0, 1.0, 0.0]),
    np.ma.array([True, False, True, False], mask=False),
])
def test_vetting_summary_preserves_valid_binary_encodings(predictions):
    original = predictions.copy()
    summary = vetting.summarize_vetting_result({
        'typesyst': 'PlanetarySystem', 'listnameclastruetype': ['PlanetarySystem', 'StellarBinary'],
        'indxclastruetypetarg': np.array([0.0, 0.0, 1.0, 1.0]),
        'listlablclasdisp': ['Classifier'], 'boolpositarg': [predictions],
    })
    np.testing.assert_array_equal(summary['confusion_matrices'], [[[1, 1], [1, 1]]])
    np.testing.assert_array_equal(summary['precision'], [0.5])
    np.testing.assert_array_equal(summary['recall'], [0.5])
    np.testing.assert_array_equal(predictions, original)
    np.testing.assert_array_equal(np.ma.getmaskarray(predictions), np.ma.getmaskarray(original))


@pytest.mark.parametrize('names', [[], ['', 'Other'], ['Same', 'Same'], [1, 'Other']])
def test_vetting_summary_rejects_ambiguous_class_names(names):
    with pytest.raises(ValueError, match='Truth class names'):
        vetting.summarize_vetting_result({'listnameclastruetype': names})


def test_summarize_vetting_result_rejects_unknown_target_type():
    with pytest.raises(ValueError, match='absent from truth labels'):
        vetting.summarize_vetting_result({
            'typesyst': 'PlanetarySystem', 'listnameclastruetype': ['StellarBinary'],
            'indxclastruetypetarg': np.array([0]), 'listlablclasdisp': [], 'boolpositarg': [],
        })


def test_summary_retains_zero_denominator_convention():
    summary = vetting.summarize_vetting_result({
        'typesyst': 'PlanetarySystem', 'listnameclastruetype': ['PlanetarySystem', 'StellarBinary'],
        'indxclastruetypetarg': np.array([1, 1]), 'listlablclasdisp': ['Reject all'],
        'boolpositarg': [np.array([False, False])],
    })
    np.testing.assert_array_equal(summary['precision'], [0.0])
    np.testing.assert_array_equal(summary['recall'], [0.0])


def test_summary_retains_prediction_shape_validation():
    with pytest.raises(ValueError, match='truth-label shape'):
        vetting.summarize_vetting_result({
            'typesyst': 'PlanetarySystem', 'listnameclastruetype': ['PlanetarySystem'],
            'indxclastruetypetarg': np.array([0, 0]), 'listlablclasdisp': ['Classifier'],
            'boolpositarg': [np.array([True])],
        })


def test_default_workflow_preserves_analysis_arguments_and_output(monkeypatch, tmp_path):
    monkeypatch.setenv('PERGAMON_PATH', str(tmp_path))
    captured = {}

    def run(**arguments):
        captured.update(arguments)
        return arguments

    monkeypatch.setattr(vetting.vetting_workflow, 'init', run)
    result = pergamon.run_candidate_vetting(pathbase=str(tmp_path), strgcnfg='TransitVet')

    assert result == captured
    assert Path(captured['pathbase']) == tmp_path / 'TransitVet' / 'simulation'
    assert captured['typesyst'] == 'PlanetarySystem'
    assert captured['listlablinst'] == [['TESS'], []]
    assert captured['liststrgtypedata'] == [['simutargsynt'], []]


@pytest.mark.parametrize('output', [{'value': 1}, SimpleNamespace(value=1), 3])
def test_explicit_inputs_and_output_conventions_are_preserved(monkeypatch, tmp_path, output):
    monkeypatch.setenv('PERGAMON_PATH', str(tmp_path))
    arguments = {'typesyst': 'StellarBinary', 'custom': [1, 2]}
    captured = {}

    def run(**supplied):
        captured.update(supplied)
        return output

    monkeypatch.setattr(vetting.vetting_workflow, 'init', run)
    result = vetting.init(dictinpt=arguments)

    assert captured == arguments
    assert result == ({'result': 3} if output == 3 else {'value': 1})


def test_vetting_paths_preserve_configuration_directories_and_legacy_key(monkeypatch, tmp_path):
    monkeypatch.setenv('PERGAMON_PATH', str(tmp_path / 'repository'))
    paths = vetting.retr_pathsard(pathbase=tmp_path / 'output', strgcnfg='TransitVet')

    assert Path(paths['pathbasesard']) == tmp_path / 'repository'
    assert Path(paths['pathcnfg']) == tmp_path / 'output' / 'TransitVet'
    for name in ('pathbase', 'pathdatapipe', 'pathvisupipe', 'pathcnfg', 'pathsimu', 'pathobsd'):
        assert Path(paths[name]).is_dir()


def test_vetting_does_not_replace_existing_population_init():
    assert pergamon.init is pergamon.main.init
    assert pergamon.init is not pergamon.run_candidate_vetting


def test_vetting_module_retains_all_generic_path_helpers(monkeypatch, tmp_path):
    monkeypatch.setenv('PERGAMON_PATH', str(tmp_path))
    assert vetting.get_repository_path() == tmp_path
    assert vetting.get_data_path() == tmp_path / 'data'
    assert vetting.get_visuals_path() == tmp_path / 'visuals'


def test_vetting_paths_require_the_new_repository_environment(monkeypatch):
    monkeypatch.delenv('PERGAMON_PATH', raising=False)
    with pytest.raises(EnvironmentError, match='PERGAMON_PATH'):
        vetting.get_repository_path()


def test_full_population_workflow_with_explicit_deterministic_analysis_responses(monkeypatch, tmp_path):
    import pergamon.vetting_workflow as backend

    monkeypatch.setenv('PERGAMON_PATH', str(tmp_path))
    calls = []

    def analyze(**arguments):
        calls.append(arguments['dicttrue']['typemodl'])
        selected = calls[-1] == 'PlanetarySystem'
        return {
            'boolcalclspe': False, 'boolsrchboxsperi': False, 'boolsrchoutlperi': True,
            'dictoutlperi': {'boolposi': selected, 'minmfrddtimeoutlsort': [0.05 if selected else 0.5]},
        }

    monkeypatch.setattr(backend._get_miletos(), 'init', analyze)
    configuration = {
        'pathbase': str(tmp_path / 'run'), 'typesyst': 'PlanetarySystem',
        'listlablinst': [['TESS'], []], 'liststrgtypedata': [['simutargsynt'], []],
        'boolplot': False, 'typeverb': 0, 'dictmileinptglob': {},
    }
    random_state = np.random.get_state()
    try:
        result = vetting.run_candidate_vetting(dictinpt=configuration)
    finally:
        np.random.set_state(random_state)
    summary = vetting.summarize_vetting_result(result)

    assert len(calls) == len(result['indxclastruetypetarg']) == 100
    assert set(calls) == {'PlanetarySystem', 'StellarBinary'}
    np.testing.assert_array_equal(result['boolpositarg'][0], np.asarray(calls) == 'PlanetarySystem')
    np.testing.assert_array_equal(result['boolpositarg'][1], np.asarray(calls) != 'PlanetarySystem')
    actual_truth = np.asarray(calls) == 'PlanetarySystem'
    np.testing.assert_array_equal(summary['positive_truth'], actual_truth)
    assert summary['confusion_matrices'].shape == (2, 2, 2)
    assert np.all(summary['confusion_matrices'].sum(axis=(1, 2)) == 100)
    positive_count = np.count_nonzero(actual_truth)
    negative_count = actual_truth.size - positive_count
    np.testing.assert_array_equal(summary['confusion_matrices'][0],
                                  [[negative_count, 0], [0, positive_count]])