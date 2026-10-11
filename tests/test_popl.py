import numpy as np
import pytest
from types import SimpleNamespace
import importlib

import pergamon
from pergamon.main import retr_llik_corr


def test_correlation_likelihood_normalizes_uncertain_data():
    state = SimpleNamespace(
        tempfrst=np.array([0., 1.]), tempseco=np.zeros(2),
        tempfrststdv=np.ones(2), tempsecostdv=np.ones(2),
    )
    np.testing.assert_allclose(retr_llik_corr(np.array([np.pi / 2., 0.]), state),
                               -np.log(2. * np.pi))
    assert retr_llik_corr(np.array([np.pi / 2., 1.]), state) < -np.log(2. * np.pi)


def test_optional_correlation_search_calls_pcat_backed_sampler(monkeypatch, tmp_path):
    monkeypatch.setenv('PERGAMON_PATH', str(tmp_path))
    module = importlib.import_module('pergamon.main')
    calls = []

    def fake_sampler(*args, **kwargs):
        calls.append(args)
        return {'angle': np.full(20, np.pi / 2.), 'intercept': np.zeros(20)}

    monkeypatch.setattr(module, 'sample_posterior', fake_sampler)
    pergamon.init(
        typeanls='defa', dictpopl={'pop': {'radistar': np.array([1., 1.2, 1.4]),
                                           'massstar': np.array([1., 1.1, 1.3])}},
        boolsrchcorr=True, boolmakeplot=False, booldiag=False,
    )
    assert calls


def test_init_accepts_flattened_array_population():
    dictpopl = {
        'pop': {
            'radistar': np.array([1.0, 1.2]),
            'massstar': np.array([1.0, 1.1]),
            'amplbeammodl': np.array([1.5, 2.5]),
            'ampellimodl': np.array([3.5, 4.5]),
            'amplslenmodl': np.array([5.5, 6.5]),
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
    np.testing.assert_allclose(result['pop']['amplbeammodl'][0], np.array([1.5, 2.5]))
    np.testing.assert_allclose(result['pop']['ampellimodl'][0], np.array([3.5, 4.5]))
    np.testing.assert_allclose(result['pop']['amplslenmodl'][0], np.array([5.5, 6.5]))


def test_init_uses_repository_runtime_paths(monkeypatch, tmp_path):
    monkeypatch.setenv('PERGAMON_PATH', str(tmp_path))

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
    assert (tmp_path / 'data' / 'defa').is_dir()
    assert (tmp_path / 'visuals' / 'defa').is_dir()


def test_retr_subp_accepts_python_list_indices():
    dictpopl = {'pop': {'radistar': [np.array([1.0, 2.0, 3.0]), '']}}
    dictnumbsamp = {}
    dictindxsamp = {}

    pergamon.retr_subp(dictpopl, 'pop', 'small', [0, 2], dictnumbsamp, dictindxsamp)

    np.testing.assert_array_equal(dictpopl['small']['radistar'][0], np.array([1.0, 3.0]))
    np.testing.assert_array_equal(dictindxsamp['pop']['small'], np.array([0, 2]))
    assert dictnumbsamp['small'] == 2


def test_classified_population_partition_preserves_target_ids():
    relevant = [7, 2]
    groups = pergamon.partition_classified_population(
        relevant=relevant, irrelevant=[5, 9], positive=[7, 5], negative=[2, 9],
    )
    expected = {
        're': [7, 2], 'ir': [5, 9], 'po': [7, 5], 'ne': [2, 9],
        'trpo': [7], 'trne': [9], 'flpo': [5], 'flne': [2],
    }
    assert set(groups) == set(expected)
    assert groups['re'] is relevant
    for name, indices in expected.items():
        np.testing.assert_array_equal(groups[name], indices)


def test_occurrence_rate_accounts_for_target_detection_efficiency():
    detections = [1, 0, 0]
    efficiencies = [1.0, 0.5, 0.5]
    assert pergamon.estimate_occurrence_rate(detections, efficiencies) == pytest.approx(2 / 3, abs=1e-5)
    assert pergamon.log_likelihood_occurrence_rate(0.5, [1, 0], [1.0, 0.5]) == pytest.approx(
        np.log(0.5 * 0.75)
    )
    assert pergamon.estimate_occurrence_rate([0, 0], [1.0, 0.5]) == 0.0
    assert pergamon.estimate_occurrence_rate([1, 1], [1.0, 0.5]) == 1.0


def test_occurrence_rate_accounts_for_inverted_control_false_alarms():
    detections = [1, 0, 0]
    efficiencies = [0.8, 0.5, 0.4]
    false_alarm_probabilities = [0.1, 0.1, 0.1]

    likelihood_with_control = pergamon.log_likelihood_occurrence_rate(
        0.5, detections, efficiencies, false_alarm_probabilities
    )
    likelihood_without_control = pergamon.log_likelihood_occurrence_rate(
        0.5, detections, efficiencies
    )

    assert np.isfinite(likelihood_with_control)
    assert likelihood_with_control != likelihood_without_control
    assert pergamon.estimate_occurrence_rate(
        detections, efficiencies, false_alarm_probabilities
    ) < pergamon.estimate_occurrence_rate(detections, efficiencies)
    rate_grid = np.array([0.0, 0.5, 1.0])
    grid_likelihood = pergamon.log_likelihood_occurrence_rate(
        rate_grid, detections, efficiencies, false_alarm_probabilities
    )
    scalar_likelihood = np.array([
        pergamon.log_likelihood_occurrence_rate(
            rate, detections, efficiencies, false_alarm_probabilities
        )
        for rate in rate_grid
    ])
    np.testing.assert_allclose(grid_likelihood, scalar_likelihood)


def test_compact_object_occurrence_rate_reads_precomputed_miletos_signatures():
    signatures = {
        'amplbeammodl': np.array([1.0, 2.0, 3.0]),
        'ampellimodl': np.array([4.0, 5.0, 6.0]),
        'amplslenmodl': np.array([7.0, 8.0, 9.0]),
    }
    received = {}

    def recovery_model(population_signatures):
        received.update(population_signatures)
        return np.array([1.0, 0.5, 0.5])

    occurrence_rate = pergamon.estimate_compact_object_occurrence_rate(
        [1, 0, 0], signatures, recovery_model
    )

    assert occurrence_rate == pytest.approx(2 / 3, abs=1e-5)
    for name, values in signatures.items():
        np.testing.assert_array_equal(received[name], values)


def test_compact_object_occurrence_rate_rejects_missing_or_misaligned_signatures():
    signatures = {
        'amplbeammodl': np.array([1.0, 2.0]),
        'ampellimodl': np.array([3.0, 4.0]),
    }
    with pytest.raises(ValueError, match="missing compact-object amplitudes"):
        pergamon.estimate_compact_object_occurrence_rate(
            [1, 0], signatures, lambda values: np.ones(2)
        )

    signatures['amplslenmodl'] = np.array([5.0])
    with pytest.raises(ValueError, match="finite per-target arrays"):
        pergamon.estimate_compact_object_occurrence_rate(
            [1, 0], signatures, lambda values: np.ones(2)
        )


def test_occurrence_rate_rejects_uninformative_or_invalid_survey_data():
    assert pergamon.log_likelihood_occurrence_rate(1.1, [0], [1.0]) == -np.inf
    with pytest.raises(ValueError, match="probabilities"):
        pergamon.estimate_occurrence_rate([0], [1.1])
    with pytest.raises(ValueError, match="positive detection efficiency"):
        pergamon.estimate_occurrence_rate([0], [0.0])
    with pytest.raises(ValueError, match="zero detection efficiency"):
        pergamon.estimate_occurrence_rate([1], [0.0])
    with pytest.raises(ValueError, match="zeros and ones"):
        pergamon.estimate_occurrence_rate([2], [1.0])


def test_compact_object_model_functions_are_not_population_package_exports():
    assert not hasattr(pergamon, 'compute_photometric_signatures')
    assert not hasattr(pergamon, 'derive_compact_object_features')
