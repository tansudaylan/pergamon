"""Candidate-vetting workflows and summaries absorbed from Sardis."""

import datetime
import os

import numpy as np
import tdpy
from tdpy.verbosity import print

from .paths import get_data_path, get_repository_path, get_visuals_path
from . import vetting_workflow
from .vetting_workflow import (
    mile_work,
    retr_dictderi_effe,
    retr_pathvetting,
    synthetic_target_truth,
    validated_classification,
)


def vetting_paths(pathbase=None, strgcnfg=None):
    """Return normalized workflow paths rooted in PERGAMON_PATH or an explicit base.

    The historical pathbasesard key is retained for existing callers. Generic
    repository path helpers are supplied by pergamon.paths.
    """

    paths = retr_pathvetting(pathbase=pathbase)
    paths['pathbasesard'] = paths['pathbasepergamon']
    if strgcnfg is not None:
        config_paths = retr_pathvetting(pathbase=paths['pathbase'], strgextn=strgcnfg)
        configuration_path = config_paths['pathpopl']
        paths['pathcnfg'] = configuration_path
        paths['pathsimu'] = tdpy.ensr_path(os.path.join(configuration_path, 'simulation'))
        paths['pathobsd'] = tdpy.ensr_path(os.path.join(configuration_path, 'observed'))
    return paths


retr_pathsard = vetting_paths


def run_candidate_vetting(
    typesyst='PlanetarySystem',
    dictinpt=None,
    strgcnfg=None,
    pathbase=None,
    *,
    dicttroiinpt=None,
):
    """Run candidate classification without inferring completeness or occurrence."""

    print('Candidate vetting initialized at %s...' % datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))
    if dictinpt is not None and dicttroiinpt is not None:
        raise ValueError('Supply only one of dictinpt and dicttroiinpt.')
    if dicttroiinpt is not None:
        dictinpt = dicttroiinpt
    configuration = typesyst if strgcnfg is None else strgcnfg
    paths = vetting_paths(pathbase=pathbase, strgcnfg=configuration)
    if dictinpt is None:
        dictinpt = {
            'pathbase': paths['pathsimu'],
            'typesyst': typesyst,
            'listlablinst': [['TESS'], []],
            'liststrgtypedata': [['simutargsynt'], []],
        }
    result = vetting_workflow.init(**dictinpt)
    if hasattr(result, '__dict__'):
        return result.__dict__
    if isinstance(result, dict):
        return result
    return {'result': result}


init = run_candidate_vetting


def summarize_vetting_result(result: dict) -> dict[str, object]:
    """Summarize classifier performance against simulated truth labels."""

    true_type_names = list(result["listnameclastruetype"])
    if (not true_type_names or any(not isinstance(name, str) or not name for name in true_type_names)
            or len(set(true_type_names)) != len(true_type_names)):
        raise ValueError("Truth class names must be distinct nonempty strings.")
    target_type = result["typesyst"]
    if target_type not in true_type_names:
        raise ValueError(f"Target type {target_type!r} is absent from truth labels.")
    raw_indices = result["indxclastruetypetarg"]
    if np.ma.is_masked(raw_indices):
        raise ValueError("Truth class indices must be complete, finite integer labels.")
    true_type_index = np.asarray(raw_indices)
    if (true_type_index.ndim != 1 or true_type_index.dtype.kind not in 'biuf'
            or not np.all(np.isfinite(true_type_index))
            or np.any(true_type_index != np.floor(true_type_index))
            or np.any((true_type_index < 0) | (true_type_index >= len(true_type_names)))):
        raise ValueError("Truth class indices must be valid finite integer labels.")
    positive_truth = true_type_index == true_type_names.index(target_type)
    classifier_labels = list(result["listlablclasdisp"])
    predictions = []
    for values in result["boolpositarg"]:
        if np.ma.is_masked(values):
            raise ValueError("Prediction arrays must contain complete binary decisions.")
        array = np.asarray(values)
        if (array.dtype.kind not in 'biuf' or not np.all(np.isfinite(array))
                or not np.all((array == 0) | (array == 1))):
            raise ValueError("Prediction arrays must contain finite binary decisions.")
        predictions.append(array.astype(bool))
    if len(predictions) != len(classifier_labels):
        raise ValueError("Classifier labels and prediction arrays must have equal length.")
    if any(values.shape != positive_truth.shape for values in predictions):
        raise ValueError("Each prediction array must match the truth-label shape.")
    confusion_matrices = np.empty((len(predictions), 2, 2), dtype=int)
    precision = np.zeros(len(predictions))
    recall = np.zeros(len(predictions))
    for index, predicted_positive in enumerate(predictions):
        true_positive = np.count_nonzero(predicted_positive & positive_truth)
        false_positive = np.count_nonzero(predicted_positive & ~positive_truth)
        false_negative = np.count_nonzero(~predicted_positive & positive_truth)
        true_negative = np.count_nonzero(~predicted_positive & ~positive_truth)
        confusion_matrices[index] = [[true_negative, false_positive], [false_negative, true_positive]]
        if true_positive + false_positive > 0:
            precision[index] = true_positive / (true_positive + false_positive)
        if true_positive + false_negative > 0:
            recall[index] = true_positive / (true_positive + false_negative)
    return {
        "classifier_labels": classifier_labels,
        "confusion_matrices": confusion_matrices,
        "precision": precision,
        "recall": recall,
        "positive_truth": positive_truth,
    }