"""Occurrence fractions from independently surveyed targets and known completeness."""

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import xlogy


def _validate_survey(detected, detection_efficiency):
    detected = np.asarray(detected)
    efficiency = np.asarray(detection_efficiency, dtype=float)
    if detected.ndim != 1 or detected.size == 0 or detected.shape != efficiency.shape:
        raise ValueError("detections and efficiencies must be nonempty matching one-dimensional arrays")
    if not np.isin(detected, (0, 1)).all():
        raise ValueError("detections must contain only zeros and ones")
    if not np.isfinite(efficiency).all() or np.any((efficiency < 0) | (efficiency > 1)):
        raise ValueError("detection efficiencies must be finite probabilities")
    return detected, efficiency


def log_likelihood_occurrence_rate(occurrence_fraction, detected, detection_efficiency):
    """Bernoulli log likelihood for at most one event per independently surveyed target.

    Each target has detection probability occurrence_fraction * detection_efficiency.
    False positives are assumed absent.
    """

    detected, efficiency = _validate_survey(detected, detection_efficiency)
    if not np.isfinite(occurrence_fraction) or not 0 <= occurrence_fraction <= 1:
        return -np.inf
    probability = occurrence_fraction * efficiency
    return float(np.sum(xlogy(detected, probability) + xlogy(1 - detected, 1 - probability)))


def estimate_occurrence_rate(detected, detection_efficiency):
    """Return the maximum-likelihood occurrence fraction between zero and one."""

    detected, efficiency = _validate_survey(detected, detection_efficiency)
    if np.any((detected == 1) & (efficiency == 0)):
        raise ValueError("a detection cannot have zero detection efficiency")
    if not np.any(efficiency > 0):
        raise ValueError("at least one surveyed target must have positive detection efficiency")

    def log_likelihood(rate):
        return log_likelihood_occurrence_rate(rate, detected, efficiency)

    result = minimize_scalar(lambda rate: -log_likelihood(rate), bounds=(0, 1), method="bounded")
    return float(max((0.0, result.x, 1.0), key=log_likelihood))