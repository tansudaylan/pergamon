"""Occurrence fractions from independently surveyed targets and known completeness."""

from collections.abc import Callable, Mapping

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import xlogy


def _validate_survey(detected, detection_efficiency, false_positive_probability=0.0):
    detected = np.asarray(detected)
    efficiency = np.asarray(detection_efficiency, dtype=float)
    false_positive_probability = np.broadcast_to(
        np.asarray(false_positive_probability, dtype=float), efficiency.shape
    )
    if detected.ndim != 1 or detected.size == 0 or detected.shape != efficiency.shape:
        raise ValueError("detections and efficiencies must be nonempty matching one-dimensional arrays")
    if not np.isin(detected, (0, 1)).all():
        raise ValueError("detections must contain only zeros and ones")
    if not np.isfinite(efficiency).all() or np.any((efficiency < 0) | (efficiency > 1)):
        raise ValueError("detection efficiencies must be finite probabilities")
    if (not np.isfinite(false_positive_probability).all()
            or np.any((false_positive_probability < 0) | (false_positive_probability > 1))):
        raise ValueError("false-positive probabilities must be finite probabilities")
    return detected, efficiency, false_positive_probability


def log_likelihood_occurrence_rate(
    occurrence_fraction, detected, detection_efficiency, false_positive_probability=0.0
):
    """Bernoulli log likelihood for at most one event per independently surveyed target.

    Each target has detection probability ``1 - (1 - f * efficiency) * (1 - alpha)``,
    where ``alpha`` is the supplied per-target false-positive probability.
    """

    detected, efficiency, false_positive_probability = _validate_survey(
        detected, detection_efficiency, false_positive_probability
    )
    occurrence_fraction = np.asarray(occurrence_fraction, dtype=float)
    if occurrence_fraction.ndim == 0:
        if not np.isfinite(occurrence_fraction) or not 0 <= occurrence_fraction <= 1:
            return -np.inf
        probability = 1.0 - (1.0 - occurrence_fraction * efficiency) * (1.0 - false_positive_probability)
        return float(np.sum(xlogy(detected, probability) + xlogy(1 - detected, 1 - probability)))
    if occurrence_fraction.ndim != 1 or not np.isfinite(occurrence_fraction).all() or np.any(
        (occurrence_fraction < 0) | (occurrence_fraction > 1)
    ):
        raise ValueError("occurrence_fraction must be a scalar or a one-dimensional probability grid")
    probability = 1.0 - (
        1.0 - occurrence_fraction[:, None] * efficiency[None, :]
    ) * (1.0 - false_positive_probability[None, :])
    return np.sum(
        xlogy(detected[None, :], probability)
        + xlogy(1 - detected[None, :], 1.0 - probability),
        axis=1,
    )


def estimate_occurrence_rate(
    detected, detection_efficiency, false_positive_probability=0.0
):
    """Estimate occurrence from detections and caller-supplied efficiencies.

    The efficiencies may come from population recovery simulations that use
    precomputed Miletos/Ephesos signatures. This function consumes those values;
    it does not predict photometric signatures or detection completeness.
    """

    detected, efficiency, false_positive_probability = _validate_survey(
        detected, detection_efficiency, false_positive_probability
    )
    detection_probability_at_zero_occurrence = false_positive_probability
    if np.any((detected == 1) & (efficiency == 0)
              & (detection_probability_at_zero_occurrence == 0)):
        raise ValueError(
            "a detection cannot have zero detection efficiency and zero false-positive probability"
        )
    if not np.any(efficiency > 0):
        raise ValueError("at least one surveyed target must have positive detection efficiency")

    def log_likelihood(rate):
        return log_likelihood_occurrence_rate(
            rate, detected, efficiency, false_positive_probability
        )

    result = minimize_scalar(lambda rate: -log_likelihood(rate), bounds=(0, 1), method="bounded")
    return float(max((0.0, result.x, 1.0), key=log_likelihood))


def estimate_compact_object_occurrence_rate(
    detected,
    photometric_signatures: Mapping[str, np.ndarray],
    recovery_model: Callable[[Mapping[str, np.ndarray]], np.ndarray],
) -> float:
    """Estimate occurrence using Miletos amplitudes and a supplied recovery model.

    The recovery model receives per-target `amplbeammodl`, `ampellimodl`, and
    `amplslenmodl` arrays and returns detection efficiencies. Pergamon reads and
    forwards these predictions but does not calculate amplitudes or completeness.
    """

    detected = np.asarray(detected)
    if detected.ndim != 1 or detected.size == 0:
        raise ValueError("detections must be a nonempty one-dimensional array")

    names = ("amplbeammodl", "ampellimodl", "amplslenmodl")
    missing = set(names).difference(photometric_signatures)
    if missing:
        raise ValueError(f"missing compact-object amplitudes: {sorted(missing)}")

    signatures = {
        name: np.asarray(photometric_signatures[name], dtype=float)
        for name in names
    }
    if any(
        values.shape != detected.shape or not np.isfinite(values).all()
        for values in signatures.values()
    ):
        raise ValueError("compact-object amplitudes must be finite per-target arrays")

    detection_efficiency = recovery_model(signatures)
    return estimate_occurrence_rate(detected, detection_efficiency)