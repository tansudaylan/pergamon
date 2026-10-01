"""Photometric features of compact-object companions across target populations."""

import numpy as np

import chalcedon
import nicomedia


def compute_photometric_signatures(
    period_days,
    companion_mass_solar,
    stellar_radius_solar=1.0,
    stellar_mass_solar=1.0,
    stellar_density_cgs=1.41,
) -> dict[str, np.ndarray]:
    """Return beaming, ellipsoidal, and self-lensing amplitudes in ppt."""

    period_days = np.asarray(period_days, dtype=float)
    companion_mass_solar = np.asarray(companion_mass_solar, dtype=float)
    if np.any(period_days <= 0.0):
        raise ValueError("Orbital periods must be positive.")
    if np.any(companion_mass_solar <= 0.0):
        raise ValueError("Companion masses must be positive.")
    if stellar_radius_solar <= 0.0 or stellar_mass_solar <= 0.0:
        raise ValueError("Stellar radius and mass must be positive.")
    if stellar_density_cgs <= 0.0:
        raise ValueError("Stellar density must be positive.")

    return {
        "beaming": np.asarray(
            nicomedia.retr_deptbeam(period_days, stellar_mass_solar, companion_mass_solar)
        ),
        "ellipsoidal": np.asarray(
            nicomedia.retr_deptelli(
                period_days, stellar_density_cgs, stellar_mass_solar, companion_mass_solar
            )
        ),
        "self_lensing": np.asarray(
            chalcedon.retr_amplslen(
                period_days, stellar_radius_solar, companion_mass_solar, stellar_mass_solar
            )
        ),
    }


def derive_compact_object_features(stellar_radius_solar, period_days,
                                   companion_mass_solar, stellar_mass_solar):
    """Return self-lensing, transit duration, orbit size, and Schwarzschild radius."""

    signatures = compute_photometric_signatures(
        period_days,
        companion_mass_solar,
        stellar_radius_solar=stellar_radius_solar,
        stellar_mass_solar=stellar_mass_solar,
    )
    semimajor_axis_solar = nicomedia.retr_smaxkepl(
        period_days, stellar_mass_solar + companion_mass_solar
    ) * 215.0  # [R_Sun]
    duration_hours = nicomedia.retr_duratrantotl(
        np.atleast_1d(period_days),
        np.atleast_1d(stellar_radius_solar / semimajor_axis_solar),
        np.zeros(1),
    )  # [hour]
    schwarzschild_radius_solar = 4.24e-6 * companion_mass_solar  # [R_Sun]
    return {
        'amplslenmodl': np.atleast_1d(signatures['self_lensing']),
        'duratrantotlmodl': duration_hours,
        'smaxmodl': np.atleast_1d(semimajor_axis_solar),
        'radischw': np.atleast_1d(schwarzschild_radius_solar),
    }