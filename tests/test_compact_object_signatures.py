import numpy as np
import pytest

from pergamon.signatures import compute_photometric_signatures
from pergamon.vetting_workflow import retr_dictderi_effe


def test_compute_photometric_signatures_matches_reference_values():
    period_days = np.array([0.3, 30.0])[:, None]  # [day]
    companion_mass_solar = np.array([5.0, 180.0])[None, :]  # [M_Sun]

    signatures = compute_photometric_signatures(
        period_days,
        companion_mass_solar,
    )

    np.testing.assert_allclose(
        signatures["beaming"],
        np.array([[6.333641, 23.529050], [1.364542, 5.069180]]),
        rtol=1e-6,
    )
    np.testing.assert_allclose(
        signatures["ellipsoidal"],
        np.array([[124.113475, 148.113266], [0.012411, 0.014811]]),
        rtol=5e-5,
    )
    np.testing.assert_allclose(
        signatures["self_lensing"],
        np.array([[0.291121, 32.625129], [6.272018, 702.887070]]),
        rtol=1e-6,
    )


def test_retr_dictderi_effe_returns_finite_edge_on_quantities():
    derived, variable_derived = retr_dictderi_effe(
        [1.0, 10.0, 5.0, 1.0],
        None,
    )

    assert variable_derived is None
    assert set(derived) == {
        "amplbeammodl",
        "ampellimodl",
        "amplslenmodl",
        "duratrantotlmodl",
        "smaxmodl",
        "radischw",
    }
    assert all(np.all(np.isfinite(values)) for values in derived.values())
    signatures = compute_photometric_signatures(10.0, 5.0)
    np.testing.assert_allclose(derived["amplbeammodl"], [signatures["beaming"]])
    np.testing.assert_allclose(derived["ampellimodl"], [signatures["ellipsoidal"]])
    np.testing.assert_allclose(derived["amplslenmodl"], np.array([3.015272]), rtol=1e-6)
    assert derived["duratrantotlmodl"][0] > 0.0


@pytest.mark.parametrize(
    "period_days, companion_mass_solar",
    [(0.0, 5.0), (1.0, -5.0)],
)
def test_compute_photometric_signatures_rejects_nonphysical_inputs(
    period_days,
    companion_mass_solar,
):
    with pytest.raises(ValueError, match="must be positive"):
        compute_photometric_signatures(period_days, companion_mass_solar)