import numpy as np
import ephesos
import miletos

from pergamon.vetting_workflow import retr_dictderi_effe


def test_miletos_compact_object_signatures_match_ephesos():
    period_days = np.array([0.3, 30.0])[:, None]  # [day]
    companion_mass_solar = np.array([5.0, 180.0])[None, :]  # [M_Sun]

    signatures = miletos.compute_photometric_signatures(
        period_days,
        companion_mass_solar,
    )
    expected = ephesos.evaluate_compact_object_signatures(
        period_days, companion_mass_solar
    )
    for name in expected:
        np.testing.assert_array_equal(signatures[name], expected[name])


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
    signatures = ephesos.evaluate_compact_object_signatures(10.0, 5.0)
    np.testing.assert_allclose(derived["amplbeammodl"], [signatures["beaming"]])
    np.testing.assert_allclose(derived["ampellimodl"], [signatures["ellipsoidal"]])
    np.testing.assert_allclose(derived["amplslenmodl"], np.array([3.015272]), rtol=1e-6)
    assert derived["duratrantotlmodl"][0] > 0.0
