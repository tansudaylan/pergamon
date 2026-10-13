import io
from urllib.parse import parse_qs, urlparse

import matplotlib.image as mpimg
import numpy as np
import pytest

from pergamon import emission_spectroscopy as catalog


def test_fetch_esm_catalog_uses_archive_values_and_filters_invalid_metrics(monkeypatch):
    csv_data = (
        "pl_name,hostname,pl_eqt,st_teff,pl_rade,st_rad,sy_kmag\n"
        "Published planet,Published host,800,5000,2.5,0.8,8.0\n"
        "Invalid planet,Other host,900,5100,3.0,0.9,9.0\n"
        "Zero metric,Third host,1000,5200,3.5,1.0,10.0\n"
    )
    requests = []
    metric_inputs = []

    def fetch(url, timeout):
        requests.append((url, timeout))
        return io.BytesIO(csv_data.encode())

    def calculate_metric(*args):
        metric_inputs.extend(args)
        return np.array([4.2, np.nan, 0.0])

    monkeypatch.setattr(catalog, "urlopen", fetch)
    monkeypatch.setattr(catalog.tdpy, "retr_esmm", calculate_metric)

    records = catalog.fetch_esm_catalog(timeout_seconds=17)

    assert len(records) == 1
    assert records[0] == {
        "planet_name": "Published planet",
        "host_name": "Published host",
        "equilibrium_temperature_kelvin": 800.0,
        "k_magnitude": 8.0,
        "esm": 4.2,
    }
    assert requests[0][1] == 17
    assert parse_qs(urlparse(requests[0][0]).query) == {
        "query": [catalog.ARCHIVE_QUERY], "format": ["csv"],
    }
    for actual, expected in zip(metric_inputs, (
        [800.0, 900.0, 1000.0], [5000.0, 5100.0, 5200.0],
        [2.5, 3.0, 3.5], [0.8, 0.9, 1.0], [8.0, 9.0, 10.0],
    )):
        np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("function_name", [
    "generate_esm_catalog_figure", "generate_esm_distribution_figure",
])
def test_generate_esm_figures(tmp_path, monkeypatch, function_name):
    requested_timeouts = []

    def fetch(timeout_seconds):
        requested_timeouts.append(timeout_seconds)
        return [{
            "planet_name": "Published planet", "host_name": "Published host",
            "equilibrium_temperature_kelvin": 800.0, "k_magnitude": 8.0, "esm": 4.2,
        }]

    monkeypatch.setattr(catalog, "fetch_esm_catalog", fetch)
    output_path = tmp_path / "figures" / "archive_esm.png"

    generate = getattr(catalog, function_name)
    assert generate(output_path, timeout_seconds=17) == output_path
    assert requested_timeouts == [17]
    image = mpimg.imread(output_path)
    assert image.shape[0] > 100 and image.shape[1] > 100
    assert image[..., :3].min() < 0.8