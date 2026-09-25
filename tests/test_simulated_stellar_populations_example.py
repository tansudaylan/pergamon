import importlib.util
from pathlib import Path

import matplotlib.image as mpimg
import numpy as np


EXAMPLE_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "examples"
    / "simulated_stellar_populations.py"
)
SPECIFICATION = importlib.util.spec_from_file_location(
    "pergamon_simulated_stellar_populations", EXAMPLE_SCRIPT
)
example = importlib.util.module_from_spec(SPECIFICATION)
SPECIFICATION.loader.exec_module(example)


def test_simulated_stellar_populations_example_runs_pipeline(tmp_path, capsys):
    output_path = tmp_path / "simulated_stellar_populations.png"

    result = example.run_example(output_path)

    image = mpimg.imread(output_path)
    assert image.shape[0] > 100
    assert image.shape[1] > 100
    assert image[..., :3].min() < 0.8
    np.testing.assert_allclose(
        result["reference"]["massstar"][0],
        np.linspace(0.75, 1.15, 40),
    )
    output = capsys.readouterr().out
    assert "Reading from" in output
    assert f"Writing to {output_path}..." in output
