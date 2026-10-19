import importlib.util
from pathlib import Path

import matplotlib.image as mpimg
import numpy as np
import pytest


EXAMPLE_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "examples"
    / "compact_object_signatures.py"
)
SPECIFICATION = importlib.util.spec_from_file_location(
    "pergamon_compact_object_signatures_example", EXAMPLE_SCRIPT
)
example = importlib.util.module_from_spec(SPECIFICATION)
SPECIFICATION.loader.exec_module(example)


@pytest.mark.parametrize('typeplotback', ['white', 'dark'])
@pytest.mark.parametrize('inherit', [False, True])
def test_compact_object_signatures_example_writes_nonblank_png(tmp_path, capsys, typeplotback, inherit):
    from tdpy.plotting import plot_style

    output_path = tmp_path / "compact_object_signatures.png"

    options = {} if inherit else {'typeplotback': typeplotback}
    outer_theme = typeplotback if inherit else ('dark' if typeplotback == 'white' else 'white')
    with plot_style(outer_theme):
        original = dict(example.plt.rcParams)
        signatures = example.run_example(output_path, **options)
        assert dict(example.plt.rcParams) == original

    image = mpimg.imread(output_path)
    assert image.shape[0] > 100
    assert image.shape[1] > 100
    assert image[..., :3].min() < 0.8
    assert set(signatures) == {"beaming", "ellipsoidal", "self_lensing"}
    np.testing.assert_array_equal(image[0, 0], (0, 0, 0, 1) if typeplotback == 'dark' else (1, 1, 1, 1))
    assert f"Writing to {output_path}..." in capsys.readouterr().out