import importlib.util
from pathlib import Path

import matplotlib.image as mpimg


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


def test_compact_object_signatures_example_writes_nonblank_png(tmp_path, capsys):
    output_path = tmp_path / "compact_object_signatures.png"

    signatures = example.run_example(output_path)

    image = mpimg.imread(output_path)
    assert image.shape[0] > 100
    assert image.shape[1] > 100
    assert image[..., :3].min() < 0.8
    assert set(signatures) == {"beaming", "ellipsoidal", "self_lensing"}
    assert f"Writing to {output_path}..." in capsys.readouterr().out