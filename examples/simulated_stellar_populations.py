#!/usr/bin/env python3
"""Compare two clearly labeled simulated stellar populations with Pergamon."""

import argparse
import contextlib
import io
import shutil
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("agg")

import matplotlib.pyplot as plt
import numpy as np

import pergamon


def make_simulated_populations() -> dict:
    """Return deterministic reference and radius-enhanced stellar samples."""

    mass_reference = np.linspace(0.75, 1.15, 40)  # [M_Sun]
    radius_reference = 0.95 * mass_reference + 0.08  # [R_Sun]
    mass_enhanced = np.linspace(0.95, 1.45, 40)  # [M_Sun]
    radius_enhanced = 1.08 * mass_enhanced + 0.03  # [R_Sun]
    return {
        "reference": {
            "massstar": [mass_reference, "M_Sun"],
            "radistar": [radius_reference, "R_Sun"],
        },
        "enhanced": {
            "massstar": [mass_enhanced, "M_Sun"],
            "radistar": [radius_enhanced, "R_Sun"],
        },
    }


def run_example(output_path: Path) -> dict:
    """Run Pergamon and retain its focused mass-radius comparison figure."""

    populations = make_simulated_populations()
    with tempfile.TemporaryDirectory() as directory:
        legacy_output = io.StringIO()
        with plt.rc_context({"savefig.dpi": 300}), contextlib.redirect_stdout(
            legacy_output
        ), contextlib.redirect_stderr(legacy_output):
            result = pergamon.init(
                typeanls="simulated_stellar_populations",
                dictpopl=populations,
                pathbase=directory,
                typefileplot=output_path.suffix.removeprefix("."),
                typeverb=0,
                booldiag=False,
                lablsampgene="star",
                listdictlablcolrpopl=[
                    {
                        "reference": ["Reference", "#1B6CA8"],
                        "enhanced": ["Radius-enhanced", "#D18B00"],
                    }
                ],
                listboolcompexcl=[True],
                listtitlcomp=["Simulated stellar populations"],
                listnameparaonlyincl=["massstar", "radistar"],
            )
        pipeline_path = (
            Path(directory)
            / "visuals"
            / (
                "pmar_radistar_massstar_referenceenhanced_scatscat"
                f"{output_path.suffix}"
            )
        )
        print(f"Reading from {pipeline_path}...")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"Writing to {output_path}...")
        shutil.copyfile(pipeline_path, output_path)
    return result


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--typefileplot",
        choices=("png", "pdf"),
        default="png",
    )
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    output_path = Path(__file__).with_name(
        f"simulated_stellar_populations.{arguments.typefileplot}"
    )
    run_example(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
