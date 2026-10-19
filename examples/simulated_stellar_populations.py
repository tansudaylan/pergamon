#!/usr/bin/env python3
"""Compare two clearly labeled simulated stellar populations with Pergamon."""

from tdpy.verbosity import print

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
from tdpy.cli import add_plot_arguments
from tdpy.plotting import plot_background_colors, themed_plot


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


@themed_plot
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
                "pmar_radistar_massstar_referenceenhanced"
                + ("_dark" if plot_background_colors()[0] == "black" else "")
                + "_scatscat"
                f"{output_path.suffix}"
            )
        )
        print(f"Reading from {pipeline_path}...")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"Writing to {output_path}...")
        shutil.copyfile(pipeline_path, output_path)
    return result


@themed_plot
def run_bokeh_example(output_path: Path) -> Path:
    """Write the deterministic mass-radius comparison as interactive HTML."""
    from bokeh.models import ColumnDataSource, HoverTool
    from bokeh.plotting import figure
    from tdpy.plotting import save_bokeh_html

    populations = make_simulated_populations()
    plot = figure(
        title="Simulated stellar populations",
        width=760,
        height=480,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,reset,save",
    )
    for name, label, color in (
        ("reference", "Reference", "#1B6CA8"),
        ("enhanced", "Radius-enhanced", "#D18B00"),
    ):
        mass = populations[name]["massstar"][0]
        radius = populations[name]["radistar"][0]
        source = ColumnDataSource({"mass": mass, "radius": radius})
        plot.scatter(
            "mass", "radius", source=source, size=8, alpha=0.75,
            color=color, legend_label=label,
        )
    plot.add_tools(HoverTool(tooltips=[("Mass [M_sun]", "@mass{0.000}"), ("Radius [R_sun]", "@radius{0.000}")]))
    plot.xaxis.axis_label = "Stellar mass [M_sun]"
    plot.yaxis.axis_label = "Stellar radius [R_sun]"
    plot.grid.visible = False
    plot.legend.location = "top_left"
    return Path(save_bokeh_html(plot, output_path, title="Simulated stellar populations"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_plot_arguments(parser)
    parser.add_argument('--typeplotback', choices=('white', 'dark'), default='white')
    parser.add_argument("--bokeh", action="store_true", help="also write interactive HTML")
    arguments = parser.parse_args()
    output_path = Path(__file__).with_name(
        f"simulated_stellar_populations.{arguments.typefileplot}"
    )
    run_example(output_path, typeplotback=arguments.typeplotback)
    if arguments.bokeh:
        run_bokeh_example(Path(__file__).with_name("simulated_stellar_populations.html"), typeplotback=arguments.typeplotback)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
