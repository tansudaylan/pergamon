#!/usr/bin/env python3
"""Evaluate Pergamon candidate vetting on its seeded synthetic population."""

from tdpy.verbosity import print

import contextlib
import io
import argparse
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("agg")

import matplotlib.pyplot as plt
import numpy as np

from pergamon import vetting
from tdpy.cli import add_plot_arguments


def run_bokeh_example(result: dict, summary: dict, output_path: Path) -> Path:
    """Write interactive truth-count, confusion-matrix, and score panels."""
    from bokeh.layouts import row
    from bokeh.models import ColorBar, ColumnDataSource, FactorRange, HoverTool, LinearColorMapper
    from bokeh.palettes import Blues256
    from bokeh.plotting import figure
    from tdpy.plotting import save_bokeh_html

    truth_types = result["listnameclastruetype"]
    truth_indices = np.asarray(result["indxclastruetypetarg"], dtype=int)
    truth_counts = np.bincount(truth_indices, minlength=len(truth_types))
    count_plot = figure(
        title="Input population",
        x_range=[name.replace("System", " System") for name in truth_types],
        width=370,
        height=360,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,reset,save",
    )
    count_source = ColumnDataSource({
        "type": [name.replace("System", " System") for name in truth_types],
        "count": truth_counts,
    })
    count_plot.vbar(x="type", top="count", width=0.7, source=count_source, color="#2E7D32")
    count_plot.add_tools(HoverTool(tooltips=[("Truth class", "@type"), ("Targets", "@count")]))
    count_plot.yaxis.axis_label = "Simulated targets"

    confusion = np.asarray(summary["confusion_matrices"][0])
    decision_labels = ["Rejected", "Selected"]
    truth_labels = ["Nonplanet", "Planet"]
    cells = [
        {"decision": decision_labels[column], "truth": truth_labels[row], "count": int(confusion[row, column])}
        for row in range(confusion.shape[0])
        for column in range(confusion.shape[1])
    ]
    confusion_source = ColumnDataSource({
        "decision": [cell["decision"] for cell in cells],
        "truth": [cell["truth"] for cell in cells],
        "count": [cell["count"] for cell in cells],
    })
    maximum = max(1, int(confusion.max()))
    mapper = LinearColorMapper(palette=Blues256, low=0, high=maximum)
    confusion_plot = figure(
        title="Low-threshold confusion matrix",
        x_range=decision_labels,
        y_range=list(reversed(truth_labels)),
        width=400,
        height=360,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,reset,save",
    )
    confusion_plot.rect(
        x="decision", y="truth", width=1, height=1, source=confusion_source,
        fill_color={"field": "count", "transform": mapper}, line_color="white",
    )
    confusion_plot.text(x="decision", y="truth", text="count", source=confusion_source,
                        text_align="center", text_baseline="middle", text_color="black")
    confusion_plot.add_layout(ColorBar(color_mapper=mapper, title="Targets"), "right")
    confusion_plot.xaxis.axis_label = "Vetting decision"
    confusion_plot.yaxis.axis_label = "Simulated truth"
    confusion_plot.add_tools(HoverTool(tooltips=[("Truth", "@truth"), ("Decision", "@decision"), ("Targets", "@count")]))

    classifier_labels = summary["classifier_labels"]
    factors = [(label, metric) for label in classifier_labels for metric in ("Precision", "Recall")]
    metric_values = [
        value
        for precision, recall in zip(summary["precision"], summary["recall"])
        for value in (float(precision), float(recall))
    ]
    colors = ["#1B6CA8" if metric == "Precision" else "#D18B00" for _, metric in factors]
    metric_source = ColumnDataSource({"factor": factors, "value": metric_values, "color": colors})
    metric_plot = figure(
        title="Vetting performance",
        x_range=FactorRange(*factors),
        y_range=(0.0, 1.05),
        width=450,
        height=360,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,reset,save",
    )
    metric_plot.vbar(x="factor", top="value", width=0.8, source=metric_source, color="color")
    metric_plot.yaxis.axis_label = "Classification metric"
    metric_plot.x_range.group_padding = 0.15
    metric_plot.add_tools(HoverTool(tooltips=[("Classifier", "@factor"), ("Score", "@value{0.00}")]))
    for plot in (count_plot, confusion_plot, metric_plot):
        plot.grid.visible = False
    return Path(save_bokeh_html(
        row(count_plot, confusion_plot, metric_plot, sizing_mode="stretch_width"),
        output_path,
        title="Pergamon synthetic vetting performance",
    ))


def run_synthetic_pipeline(pathbase: Path) -> dict:
    """Run Pergamon vetting on its deterministic synthetic target population."""

    configuration = {
        "pathbase": str(pathbase),
        "typesyst": "PlanetarySystem",
        "listlablinst": [["TESS"], []],
        "liststrgtypedata": [["simutargsynt"], []],
        "boolplot": False,
        "typeverb": 0,
    }
    legacy_output = io.StringIO()
    with contextlib.redirect_stdout(legacy_output), contextlib.redirect_stderr(
        legacy_output
    ):
        return vetting.run_candidate_vetting(
            dicttroiinpt=configuration,
            pathbase=str(pathbase),
            strgcnfg="SyntheticVetting",
        )


def run_example(
    output_path: Path,
    result: dict | None = None,
    bokeh_output_path: Path | None = None,
) -> dict[str, object]:
    """Plot truth composition and vetting performance from a Pergamon run."""

    if result is None:
        with tempfile.TemporaryDirectory() as directory:
            result = run_synthetic_pipeline(Path(directory))
    summary = vetting.summarize_vetting_result(result)

    truth_type_index = np.asarray(result["indxclastruetypetarg"], dtype=int)
    truth_counts = np.bincount(
        truth_type_index,
        minlength=len(result["listnameclastruetype"]),
    )
    classifier_labels = [
        label.replace(" ", "\n", 1) for label in summary["classifier_labels"]
    ]

    figure, axes = plt.subplots(
        1,
        3,
        figsize=(12.5, 4.2),
        facecolor="white",
        constrained_layout=True,
    )
    truth_positions = np.arange(len(truth_counts))
    axes[0].bar(truth_positions, truth_counts, color=["#2E7D32", "#6A6A6A"])
    axes[0].set_xticks(
        truth_positions,
        [name.replace("System", "\nSystem") for name in result["listnameclastruetype"]],
    )
    axes[0].set_ylabel("Simulated targets")
    axes[0].set_title("Input population")

    confusion_matrix = summary["confusion_matrices"][0]
    image = axes[1].imshow(confusion_matrix, cmap="Blues", vmin=0)
    figure.colorbar(image, ax=axes[1], label="Targets")
    axes[1].set_xticks([0, 1], ["Rejected", "Selected"])
    axes[1].set_yticks([0, 1], ["Nonplanet", "Planet"])
    axes[1].set_xlabel("Vetting decision")
    axes[1].set_ylabel("Simulated truth")
    axes[1].set_title("Low-threshold confusion matrix")
    for row_index in range(2):
        for column_index in range(2):
            axes[1].text(
                column_index,
                row_index,
                str(confusion_matrix[row_index, column_index]),
                ha="center",
                va="center",
                color="black",
                fontweight="bold",
            )

    classifier_positions = np.arange(len(classifier_labels))
    bar_width = 0.36
    axes[2].bar(
        classifier_positions - bar_width / 2.0,
        summary["precision"],
        width=bar_width,
        color="#1B6CA8",
        label="Precision",
    )
    axes[2].bar(
        classifier_positions + bar_width / 2.0,
        summary["recall"],
        width=bar_width,
        color="#D18B00",
        label="Recall",
    )
    axes[2].set_xticks(classifier_positions, classifier_labels)
    axes[2].set_ylim(0.0, 1.05)
    axes[2].set_ylabel("Classification metric")
    axes[2].set_title("Threshold tradeoff")
    legend = axes[2].legend(frameon=True, fancybox=True, framealpha=1.0)
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_edgecolor("black")

    for axis in axes:
        axis.grid(False)
    figure.suptitle(
        "Seeded simulation of TESS planetary-system vetting",
        fontweight="bold",
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing to {output_path}...")
    figure.savefig(
        output_path,
        dpi=300 if output_path.suffix == ".png" else None,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(figure)
    if bokeh_output_path is not None:
        run_bokeh_example(result, summary, bokeh_output_path)
    return summary


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plot Pergamon vetting performance on its seeded synthetic population."
    )
    add_plot_arguments(parser)
    parser.add_argument("--bokeh", action="store_true", help="also write interactive HTML")
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    output_path = Path(__file__).with_name(
        f"synthetic_vetting_performance.{arguments.typefileplot}"
    )
    html_path = Path(__file__).with_name("synthetic_vetting_performance.html") if arguments.bokeh else None
    run_example(output_path, bokeh_output_path=html_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())