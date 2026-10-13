#!/usr/bin/env python3
"""Plot emission spectroscopy metrics from current public archive data."""

import argparse
from tdpy.cli import add_plot_arguments
from pathlib import Path

from pergamon import generate_esm_catalog_figure, generate_esm_distribution_figure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_plot_arguments(parser)
    parser.add_argument("--distribution", action="store_true")
    arguments = parser.parse_args()
    basename = "archive_esm_distribution" if arguments.distribution else "archive_esm"
    output_path = Path(__file__).with_name(f"{basename}.{arguments.typefileplot}")
    generate = generate_esm_distribution_figure if arguments.distribution else generate_esm_catalog_figure
    generate(output_path)


if __name__ == "__main__":
    main()