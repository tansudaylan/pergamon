"""Run the shared matched-filter pipeline for northern TESS sectors."""

from tdpy.verbosity import print

import runpy
from pathlib import Path


def main():
    """Run the canonical pipeline for TESS sectors 15 through 28."""

    pipeline_path = Path(__file__).with_name("mf_pipeline.py")
    print(f"Reading from {pipeline_path}...")
    runpy.run_path(str(pipeline_path), init_globals={"SECTORS": range(15, 29)})


if __name__ == "__main__":
    main()
