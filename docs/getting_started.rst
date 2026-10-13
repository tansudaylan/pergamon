Getting started
===============

Pergamon summarizes astrophysical populations and estimates occurrence only
when per-target detection efficiencies are supplied. It also builds confirmed
exoplanet and planet-pair catalogs from the NASA Exoplanet Archive.

Installation
------------

Install Pergamon from its repository root:

.. code-block:: console

   python -m pip install -e .

Set ``PERGAMON_PATH`` to the repository root. Downloaded catalog snapshots are
stored under ``data/`` and generated figures under ``visuals/``.

Candidate vetting
-----------------

Sardis's functionality now lives in ``pergamon.vetting``. Use
``run_candidate_vetting`` instead of ``sardis.init``, and set ``PERGAMON_PATH``
instead of ``SARDIS_PATH``. Explicit base paths, configuration directories,
backend arguments, classification summaries, and return conventions are
preserved. The module retains ``init`` and ``retr_pathsard`` aliases.

.. code-block:: console

   python -m pip install -e '.[vetting,bokeh]'
   python examples/synthetic_vetting_performance.py --typefileplot png
   python examples/synthetic_vetting_performance.py --typefileplot png --bokeh

This preserves the full Troia workflow, including its configured Miletos
analysis. The example can perform posterior inference for 100 targets.
Rendering can also use an explicitly supplied result dictionary through
``run_example(output_path, result=...)`` without running inference.

``summarize_vetting_result`` returns classifier labels, confusion matrices,
precision, recall, and a positive-truth mask. Prediction vectors must match
the truth labels. Classification alone does not measure completeness or
occurrence. The generic population ``pergamon.init`` function is unchanged.

Occurrence-rate calculation
---------------------------

Supply one binary detection flag and one independently measured efficiency per
surveyed target:

.. code-block:: python

   from pergamon import estimate_occurrence_rate

   occurrence = estimate_occurrence_rate(
       detected=[1, 0, 0],
       detection_efficiency=[1.0, 0.5, 0.5],
   )

Efficiencies describe the survey selection function and should come from
injection recovery or another explicit completeness analysis. The likelihood
assumes independent targets and at most one event per target. See :doc:`api`
for false-positive handling and compact-object signatures.

Exoplanet catalogs
------------------

Retrieve the current confirmed-planet table and write derived planet, host,
pair, and neighboring-pair tables with:

.. code-block:: console

   python -m pergamon.exoplanet_catalog --help

The catalog builder retains published values and limit flags. Derived true-mass
orbital quantities remain blank when the required measured values are absent.
Use a retained CSV snapshot with the demographics command when a reproducible
offline plot is needed.

Atmospheric target prioritization
---------------------------------

The emission spectroscopy metric (ESM) workflow retrieves confirmed transiting
planets from the archive's composite planetary-system table and calculates
metrics with TDpy. The composite table can combine parameters from different
published solutions. These functions were transferred from Abydos.

.. code-block:: python

   from pergamon import fetch_esm_catalog

   records = fetch_esm_catalog()

The returned records contain planet and host names, equilibrium temperatures,
host K magnitudes, and finite positive metrics. Plot the catalog or its ESM
distribution with:

.. code-block:: console

   python examples/plot_archive_esm.py --typefileplot png
   python examples/plot_archive_esm.py --distribution --typefileplot png

Use ``--typefileplot pdf`` for vector figures. The metric ranks relative
observability under the assumed thermal-emission model. It is not an
atmospheric detection prediction, and the histogram is not corrected for
survey completeness. Target-level photometry and radial-velocity fitting
remain in Miletos.