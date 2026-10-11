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