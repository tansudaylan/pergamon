API reference
=============

Candidate vetting
-----------------

.. automodule:: pergamon.vetting
   :members: run_candidate_vetting, summarize_vetting_result, vetting_paths, retr_pathsard, init

``pergamon.init`` remains the population-analysis entry point; use
``pergamon.run_candidate_vetting`` or ``pergamon.vetting.init`` for vetting.

.. automodule:: pergamon.vetting_workflow
   :members: synthetic_target_truth, validated_classification, retr_pathvetting,
             retr_dictderi_effe, mile_work, init

Compact-object signatures
-------------------------

.. automodule:: pergamon.signatures
   :members: compute_photometric_signatures

Occurrence inference
--------------------

.. automodule:: pergamon.occurrence
   :members: log_likelihood_occurrence_rate, estimate_occurrence_rate, estimate_compact_object_occurrence_rate

Exoplanet catalog products
--------------------------

.. automodule:: pergamon.exoplanet_catalog
   :members: load_confirmed_exoplanets, enrich_exoplanets, add_measured_host_features, add_measured_orbit_orientation, load_mast_spectra, attach_atmosphere_observations, build_system_catalog, build_planet_pair_catalog, build_exoplanet_catalogs

Emission spectroscopy catalog
-----------------------------

.. automodule:: pergamon.emission_spectroscopy
   :members: fetch_esm_catalog, generate_esm_catalog_figure, generate_esm_distribution_figure, calc_tsmmesmm

Dictionary catalog adapters
---------------------------

.. automodule:: pergamon.catalogs
   :members: retr_toiitici, retr_dictpopltic8, retr_dictpoplrvel, retr_dicthostplan, retr_dicttoii, retr_toiifstr, retr_dictexar

Inherited population synthesis
------------------------------

.. automodule:: pergamon.population_synthesis
   :members: retr_dictpoplstarcomp

Population plots
----------------

.. automodule:: pergamon.exoplanet_demographics
   :members: add_demographic_quantities, plot_exoplanet_demographics

Paths
-----

.. automodule:: pergamon.paths
   :members: get_repository_path, get_data_path, get_visuals_path, retr_pathgaia