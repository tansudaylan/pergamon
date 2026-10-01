# Pergamon

Pergamon is a population-analysis and visualization package for astrophysical samples. Its scientific role is to take a population dictionary of sample properties, summarize the key features, visualize the one- and two-dimensional marginal distributions, and diagnose how those populations compare to one another.

## Purpose

The package is intended for population-level analysis of astrophysical systems, including comparisons across samples, feature extraction, and diagnostic visualization of multi-dimensional parameter distributions. It is especially suited to checking whether samples differ in physically meaningful ways and produces summary plots that make the underlying population structure visible.

## Population analysis

Pergamon compares astrophysical samples across measured or inferred properties, examines one- and two-dimensional distributions, identifies feature relationships, and generates summary figures that expose population differences.

## Exoplanet catalogs

Run `python -m pergamon.exoplanet_catalog` after setting `PERGAMON_PATH` to download current NASA Exoplanet Archive default planetary-system solutions. Pergamon writes UTC-tagged planet, host-system, all-pair, and neighboring-pair CSVs under `$PERGAMON_PATH/data/`. Existing snapshots are retained. Pass `--limit N` for a partial download; a partial host cannot establish the absence of an inner or outer companion. Archive row URLs and retrieval times are included with each planet.

The planet table keeps published values and their limit flags. It adds Galactic coordinates from the measured sky position, a Keplerian semimajor axis, planetary Hill radius, Safronov number, and an eccentricity-only lower bound on angular-momentum deficit. These orbital quantities require reported true mass (`pl_masse`), not the archive's best mass (`pl_bmasse`, which may be a minimum mass). Null or upper/lower-limited inputs leave derived values blank. Invariant-plane inclinations are not supplied by the archive, so full angular-momentum deficit requires measured orbital orientations.

For host activity and stellar-wind quantities, supply independently measured CSV files using `--host-measurements` (`hostname`, `source_reference`, plus any of `host_xuv_luminosity_erg_s`, `host_convective_turnover_days`, `host_alfven_radius_stellar_radii`, and `host_flare_rate_per_day`). The supplied XUV luminosity and orbital scale yield XUV flux; measured rotation and turnover periods yield Rossby number; an independently determined Alfven radius is converted to astronomical units. Use `--orbital-orientations` (`pl_name`, `inclination_invariable_plane_deg`, `source_reference`) for full angular-momentum deficit. These values remain blank unless supported by measured inputs.

The system table lists confirmed planets by host. The pair table preserves every same-host pair and calculates mutual-Hill separations with Ephesos when true masses and orbital periods are available. The neighbor table includes only adjacent planets in complete, period-ordered systems; planets with unknown order remain in the pair table without invented inner/outer labels. The reported stellar mass used for each separation is retained in the pair table.

The archive supplies published transmission, eclipse, and direct-imaging spectrum counts. `pergamon.exoplanet_catalog.load_mast_spectra(planet_name)` retrieves a real Exo.MAST curated-spectra file list with file URLs. Additional source-identified atmosphere observations can be joined with `--atmosphere-observations` using `pl_name`, `source`, `observation_id`, and `source_url`; unmatched names remain visible in the long-form observation table. No automatic ExoAtmospheres Database or MAST EAOT download is claimed because a stable, verified machine-readable export for either was not identified. Counts from such sources remain blank unless their records are supplied.

## Occurrence rates

`pergamon.estimate_occurrence_rate(detected, detection_efficiency)` estimates the fraction of surveyed targets hosting an event. Provide one zero-or-one detection flag and one independently determined detection probability per target. For example, detections `[1, 0, 0]` and efficiencies `[1, 0.5, 0.5]` yield a maximum-likelihood occurrence fraction of about 0.67. These numbers illustrate the calculation and are not observational results. `pergamon.log_likelihood_occurrence_rate` evaluates the same model for posterior inference.

The calculation assumes at most one event per target, independent targets, known detection efficiencies, and no false positives. Classification alone does not establish detection efficiency. `pergamon.partition_classified_population` groups relevant and irrelevant targets by positive and negative classifications, including true and false positives and negatives, without inferring an occurrence rate from those groups.

For compact-object companion populations, Miletos derives Doppler-beaming, ellipsoidal-variation, and self-lensing amplitudes, along with transit duration, orbital scale, and Schwarzschild radius, through Ephesos. Pergamon reads these model-derived values for population-level comparisons and occurrence-rate estimation. They are not detections or measurements of survey completeness.

## Installation

```bash
cd pergamon
python -m pip install -e .
export PERGAMON_PATH=/path/to/pergamon
```

`PERGAMON_PATH` identifies the repository root. Keep runtime inputs in `data/` and generated pipeline outputs in `visuals/`; both directories are ignored by Git.

## Example workflow

The runnable example passes two deterministic simulated stellar samples through `pergamon.init(...)` and retains the pipeline's mass-radius comparison:

```bash
python examples/simulated_stellar_populations.py --typefileplot png
```

![Pergamon simulated stellar population comparison](examples/simulated_stellar_populations.png)

The reference sample spans 0.75 to 1.15 solar masses. The radius-enhanced sample spans 0.95 to 1.45 solar masses and follows a steeper deterministic mass-radius relation. Both populations are simulated inputs rather than observed stars. Pergamon constructs the common feature matrix and produces the displayed comparison through its standard population-plotting pipeline.

## Dependencies

The package uses:

- `numpy`
- `pandas`
- `matplotlib`
- `tdpy`
- `ephesos`
- `miletos`
- `nicomedia`

## Output behavior

Pergamon produces marginal distributions, feature-pair comparisons, and population summary figures from measured or inferred sample properties.

