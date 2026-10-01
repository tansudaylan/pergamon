"""Confirmed-planet tables from the NASA Exoplanet Archive."""

import argparse
from datetime import datetime, timezone
from itertools import combinations
import json
from pathlib import Path
from urllib.parse import quote, urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.constants import G
from astropy.coordinates import SkyCoord
from ephesos.models import mutual_hill_separations

from tdpy.verbosity import print
from .paths import get_data_path


TAP_URL = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
MAST_SPECTRA_URL = "https://exo.mast.stsci.edu/api/v0.1/spectra"
PLANET_COLUMNS = (
    "pl_name", "hostname", "gaia_dr3_id", "ra", "dec", "glon", "glat",
    "sy_pnum", "discoverymethod", "pl_orbper", "pl_orbsmax", "pl_orbeccen",
    "pl_bmasse", "pl_bmassprov", "pl_masse", "pl_rade", "st_mass", "st_rad",
    "st_rotp", "st_lum", "pl_insol", "pl_ntranspec", "pl_nespec", "pl_ndispec",
    "pl_orbperlim", "pl_masselim", "pl_radelim", "st_masslim", "pl_orbeccenlim",
)


def _measured_values(frame: pd.DataFrame, name: str, limit_name: str) -> pd.Series:
    """Exclude reported bounds while retaining the original archive columns."""

    values = pd.to_numeric(frame[name], errors="coerce")
    if limit_name in frame:
        values = values.where(pd.to_numeric(frame[limit_name], errors="coerce").eq(0))
    return values


def load_confirmed_exoplanets(limit: int | None = None) -> pd.DataFrame:
    """Read default published solutions from the NASA Exoplanet Archive PS table."""

    if limit is not None and (isinstance(limit, bool) or not isinstance(limit, int) or limit < 1):
        raise ValueError("limit must be a positive integer")
    prefix = f"top {limit} " if limit is not None else ""
    query = f"select {prefix}{','.join(PLANET_COLUMNS)} from ps where default_flag=1"
    url = f"{TAP_URL}?{urlencode({'query': query, 'format': 'csv'})}"
    print(f"Reading from {url}...")
    catalog = pd.read_csv(url)
    catalog["archive_source_url"] = url
    catalog["archive_retrieved_at_utc"] = datetime.now(timezone.utc).isoformat()
    return catalog


def enrich_exoplanets(planets: pd.DataFrame) -> pd.DataFrame:
    """Add measured-coordinate and true-mass orbital quantities to archive rows."""

    required = {"pl_name", "hostname", "ra", "dec", "pl_orbper", "pl_masse",
                "pl_rade", "st_mass", "pl_orbeccen"}
    missing = required - set(planets)
    if missing:
        raise ValueError(f"Missing archive columns: {sorted(missing)}")
    catalog = planets.copy()
    ra = pd.to_numeric(catalog["ra"], errors="coerce").to_numpy(dtype=float)
    dec = pd.to_numeric(catalog["dec"], errors="coerce").to_numpy(dtype=float)
    valid_position = np.isfinite(ra) & np.isfinite(dec) & (dec >= -90) & (dec <= 90)
    catalog["galactic_longitude_deg"] = np.nan
    catalog["galactic_latitude_deg"] = np.nan
    if valid_position.any():
        coordinates = SkyCoord(ra=ra[valid_position] * u.deg, dec=dec[valid_position] * u.deg)
        catalog.loc[valid_position, "galactic_longitude_deg"] = coordinates.galactic.l.deg
        catalog.loc[valid_position, "galactic_latitude_deg"] = coordinates.galactic.b.deg

    period = _measured_values(catalog, "pl_orbper", "pl_orbperlim").to_numpy(dtype=float) * u.day
    planet_mass = _measured_values(catalog, "pl_masse", "pl_masselim").to_numpy(dtype=float) * u.M_earth
    stellar_mass = _measured_values(catalog, "st_mass", "st_masslim").to_numpy(dtype=float) * u.M_sun
    radius = _measured_values(catalog, "pl_rade", "pl_radelim").to_numpy(dtype=float) * u.R_earth
    eccentricity = _measured_values(catalog, "pl_orbeccen", "pl_orbeccenlim").to_numpy(dtype=float)
    valid_orbit = (np.isfinite(period.value) & (period.value > 0)
                   & np.isfinite(planet_mass.value) & (planet_mass.value > 0)
                   & np.isfinite(stellar_mass.value) & (stellar_mass.value > 0))
    catalog["kepler_semimajor_axis_au"] = np.nan
    catalog["hill_radius_au"] = np.nan
    catalog["safronov_number"] = np.nan
    catalog["eccentricity_amd_lower_bound_kg_m2_s"] = np.nan
    if valid_orbit.any():
        orbit = (G * (stellar_mass[valid_orbit] + planet_mass[valid_orbit])
                 * period[valid_orbit] ** 2 / (4 * np.pi**2)) ** (1 / 3)
        catalog.loc[valid_orbit, "kepler_semimajor_axis_au"] = orbit.to_value(u.au)
        catalog.loc[valid_orbit, "hill_radius_au"] = (
            orbit * (planet_mass[valid_orbit] / (3 * stellar_mass[valid_orbit])) ** (1 / 3)
        ).to_value(u.au)
        valid_radius = valid_orbit & np.isfinite(radius.value) & (radius.value > 0)
        catalog.loc[valid_radius, "safronov_number"] = (
            (catalog.loc[valid_radius, "kepler_semimajor_axis_au"].to_numpy() * u.au)
            / radius[valid_radius] * planet_mass[valid_radius] / stellar_mass[valid_radius]
        ).to_value(u.dimensionless_unscaled)
        valid_eccentricity = valid_orbit & np.isfinite(eccentricity) & (eccentricity >= 0) & (eccentricity < 1)
        if valid_eccentricity.any():
            semimajor_axis = catalog.loc[valid_eccentricity, "kepler_semimajor_axis_au"].to_numpy() * u.au
            angular_momentum = planet_mass[valid_eccentricity] * np.sqrt(
                G * stellar_mass[valid_eccentricity] * semimajor_axis
            )
            catalog.loc[valid_eccentricity, "eccentricity_amd_lower_bound_kg_m2_s"] = (
                angular_momentum * (1 - np.sqrt(1 - eccentricity[valid_eccentricity] ** 2))
            ).to_value(u.kg * u.m**2 / u.s)
    catalog["has_inner_companion"] = pd.Series(pd.NA, index=catalog.index, dtype="boolean")
    catalog["has_outer_companion"] = pd.Series(pd.NA, index=catalog.index, dtype="boolean")
    if "sy_pnum" in catalog:
        for _, group in catalog.groupby("hostname"):
            periods = _measured_values(group, "pl_orbper", "pl_orbperlim")
            reported = pd.to_numeric(group["sy_pnum"], errors="coerce")
            valid = periods.notna() & (periods > 0)
            complete = (reported.notna().all() and reported.nunique() == 1
                        and reported.iloc[0] == len(group) and valid.all()
                        and periods.is_unique)
            for index in group.index[valid]:
                for direction, found in (("inner", (periods[valid] < periods[index]).any()),
                                         ("outer", (periods[valid] > periods[index]).any())):
                    if found or complete:
                        catalog.at[index, f"has_{direction}_companion"] = bool(found)
    return catalog


def add_measured_host_features(catalog: pd.DataFrame, measurements: pd.DataFrame) -> pd.DataFrame:
    """Attach published host measurements and calculate dependent planet features."""

    required = {"hostname", "source_reference"}
    if not required.issubset(measurements) or measurements["hostname"].isna().any():
        raise ValueError("Host measurements need hostname and source_reference")
    if measurements["hostname"].duplicated().any() or measurements["source_reference"].isna().any():
        raise ValueError("Host measurements need unique names and nonnull source references")
    allowed = {"host_xuv_luminosity_erg_s", "host_convective_turnover_days",
               "host_alfven_radius_stellar_radii", "host_flare_rate_per_day"}
    if not set(measurements).issubset(required | allowed):
        raise ValueError("Unsupported host measurement columns")
    result = catalog.merge(measurements.rename(columns={"source_reference": "host_source_reference"}),
                           on="hostname", how="left", validate="many_to_one")
    axis = pd.to_numeric(result["kepler_semimajor_axis_au"], errors="coerce").to_numpy(dtype=float)
    if "host_xuv_luminosity_erg_s" in result:
        luminosity = pd.to_numeric(result["host_xuv_luminosity_erg_s"], errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(luminosity) & (luminosity >= 0) & np.isfinite(axis) & (axis > 0)
        result["xuv_flux_erg_s_cm2"] = np.nan
        result.loc[valid, "xuv_flux_erg_s_cm2"] = (
            luminosity[valid] * u.erg / u.s / (4 * np.pi * (axis[valid] * u.au) ** 2)
        ).to_value(u.erg / (u.s * u.cm**2))
    if "host_convective_turnover_days" in result and "st_rotp" in result:
        turnover = pd.to_numeric(result["host_convective_turnover_days"], errors="coerce")
        rotation = pd.to_numeric(result["st_rotp"], errors="coerce")
        result["host_rossby_number"] = (rotation / turnover).where((turnover > 0) & (rotation > 0))
    if "host_alfven_radius_stellar_radii" in result and "st_rad" in result:
        alfven = pd.to_numeric(result["host_alfven_radius_stellar_radii"], errors="coerce")
        radius = pd.to_numeric(result["st_rad"], errors="coerce")
        result["host_alfven_radius_au"] = (
            (alfven * radius).where((alfven > 0) & (radius > 0)).to_numpy() * u.R_sun
        ).to_value(u.au)
    return result


def add_measured_orbit_orientation(catalog: pd.DataFrame, measurements: pd.DataFrame) -> pd.DataFrame:
    """Calculate full angular-momentum deficit from measured invariant-plane tilts."""

    required = {"pl_name", "inclination_invariable_plane_deg", "source_reference"}
    if not required.issubset(measurements) or measurements["pl_name"].duplicated().any():
        raise ValueError("Orbital orientations need unique planet names, tilts, and source references")
    if measurements[list(required)].isna().any().any():
        raise ValueError("Orbital orientation measurements must include provenance")
    result = catalog.merge(
        measurements[list(required)].rename(columns={"source_reference": "orientation_source_reference"}),
        on="pl_name", how="left", validate="one_to_one",
    )
    inclination = pd.to_numeric(result["inclination_invariable_plane_deg"], errors="coerce").to_numpy(dtype=float)
    eccentricity = _measured_values(result, "pl_orbeccen", "pl_orbeccenlim").to_numpy(dtype=float)
    semimajor = pd.to_numeric(result["kepler_semimajor_axis_au"], errors="coerce").to_numpy(dtype=float)
    mass = _measured_values(result, "pl_masse", "pl_masselim").to_numpy(dtype=float)
    stellar_mass = _measured_values(result, "st_mass", "st_masslim").to_numpy(dtype=float)
    valid = (np.isfinite(inclination) & (inclination >= 0) & (inclination <= 180)
             & np.isfinite(eccentricity) & (eccentricity >= 0) & (eccentricity < 1)
             & np.isfinite(semimajor) & (semimajor > 0)
             & np.isfinite(mass) & (mass > 0) & np.isfinite(stellar_mass) & (stellar_mass > 0))
    result["angular_momentum_deficit_kg_m2_s"] = np.nan
    if valid.any():
        angular_momentum = mass[valid] * u.M_earth * np.sqrt(
            G * stellar_mass[valid] * u.M_sun * semimajor[valid] * u.au
        )
        result.loc[valid, "angular_momentum_deficit_kg_m2_s"] = (
            angular_momentum * (1 - np.sqrt(1 - eccentricity[valid] ** 2)
                                * np.cos(np.deg2rad(inclination[valid])))
        ).to_value(u.kg * u.m**2 / u.s)
    return result


def load_mast_spectra(planet_name: str) -> pd.DataFrame:
    """List real curated spectra for one planet through the documented Exo.MAST API."""

    if not planet_name or not planet_name.strip():
        raise ValueError("planet_name must be nonempty")
    root = f"{MAST_SPECTRA_URL}/{quote(planet_name, safe='')}"
    url = f"{root}/filelist/"
    print(f"Reading from {url}...")
    with urlopen(url, timeout=30) as response:
        filenames = json.load(response)["filenames"]
    if not isinstance(filenames, list) or not all(isinstance(name, str) for name in filenames):
        raise ValueError("Exo.MAST returned an invalid spectra file list")
    return pd.DataFrame(({
        "pl_name": planet_name, "source": "Exo.MAST", "observation_id": name,
        "source_url": f"{root}/file/{quote(name, safe='')}",
    } for name in filenames), columns=("pl_name", "source", "observation_id", "source_url"))


def attach_atmosphere_observations(planets: pd.DataFrame, observations: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Match source-identified atmosphere records to planet names without losing duplicates."""

    required = {"pl_name", "source", "observation_id", "source_url"}
    if not required.issubset(observations) or observations[list(required)].isna().any().any():
        raise ValueError("Atmosphere records need planet, source, observation ID, and URL")
    if planets["pl_name"].duplicated().any() or observations.duplicated(["source", "observation_id"]).any():
        raise ValueError("Planet names and source observation IDs must be unique")
    if not observations.empty and (observations[list(required)] == "").any().any():
        raise ValueError("Atmosphere source identifiers cannot be empty")
    records = observations.merge(planets[["pl_name", "hostname"]], on="pl_name", how="left",
                                 validate="many_to_one", indicator=True)
    records["matched_archive_planet"] = records.pop("_merge").eq("both")
    counts = records.loc[records["matched_archive_planet"]].groupby("pl_name").size()
    enriched = planets.copy()
    enriched["atmosphere_observation_count"] = enriched["pl_name"].map(counts).fillna(0).astype(int)
    return enriched, records


def build_system_catalog(planets: pd.DataFrame) -> pd.DataFrame:
    """Group confirmed planet names without treating a partial download as a whole system."""

    if not {"hostname", "pl_name", "sy_pnum"}.issubset(planets):
        raise ValueError("System catalogs require hostname, pl_name, and sy_pnum")
    if planets["pl_name"].duplicated().any() or planets["hostname"].isna().any():
        raise ValueError("Planet names must be unique and host names must be known")
    groups = []
    for hostname, group in planets.groupby("hostname", sort=True):
        reported = pd.to_numeric(group["sy_pnum"], errors="coerce")
        complete = (reported.notna().all() and reported.nunique() == 1
                    and reported.iloc[0] == len(group))
        groups.append({"hostname": hostname, "planet_names": tuple(group["pl_name"]),
                       "downloaded_planet_count": len(group),
                   "archive_planet_count": reported.iloc[0] if reported.notna().all() and reported.nunique() == 1 else np.nan,
                       "system_complete": bool(complete)})
    return pd.DataFrame(groups, columns=("hostname", "planet_names", "downloaded_planet_count",
                                         "archive_planet_count", "system_complete"))


def build_planet_pair_catalog(planets: pd.DataFrame, neighbors_only: bool = False) -> pd.DataFrame:
    """List same-host pairs; mark neighbors only when the full ordered system is known."""

    required = {"hostname", "pl_name", "pl_orbper", "pl_masse", "st_mass", "sy_pnum"}
    if not required.issubset(planets):
        raise ValueError(f"Missing pair columns: {sorted(required - set(planets))}")
    systems = build_system_catalog(planets).set_index("hostname")
    pairs = []
    for hostname, group in planets.groupby("hostname", sort=True):
        ordered = group.assign(
            _period=_measured_values(group, "pl_orbper", "pl_orbperlim"),
            _mass=_measured_values(group, "pl_masse", "pl_masselim"),
            _stellar_mass=_measured_values(group, "st_mass", "st_masslim"),
        )
        ordered = ordered.sort_values(["_period", "pl_name"], na_position="last").reset_index(drop=True)
        complete = (systems.loc[hostname, "system_complete"]
                    and ordered["_period"].notna().all() and (ordered["_period"] > 0).all()
                    and ordered["_period"].is_unique)
        for inner_index, outer_index in combinations(range(len(ordered)), 2):
            inner, outer = ordered.iloc[inner_index], ordered.iloc[outer_index]
            if neighbors_only and (not complete or outer_index != inner_index + 1):
                continue
            inner_period, outer_period = inner["_period"], outer["_period"]
            valid_period = pd.notna(inner_period) and pd.notna(outer_period) and 0 < inner_period < outer_period
            masses = np.array([inner["_mass"], outer["_mass"]])
            host_masses = np.array([inner["_stellar_mass"], outer["_stellar_mass"]])
            host_mass = next((mass for mass in host_masses if np.isfinite(mass) and mass > 0), np.nan)
            separation = np.nan
            if valid_period and np.isfinite(masses).all() and (masses > 0).all() and np.isfinite(host_mass):
                separation = float(mutual_hill_separations(
                    np.array([inner_period, outer_period]), masses, host_mass
                )[0])
            pairs.append({"hostname": hostname, "inner_planet": inner["pl_name"] if valid_period else pd.NA,
                          "outer_planet": outer["pl_name"] if valid_period else pd.NA,
                          "planet_a": inner["pl_name"], "planet_b": outer["pl_name"],
                          "inner_period_days": inner_period, "outer_period_days": outer_period,
                          "period_ratio": outer_period / inner_period if valid_period else np.nan,
                          "are_neighbors": outer_index == inner_index + 1 if complete else pd.NA,
                          "mutual_hill_separation": separation,
                          "host_mass_used_solar": host_mass})
    return pd.DataFrame(pairs, columns=("hostname", "planet_a", "planet_b", "inner_planet", "outer_planet",
                                        "inner_period_days", "outer_period_days", "period_ratio",
                                        "are_neighbors", "mutual_hill_separation", "host_mass_used_solar"))


def build_exoplanet_catalogs(planets: pd.DataFrame, *, host_measurements: pd.DataFrame | None = None,
                             orbital_orientations: pd.DataFrame | None = None,
                             atmosphere_observations: pd.DataFrame | None = None) -> dict[str, pd.DataFrame]:
    """Build linked planet, host, pair, and neighbor catalogs from observed rows."""

    enriched = enrich_exoplanets(planets)
    if host_measurements is not None:
        enriched = add_measured_host_features(enriched, host_measurements)
    if orbital_orientations is not None:
        enriched = add_measured_orbit_orientation(enriched, orbital_orientations)
    tables = {}
    if atmosphere_observations is not None:
        enriched, tables["atmosphere_observations"] = attach_atmosphere_observations(
            enriched, atmosphere_observations
        )
    for name in ("xuv_flux_erg_s_cm2", "host_rossby_number", "host_alfven_radius_au",
                 "host_flare_rate_per_day", "angular_momentum_deficit_kg_m2_s",
                 "atmosphere_observation_count"):
        if name not in enriched:
            enriched[name] = np.nan
    tables.update(planets=enriched, systems=build_system_catalog(enriched),
                  pairs=build_planet_pair_catalog(enriched),
                  neighbors=build_planet_pair_catalog(enriched, neighbors_only=True))
    return tables


def main() -> None:
    """Write a dated snapshot of the live catalog and optional sourced measurements."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-directory", type=Path, default=None)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--host-measurements", type=Path)
    parser.add_argument("--orbital-orientations", type=Path)
    parser.add_argument("--atmosphere-observations", type=Path)
    arguments = parser.parse_args()

    def read_optional(path):
        if path is None:
            return None
        print(f"Reading from {path}...")
        return pd.read_csv(path)

    tables = build_exoplanet_catalogs(
        load_confirmed_exoplanets(arguments.limit),
        host_measurements=read_optional(arguments.host_measurements),
        orbital_orientations=read_optional(arguments.orbital_orientations),
        atmosphere_observations=read_optional(arguments.atmosphere_observations),
    )
    output_directory = arguments.output_directory or get_data_path()
    tag = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    paths = {name: output_directory / f"exoplanet_{name}_{tag}.csv" for name in tables}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("This exoplanet catalog snapshot already exists")
    output_directory.mkdir(parents=True, exist_ok=True)
    for name, path in paths.items():
        print(f"Writing to {path}...")
        tables[name].to_csv(path, index=False)


if __name__ == "__main__":
    main()