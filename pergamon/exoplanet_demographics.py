"""Plot observed exoplanet demographics from NASA Exoplanet Archive data."""

import argparse
from pathlib import Path

import astropy.constants as const
import astropy.units as u
import matplotlib

matplotlib.use("agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd

from tdpy.verbosity import print

from .exoplanet_catalog import _measured_values, load_confirmed_exoplanets
from .paths import get_visuals_path


EARTH_DENSITY_G_CM3 = (const.M_earth / (4 * np.pi / 3 * const.R_earth**3)).to_value(
    u.g / u.cm**3
)  # [g cm^-3]
EARTH_ESCAPE_VELOCITY_KM_S = np.sqrt(
    2 * const.G * const.M_earth / const.R_earth
).to_value(u.km / u.s)  # [km s^-1]


def add_demographic_quantities(planets: pd.DataFrame) -> pd.DataFrame:
    """Add bulk density and escape speed only for measured masses and radii."""

    if not {"pl_masse", "pl_rade"}.issubset(planets):
        raise ValueError("Planet data must include pl_masse and pl_rade")
    result = planets.copy()
    mass = _measured_values(result, "pl_masse", "pl_masselim").to_numpy(dtype=float)
    radius = _measured_values(result, "pl_rade", "pl_radelim").to_numpy(dtype=float)
    valid = (np.isfinite(mass) & (mass > 0) & np.isfinite(radius) & (radius > 0))
    result["bulk_density_g_cm3"] = np.nan
    result["escape_velocity_km_s"] = np.nan
    result.loc[valid, "bulk_density_g_cm3"] = (
        EARTH_DENSITY_G_CM3 * mass[valid] / radius[valid] ** 3
    )
    result.loc[valid, "escape_velocity_km_s"] = (
        EARTH_ESCAPE_VELOCITY_KM_S * np.sqrt(mass[valid] / radius[valid])
    )
    return result


def _positive_measurements(planets: pd.DataFrame, column: str, limit_column: str | None = None):
    if column not in planets:
        return np.full(len(planets), np.nan)
    if limit_column is not None and limit_column in planets:
        return _measured_values(planets, column, limit_column).to_numpy(dtype=float)
    return pd.to_numeric(planets[column], errors="coerce").to_numpy(dtype=float)


def _save_figure(figure, path: Path, source_note: str) -> None:
    figure.subplots_adjust(left=0.16, right=0.88, bottom=0.21, top=0.96)
    figure.text(0.5, 0.035, source_note, ha="center", va="bottom", fontsize=10)
    print(f"Writing to {path}...")
    figure.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def plot_exoplanet_demographics(
    planets: pd.DataFrame,
    output_directory: str | Path,
    typefileplot: str = "png",
) -> dict[str, Path]:
    """Write focused mass-radius, shoreline, and period-radius figures."""

    if typefileplot not in {"png", "pdf"}:
        raise ValueError("typefileplot must be 'png' or 'pdf'")
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    data = add_demographic_quantities(planets)
    source_rows = data.get("archive_retrieved_at_utc", pd.Series(dtype=str)).dropna()
    retrieved = f"; retrieved {str(source_rows.iloc[0])[:10]} UTC" if not source_rows.empty else ""
    source_note = f"NASA Exoplanet Archive, Planetary Systems default solutions{retrieved}"
    paths = {
        name: output_directory / f"{name}.{typefileplot}"
        for name in ("exoplanet_mass_radius_density", "exoplanet_cosmic_shoreline",
                     "exoplanet_period_radius")
    }

    measured_mass = _positive_measurements(data, "pl_masse", "pl_masselim")
    measured_radius = _positive_measurements(data, "pl_rade", "pl_radelim")
    measured_density = data["bulk_density_g_cm3"].to_numpy(dtype=float)
    valid_mass_radius = (np.isfinite(measured_mass) & (measured_mass > 0)
                         & np.isfinite(measured_radius) & (measured_radius > 0)
                         & np.isfinite(measured_density) & (measured_density > 0))

    with plt.rc_context({"font.size": 10, "axes.edgecolor": "black", "axes.linewidth": 1.0,
                         "figure.facecolor": "white", "axes.facecolor": "white"}):
        figure, axis = plt.subplots(figsize=(7.2, 5.4))
        points = axis.scatter(
            measured_mass[valid_mass_radius], measured_radius[valid_mass_radius],
            c=measured_density[valid_mass_radius], norm=LogNorm(0.1, 100), cmap="viridis",
            s=20, edgecolors="black", linewidths=0.25, alpha=0.85,
        )
        axis.scatter([1], [1], marker="*", s=120, color="white", edgecolor="black",
                     linewidth=0.8, label="Earth")
        axis.set(xscale="log", yscale="log", xlim=(0.1, 1e4), ylim=(0.1, 30),
                 xlabel=r"Planet mass [$M_\oplus$]", ylabel=r"Planet radius [$R_\oplus$]")
        axis.grid(False)
        axis.legend(loc="lower right", frameon=True, fancybox=True, framealpha=1.0)
        colorbar = figure.colorbar(points, ax=axis, pad=0.02)
        colorbar.set_label(r"Bulk density [g cm$^{-3}$]")
        _save_figure(figure, paths["exoplanet_mass_radius_density"], source_note)

        instellation = _positive_measurements(data, "pl_insol")
        escape_speed = data["escape_velocity_km_s"].to_numpy(dtype=float)
        valid_shoreline = (valid_mass_radius & np.isfinite(instellation) & (instellation > 0)
                           & np.isfinite(escape_speed) & (escape_speed > 0))
        figure, axis = plt.subplots(figsize=(7.2, 5.4))
        points = axis.scatter(
            instellation[valid_shoreline], escape_speed[valid_shoreline],
            c=measured_density[valid_shoreline], norm=LogNorm(0.1, 100), cmap="viridis",
            s=20, edgecolors="black", linewidths=0.25, alpha=0.8,
            label="Measured mass and radius",
        )
        flux_guide = np.geomspace(1e-4, 3e4, 200)
        speed_guide = EARTH_ESCAPE_VELOCITY_KM_S * flux_guide**0.25
        axis.plot(flux_guide, speed_guide, color="black", linestyle="--", linewidth=1.2,
                  label=r"$v_{\rm esc}\propto F^{1/4}$ guide (Zahnle & Catling 2017)")
        axis.scatter([1], [EARTH_ESCAPE_VELOCITY_KM_S], marker="*", s=120, color="white",
                     edgecolor="black", linewidth=0.8, label="Earth")
        axis.set(xscale="log", yscale="log", xlim=(1e-4, 1e6), ylim=(0.5, 150),
                 xlabel=r"Incident stellar flux [$S_\oplus$]",
                 ylabel=r"Escape velocity [km s$^{-1}$]")
        axis.text(0.02, 0.03, "Earth-normalized scaling; not a fitted retention boundary",
                  transform=axis.transAxes, ha="left", va="bottom")
        axis.grid(False)
        axis.legend(loc="upper left", frameon=True, fancybox=True, framealpha=1.0)
        colorbar = figure.colorbar(points, ax=axis, pad=0.02)
        colorbar.set_label(r"Bulk density [g cm$^{-3}$]")
        _save_figure(figure, paths["exoplanet_cosmic_shoreline"], source_note)

        period = _positive_measurements(data, "pl_orbper", "pl_orbperlim")
        valid_period_radius = (np.isfinite(period) & (period > 0)
                               & np.isfinite(measured_radius) & (measured_radius > 0))
        period_bins = np.geomspace(0.1, 1e4, 45)
        radius_bins = np.geomspace(0.5, 30, 40)
        figure, axis = plt.subplots(figsize=(7.2, 5.4))
        _, _, _, image = axis.hist2d(
            period[valid_period_radius], measured_radius[valid_period_radius],
            bins=(period_bins, radius_bins), norm=LogNorm(vmin=1), cmap="magma",
        )
        axis.set(xscale="log", yscale="log", xlim=(period_bins[0], period_bins[-1]),
                 ylim=(radius_bins[0], radius_bins[-1]),
                 xlabel=r"Orbital period [day]", ylabel=r"Planet radius [$R_\oplus$]")
        axis.grid(False)
        colorbar = figure.colorbar(image, ax=axis, pad=0.02)
        colorbar.set_label("Detected planets per bin")
        _save_figure(figure, paths["exoplanet_period_radius"], source_note)

    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path,
                        help="saved NASA Exoplanet Archive planet CSV; defaults to a live query")
    parser.add_argument("--output-directory", type=Path, default=get_visuals_path() / "exoplanet_demographics")
    parser.add_argument("--typefileplot", choices=("png", "pdf"), default="png")
    arguments = parser.parse_args()
    if arguments.catalog is None:
        planets = load_confirmed_exoplanets()
    else:
        print(f"Reading from {arguments.catalog}...")
        planets = pd.read_csv(arguments.catalog)
    plot_exoplanet_demographics(planets, arguments.output_directory, arguments.typefileplot)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())