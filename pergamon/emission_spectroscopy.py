"""NASA Exoplanet Archive workflows for atmospheric target prioritization."""

from __future__ import annotations

from tdpy.verbosity import print

import csv
import io
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import matplotlib.pyplot as plt
import numpy as np
import tdpy

from .exoplanet_catalog import TAP_URL


ARCHIVE_QUERY = (
    "select pl_name,hostname,pl_eqt,st_teff,pl_rade,st_rad,sy_kmag "
    "from pscomppars where tran_flag=1 "
    "and pl_eqt is not null and st_teff is not null "
    "and pl_rade is not null and st_rad is not null and sy_kmag is not null"
)


def fetch_esm_catalog(timeout_seconds: int = 60) -> list[dict[str, float | str]]:
    """Fetch published composite system parameters and calculate each planet's ESM."""

    url = f"{TAP_URL}?{urlencode({'query': ARCHIVE_QUERY, 'format': 'csv'})}"
    print(f"Reading from {url}...")
    with urlopen(url, timeout=timeout_seconds) as response:
        rows = list(csv.DictReader(io.StringIO(response.read().decode("utf-8"))))

    equilibrium_temperature = np.array([float(row["pl_eqt"]) for row in rows])
    stellar_temperature = np.array([float(row["st_teff"]) for row in rows])
    planet_radius = np.array([float(row["pl_rade"]) for row in rows])
    stellar_radius = np.array([float(row["st_rad"]) for row in rows])
    k_magnitude = np.array([float(row["sy_kmag"]) for row in rows])
    esm = tdpy.retr_esmm(
        equilibrium_temperature,
        stellar_temperature,
        planet_radius,
        stellar_radius,
        k_magnitude,
    )

    return [
        {
            "planet_name": row["pl_name"],
            "host_name": row["hostname"],
            "equilibrium_temperature_kelvin": equilibrium_temperature[index],
            "k_magnitude": k_magnitude[index],
            "esm": esm[index],
        }
        for index, row in enumerate(rows)
        if np.isfinite(esm[index]) and esm[index] > 0
    ]


@tdpy.themed_plot
def generate_esm_catalog_figure(output_path: Path, timeout_seconds: int = 60) -> Path:
    """Plot archive inputs against the calculated emission spectroscopy metric."""

    records = fetch_esm_catalog(timeout_seconds=timeout_seconds)
    temperatures = np.array([record["equilibrium_temperature_kelvin"] for record in records])
    magnitudes = np.array([record["k_magnitude"] for record in records])
    metrics = np.array([record["esm"] for record in records])

    figure, axis = plt.subplots(figsize=(7.2, 4.5), facecolor="white")
    points = axis.scatter(
        temperatures,
        metrics,
        c=magnitudes,
        cmap="viridis_r",
        s=16,
        alpha=0.75,
        edgecolors="none",
    )
    axis.set_yscale("log")
    axis.set_xlabel("Planet equilibrium temperature [K]")
    axis.set_ylabel("Emission spectroscopy metric")
    axis.grid(False)
    colorbar = figure.colorbar(points, ax=axis, pad=0.02)
    colorbar.set_label("Host-star K magnitude [mag]")
    figure.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing to {output_path}...")
    tdpy.save_plot(figure, output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)
    return output_path


@tdpy.themed_plot
def generate_esm_distribution_figure(output_path: Path, timeout_seconds: int = 60) -> Path:
    """Plot the distribution of finite positive archive emission metrics."""

    records = fetch_esm_catalog(timeout_seconds=timeout_seconds)
    metrics = np.array([record["esm"] for record in records])

    figure, axis = plt.subplots(figsize=(7.2, 4.5), facecolor="white")
    axis.hist(metrics, bins=100, color="#007A6F")
    axis.set_xlabel("Emission spectroscopy metric")
    axis.set_ylabel("Catalog planets per bin")
    axis.set_yscale("log")
    axis.grid(False)
    figure.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing to {output_path}...")
    tdpy.save_plot(figure, output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)
    return output_path


def calc_tsmmesmm(dictpopl, strgelem='comp', boolsamp=False):
    
    if boolsamp:
        numbsamp = 1000
    else:
        numbsamp = 1

    strgradielem = 'radi' + strgelem
    strgmasselem = 'mass' + strgelem
    
    numbcomp = dictpopl[strgmasselem][0].size
    listtsmm = np.empty((numbsamp, numbcomp)) + np.nan
    listesmm = np.empty((numbsamp, numbcomp)) + np.nan
    
    for n in range(numbcomp):
        
        if not np.isfinite(dictpopl['tmpt%s' % strgelem][0][n]):
            continue
        
        if not np.isfinite(dictpopl[strgradielem][0][n]):
            continue
        
        if boolsamp:
            if not np.isfinite(dictpopl['stdvradi' + strgelem][0][n]):
                stdv = dictpopl[strgradielem][0][n]
            else:
                stdv = dictpopl['stdvradi' + strgelem][0][n]
            listradicomp = tdpy.samp_gaustrun(numbsamp, dictpopl[strgradielem][0][n], stdv, 0., np.inf)
            
            listmassplan = tdpy.samp_gaustrun(numbsamp, dictpopl[strgmasselem][0][n], dictpopl['stdvmass' + strgelem][0][n], 0., np.inf)

            if not np.isfinite(dictpopl['stdvtmpt%s' % strgelem][0][n]):
                stdv = dictpopl['tmpt%s' % strgelem][0][n]
            else:
                stdv = dictpopl['stdvtmpt%s' % strgelem][0][n]
            listtmptplan = tdpy.samp_gaustrun(numbsamp, dictpopl['tmpt%s' % strgelem][0][n], stdv, 0., np.inf)
            
            if not np.isfinite(dictpopl['stdvradistar'][0][n]):
                stdv = dictpopl['radistar'][0][n]
            else:
                stdv = dictpopl['stdvradistar'][0][n]
            listradistar = tdpy.samp_gaustrun(numbsamp, dictpopl['radistar'][0][n], stdv, 0., np.inf)
            
            listmagtsystKbnd = tdpy.icdf_gaus(np.random.rand(numbsamp), dictpopl['magtsystKbnd'][0][n], dictpopl['stdvmagtsystKbnd'][0][n])
            listmagtsystJbnd = tdpy.icdf_gaus(np.random.rand(numbsamp), dictpopl['magtsystJbnd'][0][n], dictpopl['stdvmagtsystJbnd'][0][n])
            listtmptstar = tdpy.samp_gaustrun(numbsamp, dictpopl['tmptstar'][0][n], dictpopl['stdvtmptstar'][0][n], 0., np.inf)
        
        else:
            listradicomp = dictpopl[strgradielem][0][None, n]
            listtmptplan = dictpopl['tmpt%s' % strgelem][0][None, n]
            listmassplan = dictpopl[strgmasselem][0][None, n]
            listradistar = dictpopl['radistar'][0][None, n]
            listmagtsystKbnd = dictpopl['magtsystKbnd'][0][None, n]
            listmagtsystJbnd = dictpopl['magtsystJbnd'][0][None, n]
            listtmptstar = dictpopl['tmptstar'][0][None, n]
        
        # TSM
        listtsmm[:, n] = tdpy.retr_tsmm(listradicomp, listtmptplan, listmassplan, listradistar, listmagtsystJbnd)

        # ESM
        listesmm[:, n] = tdpy.retr_esmm(listtmptplan, listtmptstar, listradicomp, listradistar, listmagtsystKbnd)
        
        #if (listesmm[:, n] < 1e-10).any():
        #    print('listradicomp')
        #    summgene(listradicomp)
        #    print('listtmptplan')
        #    summgene(listtmptplan)
        #    print('listmassplan')
        #    summgene(listmassplan)
        #    print('listradistar')
        #    summgene(listradistar)
        #    print('listmagtsystKbnd')
        #    summgene(listmagtsystKbnd)
        #    print('listmagtsystJbnd')
        #    summgene(listmagtsystJbnd)
        #    print('listtmptstar')
        #    summgene(listtmptstar)
        #    print('listesmm[:, n]')
        #    summgene(listesmm[:, n])
        #    raise Exception('')
    for name, samples in [('tsmm', listtsmm), ('esmm', listesmm)]:
        valid = np.isfinite(samples).any(axis=0)
        median = np.full(numbcomp, np.nan)
        stdv = np.full(numbcomp, np.nan)
        median[valid] = np.nanmedian(samples[:, valid], axis=0)
        stdv[valid] = np.nanstd(samples[:, valid], axis=0)
        tdpy.setp_dict(dictpopl, name, median)
        tdpy.setp_dict(dictpopl, 'stdv' + name, stdv)
