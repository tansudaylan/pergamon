"""Gaia DR3 astrometry and barycentric Galactic host-star phase space."""

from datetime import datetime, timezone
import re
from urllib.parse import urlencode

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from tdpy.verbosity import print

GAIA_TAP_URL = 'https://gea.esac.esa.int/tap-server/tap/sync'
GAIA_COLUMNS = (
    'source_id', 'ref_epoch', 'ra', 'dec', 'parallax', 'parallax_error',
    'pmra', 'pmdec', 'radial_velocity', 'radial_velocity_error', 'ruwe',
)


def _source_id(value):
    if value is None or value is pd.NA or isinstance(value, (float, np.floating)) and np.isnan(value):
        return pd.NA
    if isinstance(value, (bool, np.bool_, float, np.floating)):
        raise ValueError('Gaia source IDs must be exact integers or decimal strings, never floats')
    text = str(value).strip()
    if not text:
        return pd.NA
    match = re.fullmatch(r'(?:Gaia\s+DR3\s+)?([0-9]+)', text, re.IGNORECASE)
    if not match or not 0 < int(match.group(1)) < 2**63:
        raise ValueError('Require a positive Gaia DR3 source ID')
    return str(int(match.group(1)))


def load_gaia_dr3_astrometry(source_ids, batch_size=50):
    """Retrieve exact requested IDs from ESA's Gaia DR3 TAP table.

    IDs are deduplicated before requests and parsed as strings in responses.
    Missing IDs remain absent from this table for a subsequent left join.
    Each returned row retains its query URL and UTC retrieval time.
    """
    if isinstance(batch_size, bool) or not isinstance(batch_size, (int, np.integer)) or not 1 <= batch_size <= 100:
        raise ValueError('batch_size must be an integer from 1 through 100')
    identifiers = sorted({normalized for value in source_ids if (normalized := _source_id(value)) is not pd.NA})
    frames = []
    for start in range(0, len(identifiers), batch_size):
        batch = identifiers[start:start + batch_size]
        query = f"select {','.join(GAIA_COLUMNS)} from gaiadr3.gaia_source where source_id in ({','.join(batch)})"
        url = f"{GAIA_TAP_URL}?{urlencode({'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv', 'QUERY': query})}"
        print(f'Reading from {url}...')
        frame = pd.read_csv(url, dtype={'source_id': 'string'})
        if not set(GAIA_COLUMNS).issubset(frame):
            raise ValueError('Gaia TAP response is missing requested astrometry columns')
        normalized = frame['source_id'].map(_source_id)
        if normalized.isna().any() or normalized.duplicated().any() or not normalized.isin(batch).all():
            raise ValueError('Gaia TAP response contains missing, duplicate, or unrequested source IDs')
        frame['source_id'] = normalized.astype('string')
        frame['gaia_source_url'] = url
        frame['gaia_retrieved_at_utc'] = datetime.now(timezone.utc).isoformat()
        frames.append(frame)
    columns = [*GAIA_COLUMNS, 'gaia_source_url', 'gaia_retrieved_at_utc']
    return pd.concat(frames, ignore_index=True)[columns] if frames else pd.DataFrame(columns=columns)


def attach_gaia_phase_space(planets, astrometry, max_fractional_parallax_error=0.2):
    """Attach DR3 measurements and complete barycentric Galactic phase space.

    Every planet row and original index is retained. Positive inverse-parallax
    distances are used only below the explicit fractional-error threshold.
    This is not a posterior distance estimate and applies no parallax
    zero-point correction. Missing velocities are never filled with zero.
    Coordinates remain at each Gaia reference epoch. Galactic axes have the
    Solar System barycentre as origin, not the Galactic centre. No solar
    motion subtraction or Galactic orbit integration is performed.
    """
    if 'gaia_dr3_id' not in planets:
        raise ValueError('Planet table requires gaia_dr3_id')
    required = {*GAIA_COLUMNS, 'gaia_source_url', 'gaia_retrieved_at_utc'}
    if not required.issubset(astrometry):
        raise ValueError('Astrometry table requires Gaia columns and retrieval provenance')
    if (isinstance(max_fractional_parallax_error, (bool, np.bool_))
            or np.ndim(max_fractional_parallax_error) != 0
            or not np.isfinite(max_fractional_parallax_error)
            or not 0 < max_fractional_parallax_error <= 0.2):
        raise ValueError('Fractional parallax-error threshold must be positive and at most 0.2')
    measurements = astrometry.copy()
    source_ids = measurements['source_id'].map(_source_id).astype('string')
    if source_ids.isna().any() or source_ids.duplicated().any():
        raise ValueError('Astrometry source IDs must be unique and nonnull')
    for name in ('gaia_source_url', 'gaia_retrieved_at_utc'):
        if measurements[name].isna().any() or measurements[name].astype(str).str.strip().eq('').any():
            raise ValueError('Astrometry rows need nonempty retrieval provenance')
    measurements.index = source_ids
    identifiers = planets['gaia_dr3_id'].map(_source_id).astype('string')
    added_columns = {'gaia_source_id', 'gaia_source_url', 'gaia_retrieved_at_utc',
                     *('gaia_' + name for name in GAIA_COLUMNS[1:]),
                     'gaia_6d_status', 'gaia_6d_available', 'gaia_fractional_parallax_error',
                     'gaia_inverse_parallax_distance_pc', 'gaia_phase_space_frame',
                     'gaia_distance_method', 'gaia_max_fractional_parallax_error',
                     *('gaia_barycentric_galactic_' + axis + '_pc' for axis in ('x', 'y', 'z')),
                     *('gaia_barycentric_galactic_v' + axis + '_km_s' for axis in ('x', 'y', 'z'))}
    if added_columns.intersection(planets.columns):
        raise ValueError('Planet table already contains attached Gaia columns')
    result = planets.copy()
    result['gaia_source_id'] = identifiers
    for name in GAIA_COLUMNS[1:]:
        result['gaia_' + name] = identifiers.map(measurements[name]).to_numpy()
    for name in ('gaia_source_url', 'gaia_retrieved_at_utc'):
        result[name] = identifiers.map(measurements[name]).to_numpy()
    fields = {name: pd.to_numeric(result['gaia_' + name], errors='coerce').to_numpy(dtype=float, na_value=np.nan)
              for name in ('ref_epoch', 'ra', 'dec', 'parallax', 'parallax_error', 'pmra', 'pmdec', 'radial_velocity')}
    matched = identifiers.isin(measurements.index).to_numpy()
    complete = matched & np.logical_and.reduce([np.isfinite(values) for values in fields.values()])
    valid_position = ((fields['ra'] >= 0) & (fields['ra'] < 360)
                      & (fields['dec'] >= -90) & (fields['dec'] <= 90))
    positive_parallax = (fields['parallax'] > 0) & (fields['parallax_error'] > 0)
    fractional_error = np.full(len(result), np.nan)
    np.divide(fields['parallax_error'], fields['parallax'], out=fractional_error,
              where=positive_parallax)
    precise = positive_parallax & (fractional_error <= max_fractional_parallax_error)
    available = complete & valid_position & precise
    status = np.full(len(result), 'unmatched', dtype=object)
    status[matched] = 'incomplete_astrometry'
    status[complete & ~valid_position] = 'invalid_position'
    status[complete & valid_position & ~positive_parallax] = 'invalid_parallax'
    status[complete & valid_position & positive_parallax & ~precise] = 'parallax_precision_cut'
    status[available] = 'complete'
    result['gaia_6d_status'] = status
    result['gaia_6d_available'] = available
    result['gaia_fractional_parallax_error'] = fractional_error
    result['gaia_inverse_parallax_distance_pc'] = np.nan
    result['gaia_phase_space_frame'] = 'Galactic axes, Solar System barycentre origin'
    result['gaia_distance_method'] = 'inverse parallax, no zero-point correction'
    result['gaia_max_fractional_parallax_error'] = float(max_fractional_parallax_error)
    positions = ['gaia_barycentric_galactic_' + axis + '_pc' for axis in ('x', 'y', 'z')]
    velocities = ['gaia_barycentric_galactic_v' + axis + '_km_s' for axis in ('x', 'y', 'z')]
    for name in positions + velocities:
        result[name] = np.nan
    if available.any():
        distance_pc = 1000 / fields['parallax'][available]  # [pc]
        coordinates = SkyCoord(
            ra=fields['ra'][available] * u.deg, dec=fields['dec'][available] * u.deg,
            distance=distance_pc * u.pc, pm_ra_cosdec=fields['pmra'][available] * u.mas / u.yr,
            pm_dec=fields['pmdec'][available] * u.mas / u.yr,
            radial_velocity=fields['radial_velocity'][available] * u.km / u.s, frame='icrs',
        ).galactic
        result.loc[available, 'gaia_inverse_parallax_distance_pc'] = distance_pc
        for name, values in zip(positions, coordinates.cartesian.xyz.to_value(u.pc)):
            result.loc[available, name] = values
        for name, values in zip(velocities, coordinates.velocity.d_xyz.to_value(u.km / u.s)):
            result.loc[available, name] = values
    return result