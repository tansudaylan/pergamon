"""Gaia DR3 astrometry and barycentric Galactic host-star phase space."""

from datetime import datetime, timezone
import re
from urllib.parse import urlencode

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from astropy import units as u
from astropy.coordinates import SkyCoord
from tdpy.verbosity import print

GAIA_TAP_URL = 'https://gea.esac.esa.int/tap-server/tap/sync'
GAIA_COLUMNS = (
    'source_id', 'ref_epoch', 'ra', 'dec', 'parallax', 'parallax_error',
    'pmra', 'pmdec', 'radial_velocity', 'radial_velocity_error', 'ruwe',
)
_PLANET_MEASUREMENTS = {
    'pl_orbper': ('pl_orbperlim', 'Orbital period [day]'),
    'pl_rade': ('pl_radelim', 'Planet radius [Earth radius]'),
    'pl_masse': ('pl_masselim', 'Planet true mass [Earth mass]'),
}
_PHASE_SPACE_LABELS = {
    'gaia_inverse_parallax_distance_pc': 'Barycentric distance [pc]',
    **{'gaia_barycentric_galactic_' + axis + '_pc': 'Barycentric Galactic ' + axis + ' [pc]'
       for axis in ('x', 'y', 'z')},
    **{'gaia_barycentric_galactic_v' + axis + '_km_s': 'Barycentric Galactic v' + axis + ' [km/s]'
       for axis in ('x', 'y', 'z')},
}


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


def host_phase_space_rank_correlation(planets, planet_column, phase_space_column):
    """Describe a measured planet property versus one Gaia coordinate per host.

    Each host contributes its median positive, non-limit planet measurement
    among the supplied rows and one complete Gaia phase-space measurement.
    No p-value, uncertainty, selection correction, or causal claim is supplied.
    Incomplete downloads need not represent a host's full planetary system.
    """
    from .exoplanet_catalog import _measured_values

    if planet_column not in _PLANET_MEASUREMENTS or phase_space_column not in _PHASE_SPACE_LABELS:
        raise ValueError('Choose a measured period, radius, or true mass and a Gaia Cartesian component or distance')
    limit_column = _PLANET_MEASUREMENTS[planet_column][0]
    required = {'pl_name', 'gaia_source_id', 'gaia_6d_available', 'gaia_source_url',
                'gaia_retrieved_at_utc', planet_column, limit_column, phase_space_column}
    if not required.issubset(planets):
        raise ValueError('Require planet names, limit flags, joined Gaia phase space, and provenance')
    if planets['pl_name'].isna().any() or planets['pl_name'].duplicated().any():
        raise ValueError('Require unique nonnull planet names')
    if not all(isinstance(value, (bool, np.bool_)) for value in planets['gaia_6d_available'].dropna()):
        raise ValueError('Gaia availability flags must be Boolean')
    values = planets.copy()
    values['gaia_source_id'] = values['gaia_source_id'].map(_source_id).astype('string')
    values['planet_value'] = _measured_values(values, planet_column, limit_column)
    values['host_value'] = pd.to_numeric(values[phase_space_column], errors='coerce')
    available = values['gaia_6d_available'].fillna(False).astype(bool)
    if values.loc[available, 'gaia_source_id'].isna().any():
        raise ValueError('Available Gaia measurements require exact source IDs')
    for name in ('gaia_source_url', 'gaia_retrieved_at_utc'):
        provenance = values.loc[available, name]
        if provenance.isna().any() or provenance.astype(str).str.strip().eq('').any():
            raise ValueError('Available Gaia measurements require retrieval provenance')
    consistent = values.groupby('gaia_source_id')[
        ['host_value', 'gaia_6d_available', 'gaia_source_url', 'gaia_retrieved_at_utc']
    ].nunique(dropna=False)
    if (consistent > 1).any().any():
        raise ValueError('Planet rows for the same host contain conflicting Gaia measurements')
    valid = (available & np.isfinite(values['host_value']) & np.isfinite(values['planet_value'])
             & (values['planet_value'] > 0))
    hosts = values.loc[valid].groupby('gaia_source_id', sort=True).agg(
        planet_value=('planet_value', 'median'), host_value=('host_value', 'first'),
        measured_planet_count=('pl_name', 'size'),
        gaia_source_url=('gaia_source_url', 'first'),
        gaia_retrieved_at_utc=('gaia_retrieved_at_utc', 'first'),
    ).reset_index()
    coefficient = np.nan
    status = 'insufficient_hosts'
    if len(hosts) >= 3:
        status = 'constant_values'
        if hosts['planet_value'].nunique() > 1 and hosts['host_value'].nunique() > 1:
            coefficient = float(spearmanr(hosts['planet_value'], hosts['host_value']).statistic)
            status = 'descriptive_only'
    return {'hosts': hosts, 'planet_column': planet_column, 'phase_space_column': phase_space_column,
            'n_hosts': len(hosts), 'n_input_planets': len(planets), 'n_used_planets': int(valid.sum()),
            'spearman_rho': coefficient, 'status': status}


def plot_host_phase_space_rank_correlation(comparison, output_directory, typefileplot='png',
                                         typeplotback='white', write_bokeh=False, time_tag=None):
    """Plot the returned host medians without an inferred regression or trend.

    Static and optional self-contained Bokeh output use the same host rows.
    Existing timestamped outputs are never replaced. The axes retain the
    barycentric frame and the median refers only to supplied planet rows.
    """
    from pathlib import Path
    import matplotlib.pyplot as plt
    from tdpy.plotting import plot_style, save_plot

    if typefileplot not in ('png', 'pdf') or typeplotback not in ('white', 'dark'):
        raise ValueError('Use png or pdf output and white or dark backgrounds')
    tag = time_tag or datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    if not re.fullmatch(r'\d{8}_\d{6}', tag):
        raise ValueError('time_tag must use YYYYMMDD_HHMMSS')
    datetime.strptime(tag, '%Y%m%d_%H%M%S')
    hosts = comparison['hosts']
    if hosts.empty:
        raise ValueError('No eligible hosts to plot')
    planet_column = comparison['planet_column']
    host_column = comparison['phase_space_column']
    xlabel = _PHASE_SPACE_LABELS[host_column]
    planet_label = _PLANET_MEASUREMENTS[planet_column][1]
    ylabel = 'Host median ' + planet_label[0].lower() + planet_label[1:]
    coefficient = comparison['spearman_rho']
    rank_label = f'Spearman rho = {coefficient:.2f}' if np.isfinite(coefficient) else 'Spearman rho undefined'
    title = f"{len(hosts)} hosts, {comparison['n_used_planets']} measured planets"
    stem = f'gaia_{planet_column}_{host_column}_{typeplotback}_{tag}'
    directory = Path(output_directory)
    paths = {'static': directory / f'{stem}.{typefileplot}'}
    if write_bokeh:
        from bokeh.embed import file_html
        from bokeh.models import ColumnDataSource, HoverTool
        from bokeh.plotting import figure
        from bokeh.resources import INLINE

        paths['bokeh'] = directory / f'{stem}.html'
    if any(path.exists() for path in paths.values()):
        raise FileExistsError('Refusing to overwrite an existing Gaia comparison output')
    directory.mkdir(parents=True, exist_ok=True)
    background = 'black' if typeplotback == 'dark' else 'white'
    foreground = 'white' if typeplotback == 'dark' else 'black'
    color = '#ff7878' if typeplotback == 'dark' else '#a51417'
    with plot_style(typeplotback=typeplotback), plt.rc_context({
        'font.size': 10, 'axes.titlesize': 10, 'axes.labelsize': 10,
        'legend.fontsize': 10, 'xtick.labelsize': 10, 'ytick.labelsize': 10,
    }):
        static, axis = plt.subplots(figsize=(7.2, 4.8), layout='constrained')
        axis.scatter(hosts['host_value'], hosts['planet_value'], s=28, color=color,
                     edgecolor=foreground, linewidth=0.4, label=rank_label)
        axis.set(xlabel=xlabel, ylabel=ylabel, yscale='log', title=title)
        axis.grid(False)
        axis.legend(frameon=True, fancybox=True, framealpha=1.0)
        static.text(0.5, 0.01, 'Supplied planet rows only. Selection and measurement errors are not corrected.',
                    ha='center', fontsize=10)
        static.get_layout_engine().set(rect=(0, 0.06, 1, 0.94))
        print(f"Writing to {paths['static']}...")
        save_plot(static, paths['static'], dpi=300)
        plt.close(static)
    if write_bokeh:
        interactive = figure(width=840, height=520, sizing_mode='scale_width',
                             title=title + '. ' + rank_label, x_axis_label=xlabel,
                             y_axis_label=ylabel, y_axis_type='log',
                             background_fill_color=background, border_fill_color=background)
        interactive.scatter('host_value', 'planet_value', source=ColumnDataSource(hosts),
                            size=7, color=color, line_color=foreground)
        interactive.add_tools(HoverTool(tooltips=[('Gaia DR3', '@gaia_source_id'),
                                                 ('Measured planets', '@measured_planet_count'),
                                                 ('Host coordinate', '@host_value{0.000}'),
                                                 ('Planet median', '@planet_value{0.000}')]))
        interactive.grid.visible = False
        interactive.title.text_color = foreground
        interactive.title.text_font_size = '10pt'
        for axis in interactive.axis:
            axis.axis_line_color = foreground
            axis.axis_label_text_color = foreground
            axis.major_label_text_color = foreground
            axis.axis_label_text_font_size = '10pt'
            axis.major_label_text_font_size = '10pt'
            axis.major_tick_line_color = foreground
            axis.minor_tick_line_color = foreground
        print(f"Writing to {paths['bokeh']}...")
        paths['bokeh'].write_text(file_html(interactive, INLINE, 'Gaia host comparison'), encoding='utf-8')
    return paths