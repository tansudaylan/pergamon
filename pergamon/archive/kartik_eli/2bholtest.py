"""Legacy BHOL/Gaia exploratory script.

This file is retained for provenance only; it is not part of the active
scientific pipeline and should not be imported as part of the supported Pergamon API.
It is intentionally archived and should only be executed explicitly by a user
who is working with the historical prototype material.
"""

import importlib.util
import sys


if __name__ != '__main__':
    raise ImportError(
        'pergamon.archive.kartik_eli is archived legacy prototype code and is not part of the '
        'supported Pergamon API. Import it only via an explicit legacy path if you are '
        'actively reproducing historical analyses.'
    )

import os

import numpy as np
import pandas as pd
from astropy.io import ascii
from astroquery.mast import Catalogs
from tdpy.paths import RepositoryPaths


def main():
    """Run the legacy workflow only when invoked as a script."""
    base_dir = str(RepositoryPaths('PERGAMON_PATH').get_data_path() / 'bhol')
    os.makedirs(base_dir, exist_ok=True)

    path1 = os.path.join(base_dir, 'GaiaSource_6714230465835878784_6917528443525529728.csv')
    path2 = os.path.join(base_dir, 'GaiaSource_5933051914143228928_6714230117939284352.csv')
    path3 = os.path.join(base_dir, 'GaiaSource_5502601873595430784_5933051501826387072.csv')
    path4 = os.path.join(base_dir, 'GaiaSource_4475722064104327936_5502601461277677696.csv')
    path5 = os.path.join(base_dir, 'GaiaSource_3650805523966057472_4475721411269270528.csv')
    path6 = os.path.join(base_dir, 'GaiaSource_2851858288640_1584379458008952960.csv')
    path7 = os.path.join(base_dir, 'GaiaSource_2200921875920933120_3650804325670415744.csv')
    path8 = os.path.join(base_dir, 'GaiaSource_1584380076484244352_2200921635402776448.csv')

    data11through18 = [path1, path2, path3, path4, path5, path6, path7, path8]
    if not all(os.path.exists(path) for path in data11through18):
        raise FileNotFoundError(
            'Legacy BHOL data files are not present in PERGAMON_PATH/data/bhol.'
        )

    gaia_data1 = pd.read_csv(path8, index_col=None, header=0)
    gaia_data1a = np.asarray(gaia_data1['radial_velocity_error'])
    gaia_data1b = np.asarray(gaia_data1['source_id'])
    rv_index = np.where(gaia_data1a >= 15)[0]
    gaia_IDs = gaia_data1b[rv_index]

    tic_ids = np.zeros(1265, dtype=int)
    for idx in range(min(1265, len(gaia_IDs))):
        mast_table = Catalogs.query_criteria(catalog='Tic', GAIA=str(gaia_IDs[idx]), objType='STAR')
        tic_ids[idx] = int(mast_table['ID'][0])

    charar = np.chararray((1, 2), itemsize=27)
    table_data = np.zeros((1265, 2))
    charar[0] = 'TIC ID', 'Radial Velocity Uncertainty'
    for idx in range(min(1265, len(tic_ids))):
        table_data[idx] = (tic_ids[idx], gaia_data1a[rv_index][idx])

    table_data = np.concatenate((charar, table_data))
    ascii.write(table_data, os.path.join(base_dir, 'Gaia_rv_error_8.csv'), format='csv', fast_writer=False)


if __name__ == '__main__':
    main()
