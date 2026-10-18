"""Legacy BHOL/TIC output helper.

This is not part of the active scientific API. It is retained only as a
historical helper for collaborator-specific output generation.
"""

import os

import matplotlib.pyplot as plt
from astroquery.mast import Catalogs
from tdpy.paths import RepositoryPaths


OUTDIR = str(RepositoryPaths('PERGAMON_PATH').get_visuals_path() / 'bhol')
os.makedirs(OUTDIR, exist_ok=True)


def tic_query(integer):
    tic_ID = 'TIC ' + str(integer)
    bhol = Catalogs.query_object(tic_ID, catalog='TIC')
    gaia_ID = bhol[0]['GAIA']
    radius = bhol[0]['rad']
    temperature = bhol[0]['Teff']
    tmag = bhol[0]['Tmag']
    mass = bhol[0]['mass']
    ra = bhol[0]['ra']
    dec = bhol[0]['dec']
    distance = bhol[0]['d']
    tess_info = [int(gaia_ID), radius, temperature, tmag, mass, ra, dec, distance]
    return tess_info


def main(tic_id=260647166):
    tess_info = tic_query(tic_id)
    text = (
        f'TIC:{tic_id}  $R_{{*}}$:{tess_info[1]:.2f}  Teff:{tess_info[2]:.0f}  \n'
        f'Tmag:{tess_info[3]:.1f}  $m_{{*}}$:{tess_info[4]}  Ra:{tess_info[5]:.2f}  '
        f'Dec:{tess_info[6]:.2f}  d:{tess_info[7]:.1f}'
    )

    plt.plot()
    plt.xlim(-5, 5)
    plt.ylim(-5, 5)
    plt.text(-4, 0, text, fontsize=8)
    plt.savefig(os.path.join(OUTDIR, 'TIC_string_test.pdf'))
    return text


if __name__ == '__main__':
    main()
