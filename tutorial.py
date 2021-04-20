import sys

import pergamon


def cnfg_tess():

    pergamon.init( \
                 )


globals().get(sys.argv[1])()
