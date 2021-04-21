import sys

import pergamon


def cnfg_tess():
    
    # Fig 1
    pergamon.init( \
                 typepopl='tic8'
                 )


globals().get(sys.argv[1])()
