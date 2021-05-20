import sys

import pergamon


def cnfg_tess():
    
    # Figures for XRP
    # Fig 1: stellar type vs. magnitude
    # atmospheric escape rate
    # Fig 2: completeness and recall for flare
    # Spot model:

    pergamon.init( \
                 listtypepopl=['2minnomi'], \
                 )
        

def cnfg_spot():
    
    # Figures for XRP
    # Fig 1: stellar type vs. magnitude
    # atmospheric escape rate
    # Fig 2: completeness and recall for flare
    # Spot model:

    pergamon.init( \
                 listtypepopl=['2minnomi', '20sc'], \
                 )
        

globals().get(sys.argv[1])()
