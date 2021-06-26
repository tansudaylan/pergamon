import sys

import pergamon


def cnfg_fstr():
    
    pergamon.init( \
                 typepopl='fstr', \
                 listtypepoplsubb=['fstrpcan'], \
                 )
        

def cnfg_flarspot():
    
    # Figures for XRP
    # Fig 1: stellar type vs. magnitude
    # atmospheric escape rate
    # Fig 2: completeness and recall for flare
    # Spot model:

    pergamon.init( \
                 #listtypepopl=['ffimm135nomi', 'fstr'], \
                 listtypepopl=['2minnomi', '20sc'], \
                 )
        

globals().get(sys.argv[1])()
