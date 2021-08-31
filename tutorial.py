import sys

import pergamon


def cnfg_tessnomi2minbholmock():
    
    # 
    pergamon.init( \
                 typeanls='tessnomi2minbholmock', \
                 )
        

def cnfg_fstrtoii():
    
    pergamon.init( \
                 typeanls='exof'
                 )
        

def cnfg_fstr():
    
    # Kepler dichotomy -- stats of multis
    pergamon.init( \
                 typeanls='tessm135nomipcanfstr'
                 )
        

def cnfg_flarspot():
    
    # Figures for XRP
    # Fig 1: stellar type vs. magnitude
    # atmospheric escape rate
    # Fig 2: completeness and recall for flare
    # Spot model:

    pergamon.init( \
                 typeanls='flarmock20scnomi'
                 )
        

# exoplanet metallicity vs mass

globals().get(sys.argv[1])()
