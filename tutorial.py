import sys

import pergamon


def cnfg_tessnomi2minbcanmocktoyy():
    
    # 
    pergamon.init( \
                 typeanly='tessnomi2minbcanmocktoyy', \
                 )
        

def cnfg_fstr():
    
    # Kepler dichotomy -- stats of multis
    pergamon.init( \
                 typeanly='pcanmockfstrnomi'
                 #typeanly='tessnomi2minbcanmocktoyy', \
                 )
        

def cnfg_flarspot():
    
    # Figures for XRP
    # Fig 1: stellar type vs. magnitude
    # atmospheric escape rate
    # Fig 2: completeness and recall for flare
    # Spot model:

    pergamon.init( \
                 typeanly='flarmock20scnomi'
                 )
        

# exoplanet metallicity vs mass

globals().get(sys.argv[1])()
