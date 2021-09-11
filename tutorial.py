import sys

import pergamon


def cnfg_tessnomi2minbholmock():
    
    # 
    pergamon.init( \
                 typeanls='tessnomi2minbholmock', \
                 )
        

def cnfg_micc():
    '''
    Merger-induced cora collapse supernovae
    '''
    
    pergamon.init( \
                 typeanls='micc'
                 )
    

def cnfg_multorbt():
    '''
    Kepler dichotomy -- stats of multis
    '''
    
    pergamon.init( \
                 typeanls='toii'
                 )
    

def cnfg_fstr():
    '''
    FaintStars vetting
    '''
    
    pergamon.init( \
                 typeanls='tesspcanfstr'
                 )
        

def cnfg_toiifstr():
    '''
    Contribution of FaintStars to the TOI catalog
    '''
    
    pergamon.init( \
                 typeanls='toiifstr'
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
