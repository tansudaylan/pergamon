import sys

import pergamon


def cnfg_tessnomi2minbholmock():
    
    # 
    pergamon.init( \
                 typeanls='tessnomi2minbholmock', \
                 )
        

def cnfg_featsupntess():
    '''
    features of supernovae in the TESS FOV
    '''
    
    pergamon.init( \
                  typeanls='featsupntess'
                 )
    

def cnfg_featexartran():
    '''
    period and duty cycle relationship of confirmed exoplanets on NASA Exoplanet Archive
    '''
    
    pergamon.init( \
                 typeanls='featexartran'
                 )
    

def cnfg_featexarmassradi():
    '''
    mass and radii relationship of confirmed exoplanets on NASA Exoplanet Archive
    '''
    
    pergamon.init( \
                 typeanls='featexarmassradi'
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
    

def cnfg_feattceefstr():
    '''
    QLP TCE feature space for
        all TCEs 
        all TOIs
        FaintStar TOIs
    '''
    
    pergamon.init( \
                 typeanls='feattceefstr'
                 )
        

def cnfg_feattoiiatmo():
    '''
    TOI TSMs binned in equatorial latitude
    '''
    
    pergamon.init( \
                 typeanls='feattoiiatmo'
                 )
        

def cnfg_featexaratmo():
    '''
    Collect relevant features for all confirmed exoplanets on the NASA Exoplanet Archive suitable for atmospheric characterization
    '''
    
    pergamon.init( \
                 typeanls='featexaratmo'
                 )
        

def cnfg_featexarweakmass():
    '''
    Collect relevant features for all confirmed exoplanets on the NASA Exoplanet Archive suitable for atmospheric characterization with weak mass measurements
    '''
    
    pergamon.init( \
                 typeanls='featexarweakmass'
                 )
        

def cnfg_featexaratmogeor():
    '''
    TSMs and ESMs for hand-picked confirmed exoplanets on the NASA Exoplanet Archive
    '''
    
    pergamon.init( \
                 typeanls='featexaratmogeor'
                 )
        

def cnfg_feattoiifstr():
    '''
    TOI feature space for
        all TOIs
        FaintStars TOIs
    '''
    
    pergamon.init( \
                 typeanls='feattoiifstr'
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
