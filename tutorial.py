import sys

import pergamon


def cnfg_featpsys_s2nr_lsstwfdsfull_lsstwfds():

    '''
    Explore features of planetary systems (underlying physical features as well as those that can be measured by the 10-year LSST WFD survey) based on signal-to-noise calculation
    '''
    
    pergamon.init( \
                 typeanls='featpsys_s2nr_lsstwfds_lsstwfdsfullm180', \
                 )
        

def cnfg_featpsys_s2nr_tessffim_tessexm2():
    '''
    Explore features of planetary systems (underlying physical features as well as those that can be measured by 200-sec cadence TESS FFI data in the second extended mission) 
    based on signal-to-noise calculation
    '''
    
    pergamon.init( \
                 typeanls='featpsys_s2nr_tessffim_tessexm2m135', \
                 )
        

def cnfg_featpsys_s2nr_tess2min_tessnomi2min():
    '''
    Explore features of planetary systems (underlying physical features as well as those that can be measured by 2-min cadence TESS data in the nominal mission)
    based on signal-to-noise calculation
    '''
    
    pergamon.init( \
                 typeanls='featpsys_s2nr_tess2min_tessnomi2min', \
                 )
        

def cnfg_featcosc_s2nr_lsst_wfds_lsstdeepfull_lsstdeep():

    '''
    Explore features of COSCs (underlying physical features as well as those that can be measured by the 10-year LSST WFD survey) based on signal-to-noise calculation
    '''
    
    pergamon.init( \
                 typeanls='featcosc_s2nr_lsstwfdsfull_lsstwfds', \
                 )
        

def cnfg_featcosc_s2nr_tessffim_tessexm2():
    '''
    Explore features of COSCs (underlying physical features as well as those that can be measured by 2-min cadence TESS data in the nominal mission) based on signal-to-noise calculation
    '''
    
    pergamon.init( \
                 typeanls='featcosc_s2nr_tessffim_tessexm2m135', \
                 )
        

def cnfg_featcosc_s2nr_tess2min_tessnomi2min():
    '''
    Explore features of COSCs (underlying physical features as well as those that can be measured by 2-min cadence TESS data in the nominal mission) based on signal-to-noise calculation
    '''
    
    pergamon.init( \
                 typeanls='featcosc_s2nr_tess2min_tessnomi2min', \
                 )
        

def cnfg_featsupntess():
    '''
    features of supernovae in the TESS FOV
    '''
    
    pergamon.init( \
                  typeanls='featsupntess'
                 )
    

def cnfg_feattoii():
    '''
    all features in the TOI Catalog with all subpopulations
    '''
    
    pergamon.init( \
                 typeanls='feattoii'
                 )
    

def cnfg_featexar():
    '''
    all features in the NASA Exoplanet Archive with all subpopulations
    '''
    
    pergamon.init( \
                 typeanls='featexar'
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
        

def cnfg_featobsvjwstexop():
    '''
    features of approved JWST exoplanet programs and their modes
    '''
    
    pergamon.init( \
                 typeanls='featobsvjwstexop'
                 )
        

def cnfg_featmult():
    '''
    Collect number of transiting planets per system for all hosts of TOIs and TESS-confirmed exoplanets on NASA Exoplanet Archive
    '''
    
    #pergamon.init( \
    #             typeanls='feattoiimult'
    #             )

    pergamon.init( \
                 typeanls='featexarmult'
                 )
        
    pergamon.init( \
                 typeanls='feathosttoiimult'
                 )

    pergamon.init( \
                 typeanls='feathostexarmult'
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
        FaintStars TOIs (PM)
        FaintStars TOIs (EM1)
        FaintStars TOIs (EM1Ecl)
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
