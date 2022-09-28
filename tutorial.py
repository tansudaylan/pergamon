import sys

import pergamon


def cnfg_toii(strgextn):
    '''
    all features in the TOI Catalog with all subpopulations
    '''
    
    typeanls = 'toii' + strgextn
    pergamon.init( \
                 typeanls=typeanls, \
                 )
    

def cnfg_exar(strgextn):
    '''
    all features in the NASA Exoplanet Archive with all subpopulations
    '''
    
    typeanls = 'exar' + strgextn
    pergamon.init( \
                 typeanls=typeanls, \
                 typelang='Turkish', \
                 )
    

def cnfg_psys():
    '''
    Explore features of planetary systems
    '''
   
    pergamon.init( \
                 typeanls='psys'
                 )


def cnfg_cosc():
    '''
    Explore features of COSCs
    '''
    
    pergamon.init( \
                 typeanls='cosc', \
                 )
        

def cnfg_ISOB():
    '''
    features of interstellar objects (ISOs)
    '''
    
    pergamon.init( \
                  typeanls='isob'
                 )
    

def cnfg_supntess():
    '''
    features of supernovae in the TESS FOV
    '''
    
    pergamon.init( \
                  typeanls='supntess'
                 )
    

def cnfg_qtce():
    '''
    Features of QLP TCEs
    '''
    
    pergamon.init( \
                 'qtce'
                 )
        

def cnfg_micc():
    '''
    Merger-induced cora collapse supernovae
    '''
    
    pergamon.init( \
                 typeanls='micc'
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
        

globals().get(sys.argv[1])(*sys.argv[2:])
