import sys
import matplotlib
matplotlib.use('Agg')
import pergamon


def cnfg_TOICatalog():
    '''
    all features in the TOI Catalog with all subpopulations
    '''
    
    listtypelang = ['English']
    #listtypelang = ['Turkish']
    
    listnamefeatcumu = ['yearaler', 'yeardisc']
    typeanls = 'toii'
    for typelang in listtypelang:
        pergamon.init( \
                      typeanls=typeanls, \
                      listnamefeatcumu=listnamefeatcumu, \
                      typelang=typelang, \
                     )
    

def cnfg_NEA():
    '''
    exoplanet features in the NASA Exoplanet Archive (NEA)
    '''
    
    #listtypelang = ['English']
    listtypelang = ['Turkish']
    
    listnamefeatcumu = ['yearaler', 'yeardisc']
    
    dictplotgrid = dict()
    dictplotgrid['listlablannosamp'] = ['HD 108236 b', 'HD 108236 c', 'HD 108236 d', 'HD 108236 e', 'HD 108236 f']

    typeanls = 'NASA_Exoplanet_Archive'
    for typelang in listtypelang:
        pergamon.init( \
                      typeanls=typeanls, \
                      listnamefeatcumu=listnamefeatcumu, \
                      typelang=typelang, \
                      dictplotgrid=dictplotgrid, \
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
    

def cnfg_vetting():
    '''
    Features of QLP TCEs
    '''
    
    pergamon.init( \
                 'vetting'
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
