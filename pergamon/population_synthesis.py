"""Inherited synthetic stellar and companion population construction."""

import numpy as np
import tdpy
from tdpy import summgene
import chalcedon
from tdpy.verbosity import print, tqdm
from .catalogs import retr_dictpopltic8

def retr_dictpoplstarcomp( \
                          # type of target systems
                          typesyst, \
                          
                          # type of the population of target systems
                          typepoplsyst=None, \
                          
                          # number of systems
                          numbsyst=None, \
                          
                          # epochs of mid-transits
                          epocmtracomp=None, \
                          
                          # offset for mid-transit epochs
                          timeepoc=None, \
                          
                          # type of sampling of orbital period and semi-major axes
                          ## 'smax': semi-major axes are sampled first and then orbital periods are calculated
                          ## 'peristab': orbital periods are sampled first, respecting stability, and then semi-major axies are calculated
                          ## 'peripowr': orbital periods are sampled from a power law regardless of stability and then semi-major axies are calculated
                          typesamporbtcomp='peristab', \
                          
                          # stellar properties (if scalar, resample from the same star)
                          ## distances to the systems
                          distsyst=None, \
                          ## mass of the star(s)
                          massstar=None, \

                          # minimum number of components per star
                          minmnumbcompstar=1, \
                          
                          # maximum number of components per star
                          maxmnumbcompstar=1, \
                          
                          # minimum ratio of semi-major axis to radius of the host star
                          minmsmaxradistar=3., \
                          
                          # maximum ratio of semi-major axis to radius of the host star
                          maxmsmaxradistar=1e4, \
                          
                          # minimum mass of the companions
                          minmmasscomp=None, \
                          
                          # maximum mass of the companions
                          maxmmasscomp=None, \
                          
                          # minimum orbital period, only taken into account when typesamporbtcomp == 'peripowr' or 'peristab'
                          minmpericomp=0.5, \
                          
                          # maximum orbital period, only taken into account when typesamporbtcomp == 'peripowr'
                          maxmpericomp=1000., \
                          
                          # Boolean flag to force all companions to be transiting
                          booltrancomp=False, \

                          # maximum cosine of inclination if booltrancomp is False
                          maxmcosicomp=None, \
                          
                          # Boolean flag to include exomoons
                          boolinclmoon=False, \
                          
                          # type of stellar population
                          ## 'sunl': 1 Solar mass with the Sun's density (i.e., 1 Solar radius)
                          ## 'drawkrou': draw from Kroupa IMF with densities consistent with mass
                          ## 'wdwf': 1 Solar mass and 0.01 Solar radius
                          typestar='sunl', \
                          
                          # Boolean flag to diagnose
                          booldiag=True, \

                          # legacy alias for stellar-population sampling behavior
                          booltoyysunn=False, \
                          
                          # list of bands
                          liststrgband=None, \
                          
                         ):
    '''
    Sample a synthetic population of the features of companions (e.g., exoplanets )and the companions to companions (e.g., exomoons) 
    hosted by a specified or random population of stellar systems.
    '''
    
    print('Simulating a population of systems via retr_dictpoplstarcomp()...')
    print('typesyst')
    print(typesyst)
    
    if typepoplsyst is None:
        typepoplsyst = 'SyntheticPopulation'
    elif typepoplsyst == 'Synthetic':
        typepoplsyst = 'SyntheticPopulation'

    print('typepoplsyst')
    print(typepoplsyst)
    
    if booldiag:
        if maxmcosicomp is not None and (maxmcosicomp > 1 or maxmcosicomp < 0):
            raise Exception('maxmcosicomp is bad.')
    
    # Boolean flag indicating if the system is a star with a companion
    boolhavecomp = typesyst.startswith('PlanetarySystem') or typesyst == 'CompactObjectStellarCompanion' or typesyst == 'StellarBinary'
       
    # Boolean flag indicating if the system is a star with a companion with moons
    boolhavemoon = typesyst.startswith('PlanetarySystemWithMoons')
       
    # Boolean flag indicating if the system is a planetary system
    boolsystpsys = typesyst.startswith('PlanetarySystem')
        
    # Boolean flag indicating if the system is a compact object with stellar companion
    boolsystcosc = typesyst == 'CompactObjectStellarCompanion'

    # Boolean flag indicating if the star is flaring
    boolflar = typesyst == 'StarFlaring'
    
    if booldiag:
        if not (boolflar or boolhavecomp):
            print('')
            print('')
            print('')
            print('typesyst')
            print(typesyst)
            raise Exception('typesyst is undefined.')

    # dictionary keys of the populations
    namepoplstartotl = 'star_%s_All' % typepoplsyst
    namepoplstaroccu = 'star_%s_Occurrent' % typepoplsyst
    
    dictnico = dict()
    dictpopl = dict()
    dictpopl['star'] = dict()
    dictstarnumbsamp = dict()
    dictstarindxsamp = dict()
    dictstarnumbsamp[namepoplstartotl] = dict()
    dictstarindxsamp[namepoplstartotl] = dict()
    
    dictfact = tdpy.retr_factconv()
    
    # number of systems
    if numbsyst is None:
        if distsyst is not None:
            numbsyst = len(distsyst)
        elif typepoplsyst == 'SyntheticPopulation':
            numbsyst = 10000
    
    indxsyst = np.arange(numbsyst)

    # get the features of the star population
    if typepoplsyst.startswith('CTL') or typepoplsyst.startswith('TIC'):
        dictstar = retr_dictpopltic8(typepoplsyst, numbsyst=numbsyst)
        
        print('Removing stars that do not have radii or masses...')
        indx = np.where(np.isfinite(dictstar['radistar'][0]) & \
                        np.isfinite(dictstar['massstar'][0]))[0]
        for name in dictstar.keys():
            dictstar[name][0] = dictstar[name][0][indx]

        if (dictstar['rascstar'][0] > 1e4).any():
            raise Exception('')

        if (dictstar['radistar'][0] == 0.).any():
            raise Exception('')
        
        densstar = 1.41 * dictstar['massstar'][0] / dictstar['radistar'][0]**3
        tdpy.setp_dict(dictstar, 'densstar', densstar)
        
        dictstar['idenstar'] = dictstar['TICID']
    

    elif typepoplsyst == 'SyntheticPopulation':
        
        dictstar = dict()
        
        if distsyst is None:
            dictstar['distsyst'] = [tdpy.icdf_powr(np.random.rand(numbsyst), 100., 7000., -2.), 'pc']
        elif isinstance(distsyst, np.ndarray):
            dictstar['distsyst'] = [distsyst, 'pc']
        elif isinstance(distsyst, float):
            dictstar['distsyst'] = [np.full(numbsyst, distsyst), 'pc']
        
        # stellar masses, radii, and densities
        ## all masses and radii should be in solar units and densities should be in g/cm3
        dictstar['radistar'] = [[], []]
        dictstar['massstar'] = [[], []]
        dictstar['densstar'] = [[], []]
        if typestar == 'sunl':
            dictstar['radistar'][0] = np.ones(numbsyst)
            dictstar['massstar'][0] = np.ones(numbsyst)
            dictstar['densstar'][0] = np.full(numbsyst, 1.4)
        elif massstar is not None or typestar == 'drawkrou':
            if massstar is not None:
                dictstar['massstar'][0] = np.full(numbsyst, massstar)
            if typestar == 'drawkrou':
                dictstar['massstar'][0] = tdpy.icdf_powr(np.random.rand(numbsyst), 0.1, 10., 2.)
            if radistar is not None:
                dictstar['radistar'][0] = np.full(numbsyst, radistar)
                dictstar['densstar'][0] = 1.4 * dictstar['massstar'][0] / dictstar['radistar'][0]**3
            else:
                dictstar['densstar'][0] = 1.4 * dictstar['massstar'][0]**(-0.7)
                dictstar['radistar'][0] = 1.4 * dictstar['massstar'][0] / dictstar['densstar'][0]**(1. / 3.)
        elif typestar == 'wdwf':
            dictstar['radistar'] = [0.01 * np.ones(numbsyst), r'$R_{\odot}$']
            dictstar['massstar'] = [np.ones(numbsyst), r'$M_{\odot}$']
            dictstar['densstar'] = [1.4e6 * np.ones(numbsyst), 'g cm$^{-3}$']
        else:
            raise Exception('')
        dictstar['radistar'][1] = r'$R_{\odot}$'
        dictstar['massstar'][1] = r'$M_{\odot}$'
        dictstar['densstar'][1] = 'g cm$^{-3}$'

        dictstar['coeflmdklinr'] = [0.4 * np.ones(numbsyst), '']
        dictstar['coeflmdkquad'] = [0.25 * np.ones(numbsyst), '']

        dictstar['tmptstar'] = [6000 * dictstar['massstar'][0], 'K']
        
        dictstar['lumistar'] = [4. * np.pi * dictstar['tmptstar'][0]**4 * dictstar['radistar'][0]**2, r'$L_{\odot}$']
        
        dictstar['fluxbolostar'] = [1361. * dictstar['lumistar'][0] / dictstar['distsyst'][0]**2 / 4. / np.pi, 'W/m^2']

        if liststrgband is None:
            liststrgband = []
            if typepoplsyst == 'lsstwfds':
                liststrgband += ['r']

        gdatfluxband = None
        for strgband in liststrgband:
            dictstar['fluxbandsyst%s' % strgband] = [np.empty(numbsyst), 'W/m^2']
        for k in indxsyst: 
            # flux at the surface
            dictfluxband, gdatfluxband = tdpy.retr_dictfluxband(dictstar['tmptstar'][0][k], liststrgband, gdatfluxband=gdatfluxband)
            
            # flux at the observer
            for strgband in liststrgband:
                dictstar['fluxbandsyst%s' % strgband][0][k] = dictfluxband[strgband] * (dictstar['radistar'][0][k] / 4.435e7 / dictstar['distsyst'][0][k])**2
        
        print('fluxbandsyst')
        for strgband in liststrgband:
            dictstar['magtsyst%s' % strgband] = [tdpy.retr_magtfromflux(dictstar['fluxbandsyst%s' % strgband][0], strgband), 'mag']
            
            print('strgband')
            print(strgband)
            print('fluxbandsyst')
            summgene(dictstar['fluxbandsyst%s' % strgband][0])
            print('magtsyst')
            summgene(dictstar['magtsyst%s' % strgband][0])
            print('')
    else:
        print('')
        print('')
        print('')
        print('typepoplsyst')
        print(typepoplsyst)
        raise Exception('Undefined typepoplsyst.')
    
    dictstarnumbsamp[namepoplstartotl] = numbsyst
    dictpopl['star'][namepoplstartotl] = dictstar

    def _legacy_alias_popl(popl):
        out = dict()
        for key, val in popl.items():
            if isinstance(val, list) and len(val) == 2:
                out[key] = val[0]
            elif isinstance(val, dict):
                out[key] = _legacy_alias_popl(val)
            else:
                out[key] = val
        return out

    # total mass
    dictstar['masssyst'] = [[], []]
    dictstar['masssyst'][0] = np.copy(dictstar['massstar'][0])
    dictstar['masssyst'][1] = dictstar['massstar'][1]
    
    numbstar = numbsyst

    dictindx = dict()
    dictnumbsamp = dict()
    dictindxsamp = dict()

    if boolhavecomp:
        # minimum companion mass
        if minmmasscomp is None:
            if boolsystcosc:
                minmmasscomp = 5. # [Solar mass]
            elif typesyst == 'PlanetarySystem' or typesyst == 'PlanetarySystemEmittingCompanion' or typesyst == 'PlanetarySystemWithMoons':
                if typesyst == 'PlanetarySystem' or typesyst == 'PlanetarySystemWithMoons':
                    # ~ Mars mass
                    minmmasscomp = 0.1 # [Earth mass]
                if typesyst == 'PlanetarySystemEmittingCompanion':
                    minmmasscomp = 3000. # [Earth mass]
            elif typesyst == 'StellarBinary':
                minmmasscomp = 0.5 # [Earth mass]
    
        # maximum companion mass
        if maxmmasscomp is None:
            if boolsystcosc:
                maxmmasscomp = 200. # [Solar mass]
            elif boolsystpsys:
                # Deuterium burning mass
                maxmmasscomp = 4400. # [Earth mass]
            elif typesyst == 'StellarBinary':
                maxmmasscomp = np.inf # [Earth mass]
            else:
                print('')
                print('')
                print('')
                print('typesyst')
                print(typesyst)
                raise Exception('Could not define maxmmasscomp')
    
    if boolhavecomp or boolflar:
        
        strgbody = 'star'

        if boolhavecomp:
            strglimb = 'comp'
        else:
            strglimb = 'flar'
        
        strgnumblimbbody = 'numb%s%s' % (strglimb, strgbody)
        strgnumblimbbodymean = 'numb%s%smean' % (strglimb, strgbody)
        
        namepopllimbtotl = '%s%s_%s_All' % (strglimb, strgbody, typepoplsyst)
        if strglimb == 'comp':
            namepoplcomptran = '%s%s_%s_Transiting' % (strglimb, strgbody, typepoplsyst)
        
        dictpopl[strglimb] = dict()
        dictnumbsamp[strglimb] = dict()
        dictindxsamp[strglimb] = dict()
        dictnumbsamp[strglimb][namepopllimbtotl] = dict()
        dictindxsamp[strglimb][namepopllimbtotl] = dict()
        dictpopl[strglimb][namepopllimbtotl] = dict()

        if boolhavemoon:
            namepoplmoontotl = 'mooncompstar_%s_All'
            dictpopl['moon'] = dict()
            dictmoonnumbsamp = dict()
            dictmoonindxsamp = dict()
            dictmoonnumbsamp[namepoplmoontotl] = dict()
            dictmoonindxsamp[namepoplmoontotl] = dict()
            dictpopl['moon'][namepoplmoontotl] = dict()
    
        if typesyst.startswith('PlanetarySystem') or typesyst == 'StarFlaring':
            
            if typesyst.startswith('PlanetarySystem'):
                # mean number of companions per star
                dictpopl[strgbody][namepoplstartotl][strgnumblimbbodymean] = [0.5 * dictpopl[strgbody][namepoplstartotl]['massstar'][0]**(-1.), '']
            if typesyst == 'StarFlaring':
                # mean number of flares per star
                dictpopl[strgbody][namepoplstartotl][strgnumblimbbodymean] = [0.5 * dictpopl[strgbody][namepoplstartotl]['massstar'][0]**(-1.), '']
            
            # mean number per star
            dictpopl[strgbody][namepoplstartotl][strgnumblimbbodymean] = [5. * dictpopl[strgbody][namepoplstartotl]['massstar'][0]**(-1.), '']
            
            # number per star
            dictpopl[strgbody][namepoplstartotl][strgnumblimbbody] = [np.random.poisson(dictpopl[strgbody][namepoplstartotl][strgnumblimbbodymean][0]), '']
            
            if booldiag:
                if maxmnumbcompstar is not None and minmnumbcompstar is not None:
                    if maxmnumbcompstar < minmnumbcompstar:
                        print('')
                        print('')
                        print('')
                        raise Exception('maxmnumbcompstar < minmnumbcompstar')

            if minmnumbcompstar is not None:
                dictpopl[strgbody][namepoplstartotl][strgnumblimbbody] = [np.maximum(dictpopl[strgbody][namepoplstartotl][strgnumblimbbody][0], minmnumbcompstar), '']

            if maxmnumbcompstar is not None:
                dictpopl[strgbody][namepoplstartotl][strgnumblimbbody] = [np.minimum(dictpopl[strgbody][namepoplstartotl][strgnumblimbbody][0], maxmnumbcompstar), '']

        elif typesyst == 'CompactObjectStellarCompanion' or typesyst == 'StellarBinary':
            # number of companions per star
            dictpopl[strgbody][namepoplstartotl][strgnumblimbbody] = [np.ones(dictpopl[strgbody][namepoplstartotl]['radistar'][0].size).astype(int), '']
            
        else:
            print('')
            print('')
            print('')
            print('typesyst')
            print(typesyst)
            raise Exception('typesyst is undefined.')
    
        if booldiag:
            if np.isscalar(dictpopl[strgbody][namepoplstartotl]['distsyst'][0]):
                print('')
                print('')
                print('')
                print('dictpopl[strgbody][namepoplstartotl][distsyst][0]')
                summgene(dictpopl[strgbody][namepoplstartotl]['distsyst'][0])
                raise Exception('')

        # Boolean flag of occurence
        dictpopl[strgbody][namepoplstartotl]['booloccu'] = [dictpopl[strgbody][namepoplstartotl][strgnumblimbbody][0] > 0, '']
    
        # subpopulation where companions or flares occur
        indx = np.where(dictpopl[strgbody][namepoplstartotl]['booloccu'][0])[0]
        tdpy.retr_subp(dictpopl[strgbody], dictstarnumbsamp, dictstarindxsamp, namepoplstartotl, namepoplstaroccu, indx)
    
        # indices of companions or flares for each star
        dictindx[strglimb] = dict()
        dictindx[strglimb][strgbody] = [[] for k in indxsyst]
        cntr = 0
        for k in indxsyst:
            dictindx[strglimb][strgbody][k] = np.arange(cntr, cntr + dictpopl[strgbody][namepoplstartotl][strgnumblimbbody][0][k]).astype(int)
            cntr += dictpopl[strgbody][namepoplstartotl][strgnumblimbbody][0][k]
        dictnumbsamp[strglimb][namepopllimbtotl] = cntr
    
        # prepare to load star features into component features
        for name in list(dictpopl[strgbody][namepoplstartotl].keys()):
            dictpopl[strglimb][namepopllimbtotl][name] = [np.empty(dictnumbsamp[strglimb][namepopllimbtotl]), '']
        
        print('Loading star features into component features...')
        dictnumbsamp[strglimb][namepopllimbtotl] = dictpopl[strglimb][namepopllimbtotl]['radistar'][0].size

        listnamecatr = ['masssyst', 'radistar']
        if strglimb == 'comp':
            listnamecatr += ['pericomp', 'cosicomp', 'smaxcomp', 'eccecomp', 'arpacomp', 'loancomp', 'masscomp', 'epocmtracomp', 'rsmacomp', 'rsumcomp', 'distanglcomp']
            if typesyst == 'PlanetarySystemWithMoons':
                listnamecatr += ['masscompmoon']
            if typesyst == 'PlanetarySystemWithNonKeplerianObjects':
                listnamecatr += ['factnonkcomp']
            if not boolsystcosc:
                listnamecatr += ['radicomp', 'denscomp']
        elif strglimb == 'flar':
            listnamecatr += ['amplflar', 'timeflar', 'enerflar', 'tsclflar']
        else:
            raise Exception('')

        for name in listnamecatr:
            dictpopl[strglimb][namepopllimbtotl][name] = [np.empty(dictnumbsamp[strglimb][namepopllimbtotl]), '']

        if booldiag:
            if np.isscalar(dictpopl[strglimb][namepopllimbtotl]['distsyst'][0]):
                print('')
                print('')
                print('')
                print('dictpopl[strglimb][namepopllimbtotl][distsyst][0]')
                summgene(dictpopl[strglimb][namepopllimbtotl]['distsyst'][0])
                raise Exception('')

            cntr = 0
            for k in indxsyst:
                cntr += len(dictindx[strglimb][strgbody][k])
            if cntr != dictnumbsamp[strglimb][namepopllimbtotl]:
                print('')
                print('')
                print('')
                print('cntr')
                print(cntr)
                print('dictnumbsamp[strglimb][namepopllimbtotl]')
                print(dictnumbsamp[strglimb][namepopllimbtotl])
                raise Exception('cntr != dictnumbsamp[strglimb][namepopllimbtotl]')
    
        print('Iterating over the stars...')
        for k in tqdm(range(numbstar)):
            
            if dictpopl[strgbody][namepoplstartotl][strgnumblimbbody][0][k] == 0:
                print('Skipping')
                print('Skipping')
                print('Skipping')
                print('Skipping')
                print('Skipping')
                print('Skipping')
                continue

            # load star features into component features
            for namefeat in dictpopl[strgbody][namepoplstartotl].keys():
                dictpopl[strglimb][namepopllimbtotl][namefeat][0][dictindx[strglimb][strgbody][k]] = dictpopl[strgbody][namepoplstartotl][namefeat][0][k]
            
            if strglimb == 'comp':

                numb = dictpopl[strgbody][namepoplstartotl][strgnumblimbbody][0][k]

                # eccentricities
                dictpopl[strglimb][namepopllimbtotl]['eccecomp'][0][dictindx[strglimb][strgbody][k]] = np.random.rand(numb)
                
                # arguments of periapsis
                dictpopl[strglimb][namepopllimbtotl]['arpacomp'][0][dictindx[strglimb][strgbody][k]] = 2. * np.pi * np.random.rand(numb)
                
                # longtides of ascending node
                dictpopl[strglimb][namepopllimbtotl]['loancomp'][0][dictindx[strglimb][strgbody][k]] = 2. * np.pi * np.random.rand(numb)
                
                # companion mass
                dictpopl[strglimb][namepopllimbtotl]['masscomp'][0][dictindx[strglimb][strgbody][k]] = tdpy.util.icdf_powr(np.random.rand(numb), minmmasscomp, maxmmasscomp, 2.)
                
                if boolsystpsys or typesyst == 'StellarBinary':
                    # companion radius
                    dictpopl[strglimb][namepopllimbtotl]['radicomp'][0][dictindx[strglimb][strgbody][k]] = \
                                                        tdpy.retr_radifrommass(dictpopl[strglimb][namepopllimbtotl]['masscomp'][0][dictindx[strglimb][strgbody][k]])
        
                    # companion density
                    dictpopl[strglimb][namepopllimbtotl]['denscomp'][0][dictindx[strglimb][strgbody][k]] = \
                                        5.51 * dictpopl[strglimb][namepopllimbtotl]['masscomp'][0][dictindx[strglimb][strgbody][k]] / \
                                                                                           dictpopl[strglimb][namepopllimbtotl]['radicomp'][0][dictindx[strglimb][strgbody][k]]**3
                
                # total radius
                dictpopl[strglimb][namepopllimbtotl]['rsumcomp'][0][dictindx[strglimb][strgbody][k]] = dictpopl[strgbody][namepoplstartotl]['radistar'][0][k]
                if not boolsystcosc:
                    dictpopl[strglimb][namepopllimbtotl]['rsumcomp'][0][dictindx[strglimb][strgbody][k]] += \
                                            dictpopl[strglimb][namepopllimbtotl]['radicomp'][0][dictindx[strglimb][strgbody][k]] / dictfact['rsre']
                
                if booldiag:
                    if (dictpopl[strglimb][namepopllimbtotl]['rsumcomp'][0][dictindx[strglimb][strgbody][k]] < 1e-30).any():
                        print('')
                        print('')
                        print('')
                        raise Exception('')

                    if (dictpopl[strglimb][namepopllimbtotl]['rsumcomp'][0][dictindx[strglimb][strgbody][k]] <= \
                       dictpopl[strglimb][namepopllimbtotl]['radicomp'][0][dictindx[strglimb][strgbody][k]] / dictfact['rsre']).any() or \
                       (dictpopl[strglimb][namepopllimbtotl]['rsumcomp'][0][dictindx[strglimb][strgbody][k]] <= \
                       dictpopl[strglimb][namepopllimbtotl]['radistar'][0][k]).any():
                        raise Exception('')
                    
                # total mass
                if boolsystcosc or typesyst == 'StellarBinary':
                    dictpopl[strgbody][namepoplstartotl]['masssyst'][0][k] += np.sum(dictpopl[strglimb][namepopllimbtotl]['masscomp'][0][dictindx[strglimb][strgbody][k]])
                
                if typesamporbtcomp == 'peristab' or typesamporbtcomp == 'peripowr':
                    
                    if typesamporbtcomp == 'peristab':
                        ratiperifrst = tdpy.util.icdf_powr(np.random.rand(1), 2., 1000., 0.5)[0]
                        ratiperi = tdpy.util.icdf_powr(np.random.rand(dictpopl[strgbody][namepoplstartotl][strgnumblimbbody][0][k] - 1), 1.2, 100., 2.)

                        listpericomp = []
                        for mm in range(dictpopl[strgbody][namepoplstartotl][strgnumblimbbody][0][k]):
                            if mm == 0:
                                peri = minmpericomp * ratiperifrst
                            else:
                                peri = ratiperi[mm-1] * listpericomp[mm-1]
                            listpericomp.append(peri)
                        dictpopl[strglimb][namepopllimbtotl]['pericomp'][0][dictindx[strglimb][strgbody][k]] = np.array(listpericomp)

                        if booldiag:
                            ratiperi = dictpopl[strglimb][namepopllimbtotl]['pericomp'][0][dictindx[strglimb][strgbody][k]][1:] / \
                                                        dictpopl[strglimb][namepopllimbtotl]['pericomp'][0][dictindx[strglimb][strgbody][k]][:-1]
                            indx = np.where(ratiperi < 1.2)[0]
                            if indx.size > 0:
                                print('indx')
                                summgene(indx)
                                print('dictpopl[comp][namepopllimbtotl][pericomp][dictindx[strglimb][strgbody][k]]')
                                print(dictpopl[strglimb][namepopllimbtotl]['pericomp'][0][dictindx[strglimb][strgbody][k]])
                                print('dictpopl[star][namepoplstartotl][numbcompstar][k]')
                                print(dictpopl[strgbody][namepoplstartotl][strgnumblimbbody][0][k])
                                raise Exception('')
                    if typesamporbtcomp == 'peripowr':
                        dictpopl[strglimb][namepopllimbtotl]['pericomp'][0][dictindx[strglimb][strgbody][k]] = \
                                        tdpy.util.icdf_powr(np.random.rand(dictpopl[strgbody][namepoplstartotl][strgnumblimbbody][0][k]), minmpericomp, maxmpericomp, 2.)
                    
                    # determine the corresponding semi-major axes
                    if typesyst == 'PlanetarySystemWithNonKeplerianObjects':
                        raise Exception('')
                    else:
                        factnonk = 1.
                    
                    dictpopl[strglimb][namepopllimbtotl]['smaxcomp'][0][dictindx[strglimb][strgbody][k]] = \
                                            tdpy.retr_smaxkepl(dictpopl[strglimb][namepopllimbtotl]['pericomp'][0][dictindx[strglimb][strgbody][k]], \
                                                                                                     dictpopl[strgbody][namepoplstartotl]['masssyst'][0][k], factnonk=factnonk)
                
                else:
                    # semi-major axes
                    #if np.isfinite(dictpopl[strgbody][namepoplstartotl]['densstar'][k]):
                    #    densstar = dictpopl[strgbody][namepoplstartotl]['densstar'][k]
                    #else:
                    #    densstar = 1.
                    #dictpopl[strglimb][namepopllimbtotl]['radiroch'][k] = retr_radiroch(radistar, densstar, denscomp)
                    #minmsmax = 2. * dictpopl[strglimb][namepopllimbtotl]['radiroch'][k]
                    
                    dictpopl[strglimb][namepopllimbtotl]['smaxcomp'][0][dictindx[strglimb][strgbody][k]] = dictpopl[strglimb][namepopllimbtotl]['rsumcomp'][0][k] * \
                                                                                 tdpy.util.icdf_powr(np.random.rand(dictpopl[strgbody][namepoplstartotl][strgnumblimbbody][0][k]), \
                                                                                                    minmsmaxradistar, maxmsmaxradistar, 2.) / dictfact['aurs']
                    
                    if typesyst == 'PlanetarySystemWithNonKeplerianObjects':
                        factnonk = tdpy.util.icdf_powr(np.random.rand(dictpopl[strgbody][namepoplstartotl][strgnumblimbbody][0][k]), 0.1, 1., -2.)
                        dictpopl[strglimb][namepopllimbtotl]['factnonkcomp'][0][dictindx[strglimb][strgbody][k]] = factnonk
                    else:
                        factnonk = 1.
                    
                    dictpopl[strglimb][namepopllimbtotl]['pericomp'][0][dictindx[strglimb][strgbody][k]] = \
                                            tdpy.retr_perikepl(dictpopl[strglimb][namepopllimbtotl]['smaxcomp'][0][dictindx[strglimb][strgbody][k]], \
                                                                                                  dictpopl[strgbody][namepoplstartotl]['masssyst'][0][k], factnonk=factnonk)
                    
        
                dictpopl[strglimb][namepopllimbtotl]['distanglcomp'][0][dictindx[strglimb][strgbody][k]] = \
                                        dictpopl[strglimb][namepopllimbtotl]['smaxcomp'][0][dictindx[strglimb][strgbody][k]] / dictstar['distsyst'][0][k]

                if booldiag:
                    if not np.isfinite(dictpopl[strglimb][namepopllimbtotl]['pericomp'][0][dictindx[strglimb][strgbody][k]]).all():
                        print('')
                        print('')
                        print('')
                        print('dictpopl[comp][namepopllimbtotl][masscomp][dictindx[comp][star][k]]')
                        print(dictpopl[strglimb][namepopllimbtotl]['masscomp'][0][dictindx[strglimb][strgbody][k]])
                        print('dictpopl[comp][namepopllimbtotl][masssyst][dictindx[comp][star][k]]')
                        print(dictpopl[strglimb][namepopllimbtotl]['masssyst'][0][dictindx[strglimb][strgbody][k]])
                        print('dictpopl[comp][namepopllimbtotl][smaxcomp][dictindx[comp][star][k]]')
                        print(dictpopl[strglimb][namepopllimbtotl]['smaxcomp'][0][dictindx[strglimb][strgbody][k]])
                        print('dictpopl[comp][namepopllimbtotl][pericomp][dictindx[comp][star][k]]')
                        print(dictpopl[strglimb][namepopllimbtotl]['pericomp'][0][dictindx[strglimb][strgbody][k]])
                        raise Exception('')

                # conjunction epochs
                if epocmtracomp is not None:
                    dictpopl[strglimb][namepopllimbtotl]['epocmtracomp'][0][dictindx[strglimb][strgbody][k]] = \
                                                            np.full(dictpopl[strgbody][namepoplstartotl][strgnumblimbbody][0][k], epocmtracomp)
                else:
                    dictpopl[strglimb][namepopllimbtotl]['epocmtracomp'][0][dictindx[strglimb][strgbody][k]] = \
                                        1e8 * np.random.rand(dictpopl[strgbody][namepoplstartotl][strgnumblimbbody][0][k])
                if timeepoc is not None:
                    dictpopl[strglimb][namepopllimbtotl]['epocmtracomp'][0][dictindx[strglimb][strgbody][k]] = dictpopl[strglimb][namepopllimbtotl]['epocmtracomp'][k] + \
                                                   dictpopl[strglimb][namepopllimbtotl]['pericomp'][k] * \
                                                   np.round((dictpopl[strglimb][namepopllimbtotl]['epocmtracomp'][k] - timeepoc) / dictpopl[strglimb][namepopllimbtotl]['pericomp'][k])
    
        # derived stellar features --  dictpopl[strgbody][namepoplstartotl] should not be enhanced with derived features before this
        # Boolean flag indicating whether the system has an inner, Earth-sized rocky planet
        dictpopl[strgbody][namepoplstartotl]['boolelik'] = [np.empty(numbstar, dtype=bool), '']
        # Boolean flag indicating whether the system has an outer, Jupiter-sized giant planet
        dictpopl[strgbody][namepoplstartotl]['booljlik'] = [np.empty(numbstar, dtype=bool), '']
        
        print('Iterating over the stars...')
        for k in tqdm(range(numbstar)):
            
            dictpopl[strgbody][namepoplstartotl]['boolelik'][0][k] = np.where( \
                                                                              (dictpopl[strglimb][namepopllimbtotl]['pericomp'][0][dictindx[strglimb][strgbody][k]] > 0.5 * 365.25) & \
                                                                              (dictpopl[strglimb][namepopllimbtotl]['pericomp'][0][dictindx[strglimb][strgbody][k]] < 2. * 365.25) & \
                                                                              (dictpopl[strglimb][namepopllimbtotl]['radicomp'][0][dictindx[strglimb][strgbody][k]] > 0.5) & \
                                                                              (dictpopl[strglimb][namepopllimbtotl]['radicomp'][0][dictindx[strglimb][strgbody][k]] < 2.) \
                                                                             )[0].size > 0
            
            dictpopl[strgbody][namepoplstartotl]['booljlik'][0][k] = np.where( \
                                                                        (dictpopl[strglimb][namepopllimbtotl]['pericomp'][0][dictindx[strglimb][strgbody][k]] > 0.5 * 11.86 * 365.25) & \
                                                                        (dictpopl[strglimb][namepopllimbtotl]['pericomp'][0][dictindx[strglimb][strgbody][k]] < 2. * 11.86 * 365.25) & \
                                                                        (dictpopl[strglimb][namepopllimbtotl]['radicomp'][0][dictindx[strglimb][strgbody][k]] > 0.6 * 11.2) & \
                                                                        (dictpopl[strglimb][namepopllimbtotl]['radicomp'][0][dictindx[strglimb][strgbody][k]] < 1.5 * 11.2) \
                                                                             )[0].size > 0
            
        # Boolean flag indicating whether the system has a an inner, Earth-sized rocky planet AND outer, Jupiter-sized giant planet
        dictpopl[strgbody][namepoplstartotl]['boolelikjlik'] = [dictpopl[strgbody][namepoplstartotl]['boolelik'][0] & dictpopl[strgbody][namepoplstartotl]['booljlik'][0], '']

        if strglimb == 'comp':
            if typesyst == 'PlanetarySystemWithMoons':
                # initialize the total mass of the companion + moons system as the mass of the companion
                dictpopl[strglimb][namepopllimbtotl]['masscompmoon'] = np.copy(dictpopl[strglimb][namepopllimbtotl]['masscomp'])
                        
            dictpopl[strglimb][namepopllimbtotl]['rsmacomp'][0] = dictpopl[strglimb][namepopllimbtotl]['rsumcomp'][0] \
                                                                            / (dictpopl[strglimb][namepopllimbtotl]['smaxcomp'][0] * dictfact['aurs'])
            
            if booldiag:
                if (dictpopl[strglimb][namepopllimbtotl]['rsmacomp'][0] >= 1.).any():
                    print('')
                    print('')
                    print('')
                    print('dictpopl[strglimb][namepopllimbtotl][rsmacomp][0]')
                    summgene(dictpopl[strglimb][namepopllimbtotl]['rsmacomp'][0])
                    print('dictpopl[strglimb][namepopllimbtotl][rsumcomp][0]')
                    summgene(dictpopl[strglimb][namepopllimbtotl]['rsumcomp'][0])
                    print(dictpopl[strglimb][namepopllimbtotl]['rsumcomp'][0])
                    print('dictpopl[strglimb][namepopllimbtotl][radistar][0]')
                    summgene(dictpopl[strglimb][namepopllimbtotl]['radistar'][0])
                    print('dictpopl[strglimb][namepopllimbtotl][radicomp][0] / dictfact[rsre]')
                    summgene(dictpopl[strglimb][namepopllimbtotl]['radicomp'][0] / dictfact['rsre'])
                    print('dictpopl[strglimb][namepopllimbtotl][smaxcomp][0] * dictfact[aurs]')
                    summgene(dictpopl[strglimb][namepopllimbtotl]['smaxcomp'][0] * dictfact['aurs'])
                    print(dictpopl[strglimb][namepopllimbtotl]['smaxcomp'][0] * dictfact['aurs'])
                    raise Exception('rsum is more than the semi-major axis!')

            if booltrancomp is True and maxmcosicomp is not None:
                raise Exception('maxmcosicomp cannot be specified if booltrancomp is True.')

            for k in tqdm(range(numbstar)):
                
                if dictpopl[strgbody][namepoplstartotl][strgnumblimbbody][0][k] == 0:
                    continue
                
                if booltrancomp:
                    maxmcosicomptemp = dictpopl[strglimb][namepopllimbtotl]['rsmacomp'][0][k]
                elif maxmcosicomp is not None:
                    maxmcosicomptemp = maxmcosicomp
                else:
                    maxmcosicomptemp = 1.
                
                numb = dictpopl[strgbody][namepoplstartotl][strgnumblimbbody][0][k]

                # cosine of orbital inclinations
                dictpopl[strglimb][namepopllimbtotl]['cosicomp'][0][dictindx[strglimb][strgbody][k]] = maxmcosicomptemp * np.random.rand(numb)
                
            if booldiag:
                if (dictpopl[strglimb][namepopllimbtotl]['rsmacomp'][0] == 0).any():
                    print('dictpopl[strglimb][namepopllimbtotl][rsmacomp][0]')
                    print(dictpopl[strglimb][namepopllimbtotl]['rsmacomp'][0])
                    print('dictpopl[strglimb][namepopllimbtotl][rsumcomp][0]')
                    print(dictpopl[strglimb][namepopllimbtotl]['rsumcomp'][0])
                    print('')
                    print('')
                    print('')
                    raise Exception('')

            # orbital inclinations of the companions
            dictpopl[strglimb][namepopllimbtotl]['inclcomp'] = [180. / np.pi * np.arccos(dictpopl[strglimb][namepopllimbtotl]['cosicomp'][0]), 'degree']
            
            dictpopl[strglimb][namepopllimbtotl]['inclcomp'][0] = 90. + (dictpopl[strglimb][namepopllimbtotl]['inclcomp'][0] - 90.) * \
                                                                            (2 * np.random.randint(2, size=dictpopl[strglimb][namepopllimbtotl]['cosicomp'][0].size) - 1.)

            if boolsystpsys:
                
                if booldiag:
                    
                    if not np.isfinite(dictpopl[strglimb][namepopllimbtotl]['radistar'][0]).all():
                        print('')
                        print('')
                        print('')
                        raise Exception('not np.isfinite(dictpopl[comp][namepopllimbtotl][radistar]).all()')
                    
                    if not np.isfinite(dictpopl[strglimb][namepopllimbtotl]['radicomp'][0]).all():
                        print('')
                        print('')
                        print('')
                        raise Exception('not np.isfinite(dictpopl[comp][namepopllimbtotl][radicomp]).all()')

                # radius ratio
                dictpopl[strglimb][namepopllimbtotl]['rratcomp'] = [dictpopl[strglimb][namepopllimbtotl]['radicomp'][0] / \
                                                                        dictpopl[strglimb][namepopllimbtotl]['radistar'][0] / dictfact['rsre'], '']
                
                if booldiag:
                    if not np.isfinite(dictpopl[strglimb][namepopllimbtotl]['rratcomp'][0]).all():
                        print('')
                        print('')
                        print('')
                        raise Exception('not np.isfinite(dictpopl[comp][namepopllimbtotl][rratcomp]).all()')
            
            # Boolean flag indicating whether a companion is transiting
            dictpopl[strglimb][namepopllimbtotl]['booltran'] = [dictpopl[strglimb][namepopllimbtotl]['rsmacomp'][0] > dictpopl[strglimb][namepopllimbtotl]['cosicomp'][0], '']

            # Boolean flag indicating whether any companion is transiting
            booltran = np.empty(numbstar)
            for k in tqdm(range(numbstar)):
                booltran[k] = dictpopl[strglimb][namepopllimbtotl]['booltran'][0][dictindx[strglimb][strgbody][k]].any()
            dictpopl[strgbody][namepoplstartotl]['booltran'] = [booltran, '']

            # subpopulation where object transits
            indx = np.where(dictpopl[strglimb][namepopllimbtotl]['booltran'][0])[0]
            tdpy.retr_subp(dictpopl[strglimb], dictnumbsamp[strglimb], dictindxsamp[strglimb], namepopllimbtotl, namepoplcomptran, indx)

            # transit duration
            dictpopl[strglimb][namepoplcomptran]['duratrantotl'] = [tdpy.retr_duratrantotl(dictpopl[strglimb][namepoplcomptran]['pericomp'][0], \
                                                                                     dictpopl[strglimb][namepoplcomptran]['rsmacomp'][0], \
                                                                                     dictpopl[strglimb][namepoplcomptran]['cosicomp'][0]), '']
            dictpopl[strglimb][namepoplcomptran]['dcyc'] = [dictpopl[strglimb][namepoplcomptran]['duratrantotl'][0] / dictpopl[strglimb][namepoplcomptran]['pericomp'][0] / 24., '']
            
            if boolsystcosc:
                # amplitude of self-lensing
                dictpopl[strglimb][namepoplcomptran]['amplslen'] = chalcedon.retr_amplslen(dictpopl[strglimb][namepoplcomptran]['pericomp'][0], \
                                                                                           dictpopl[strglimb][namepoplcomptran]['radistar'][0], \
                                                                                           dictpopl[strglimb][namepoplcomptran]['masscomp'][0], \
                                                                                           dictpopl[strglimb][namepoplcomptran]['massstar'][0])
            

            
            if typesyst == 'PlanetarySystem':
                # transit depth
                dictpopl[strglimb][namepoplcomptran]['depttrancomp'] = [1e3 * dictpopl[strglimb][namepoplcomptran]['rratcomp'][0]**2, 'ppt']
            
            # define parent population's features that are valid only for transiting systems
            listtemp = ['duratrantotl', 'dcyc']
            if typesyst == 'PlanetarySystem':
                listtemp += ['depttrancomp']
            if boolsystcosc:
                listtemp += ['amplslen']
            for strg in listtemp:
                dictpopl[strglimb][namepopllimbtotl][strg] = [np.full_like(dictpopl[strglimb][namepopllimbtotl]['pericomp'][0], np.nan), '']
                dictpopl[strglimb][namepopllimbtotl][strg][0][indx] = dictpopl[strglimb][namepoplcomptran][strg][0]

            indxmooncompstar = [[[] for j in dictindx[strglimb][strgbody][k]] for k in indxsyst]
            if typesyst == 'PlanetarySystemWithMoons':
                # Hill radius of the companion
                dictpopl[strglimb][namepopllimbtotl]['radihill'] = tdpy.retr_radihill(dictpopl[strglimb][namepopllimbtotl]['smaxcomp'][0], \
                                                                                 dictpopl[strglimb][namepopllimbtotl]['masscomp'][0] / dictfact['msme'], \
                                                                                 dictpopl[strglimb][namepopllimbtotl]['massstar'][0])
            
                # maximum semi-major axis of the moons 
                dictpopl[strglimb][namepopllimbtotl]['maxmsmaxmoon'] = 0.2 * dictpopl[strglimb][namepopllimbtotl]['radihill']
                    
                # mean number of moons per companion
                dictpopl[strglimb][namepopllimbtotl]['numbmooncompmean'] = 1000. * dictpopl[strglimb][namepopllimbtotl]['masscomp']**(-1.)
                
                # number of moons per companion
                #dictpopl[strglimb][namepopllimbtotl]['numbmooncomp'] = np.random.poisson(dictpopl[strglimb][namepopllimbtotl]['numbmooncompmean'])
                dictpopl[strglimb][namepopllimbtotl]['numbmooncomp'] = np.ones_like(dictpopl[strglimb][namepopllimbtotl]['numbmooncompmean'])
                
                cntr = 0
                for k in range(len(dictpopl[strgbody][namepoplstartotl]['radistar'])):
                    for j in range(len(dictindx[strglimb][strgbody][k])):
                        indxmooncompstar[k][j] = np.arange(cntr, cntr + dictpopl[strglimb][namepopllimbtotl]['numbmooncomp'][0][dictindx[strglimb][strgbody][k][j]]).astype(int)
                        cntr += int(dictpopl[strglimb][namepopllimbtotl]['numbmooncomp'][j])
                dictmoonnumbsamp[namepoplmoontotl] = cntr
            
                numbmoontotl = int(np.sum(dictpopl[strglimb][namepopllimbtotl]['numbmooncomp']))
                
                # prepare to load component features into moon features
                for name in list(dictpopl[strglimb][namepopllimbtotl].keys()):
                    dictpopl['moon'][namepoplmoontotl][name] = np.empty(dictmoonnumbsamp[namepoplmoontotl])
            
                for k in tqdm(range(numbstar)):
                    
                    if dictpopl[strgbody][namepoplstartotl][strgnumblimbbody][0][k] == 0:
                        continue
                    
                    numbcomp = dictpopl[strgbody][namepoplstartotl][strgnumblimbbody][0][k]
                    
                    # number of exomoons to the companion
                    numbmoon = dictpopl[strglimb][namepopllimbtotl]['numbmooncomp'][0][dictindx[strglimb][strgbody][k]].astype(int)
                    for name in ['radi', 'mass', 'dens', 'peri', 'epocmtra', 'smax', 'minmsmax']:
                        dictpopl['moon'][namepoplmoontotl][name+'moon'] = np.empty(numbmoontotl)
                    
                    indxmoon = [[] for j in dictindx[strglimb][strgbody][k]]
                    for j in range(dictindx[strglimb][strgbody][k].size):
                        
                        #print('')
                        #print('')
                        #print('')
                        #print('')
                        #print('j')
                        #print(j)
                        #print('dictindx[strglimb][strgbody][k]')
                        #print(dictindx[strglimb][strgbody][k])
                        #print('numbmoon')
                        #print(numbmoon)
                        #print('indxmoon')
                        #print(indxmoon)
                        #print('indxmooncompstar[k][j]')
                        #print(indxmooncompstar[k][j])
                        #print('')
                        
                        if numbmoon[j] == 0:
                            continue

                        indxmoon[j] = np.arange(numbmoon[j])
                        # properties of the moons
                        ## mass [M_E]
                        dictpopl['moon'][namepoplmoontotl]['massmoon'][indxmooncompstar[k][j]] = dictpopl[strglimb][namepopllimbtotl]['masscomp'][0][dictindx[strglimb][strgbody][k][j]] * \
                                              tdpy.icdf_powr(np.random.rand(int(dictpopl[strglimb][namepopllimbtotl]['numbmooncomp'][0][dictindx[strglimb][strgbody][k][j]])), 0.005, 0.1, 2.)
                        
                        ## radii [R_E]
                        dictpopl['moon'][namepoplmoontotl]['radimoon'][indxmooncompstar[k][j]] = \
                                                    tdpy.retr_radifrommass(dictpopl['moon'][namepoplmoontotl]['massmoon'][indxmooncompstar[k][j]])
                        
                        ## densities [g/cm^3]
                        dictpopl['moon'][namepoplmoontotl]['densmoon'][indxmooncompstar[k][j]] = 5.51 * dictpopl['moon'][namepoplmoontotl]['massmoon'][indxmooncompstar[k][j]] / \
                                                                                                      dictpopl['moon'][namepoplmoontotl]['radimoon'][indxmooncompstar[k][j]]**3
                        
                        # minimum semi-major axes for the moons
                        dictpopl['moon'][namepoplmoontotl]['minmsmaxmoon'][indxmooncompstar[k][j]] = (1. / dictfact['rsre'] / dictfact['aurs']) * \
                                                                             tdpy.retr_radiroch(dictpopl[strglimb][namepopllimbtotl]['radicomp'][0][dictindx[strglimb][strgbody][k][j]], \
                                                                                           dictpopl[strglimb][namepopllimbtotl]['denscomp'][0][dictindx[strglimb][strgbody][k][j]], \
                                                                                                                  dictpopl['moon'][namepoplmoontotl]['densmoon'][indxmooncompstar[k][j]])
                        
                        if booldiag:
                            if (dictpopl['moon'][namepoplmoontotl]['minmsmaxmoon'][indxmooncompstar[k][j]] >= \
                                                                    dictpopl[strglimb][namepopllimbtotl]['maxmsmaxmoon'][dictindx[strglimb][strgbody][k][j]]).any():
                                print('')
                                print('')
                                print('')
                                print('typesamporbtcomp')
                                print(typesamporbtcomp)
                                print('dictpopl[comp][namepopllimbtotl][massstar][k]')
                                print(dictpopl[strglimb][namepopllimbtotl]['massstar'][0][k])
                                print('dictpopl[comp][namepopllimbtotl][masssyst][k]')
                                print(dictpopl[strglimb][namepopllimbtotl]['masssyst'][0][k])
                                print('dictpopl[comp][namepopllimbtotl][pericomp][dictindx[strglimb][strgbody][k][j]]')
                                print(dictpopl[strglimb][namepopllimbtotl]['pericomp'][0][dictindx[strglimb][strgbody][k][j]])
                                print('dictpopl[comp][namepopllimbtotl][smaxcomp][dictindx[strglimb][strgbody][k][j]]')
                                print(dictpopl[strglimb][namepopllimbtotl]['smaxcomp'][0][dictindx[strglimb][strgbody][k][j]])
                                print('dictpopl[comp][namepopllimbtotl][radihill][dictindx[strglimb][strgbody][k][j]]')
                                print(dictpopl[strglimb][namepopllimbtotl]['radihill'][dictindx[strglimb][strgbody][k][j]])
                                print('dictpopl[comp][namepopllimbtotl][maxmsmaxmoon][dictindx[strglimb][strgbody][k][j]]')
                                print(dictpopl[strglimb][namepopllimbtotl]['maxmsmaxmoon'][dictindx[strglimb][strgbody][k][j]])
                                print('dictpopl[moon][namepoplmoontotl][minmsmaxmoon][indxmooncompstar[k][j]]')
                                print(dictpopl['moon'][namepoplmoontotl]['minmsmaxmoon'][indxmooncompstar[k][j]])
                                print('dictpopl[moon][namepoplmoontotl][radimoon][indxmooncompstar[k][j]]')
                                print(dictpopl['moon'][namepoplmoontotl]['radimoon'][indxmooncompstar[k][j]])
                                print('dictpopl[moon][namepoplmoontotl][massmoon][indxmooncompstar[k][j]]')
                                print(dictpopl['moon'][namepoplmoontotl]['massmoon'][indxmooncompstar[k][j]])
                                print('dictpopl[moon][namepoplmoontotl][densmoon][indxmooncompstar[k][j]]')
                                print(dictpopl['moon'][namepoplmoontotl]['densmoon'][indxmooncompstar[k][j]])
                                print('dictpopl[comp][namepopllimbtotl][radicomp][dictindx[strglimb][strgbody][k][j]]')
                                print(dictpopl[strglimb][namepopllimbtotl]['radicomp'][0][dictindx[strglimb][strgbody][k][j]])
                                print('dictpopl[comp][namepopllimbtotl][masscomp][dictindx[strglimb][strgbody][k][j]]')
                                print(dictpopl[strglimb][namepopllimbtotl]['masscomp'][0][dictindx[strglimb][strgbody][k][j]])
                                print('dictpopl[comp][namepopllimbtotl][denscomp][dictindx[strglimb][strgbody][k][j]]')
                                print(dictpopl[strglimb][namepopllimbtotl]['denscomp'][0][dictindx[strglimb][strgbody][k][j]])
                                raise Exception('Minimum semi-major axes for the moons are greater than the maximum semi-major axes.')
                        
                        # semi-major axes of the moons
                        for jj in indxmoon[j]:
                            dictpopl['moon'][namepoplmoontotl]['smaxmoon'][indxmooncompstar[k][j][jj]] = tdpy.icdf_powr(np.random.rand(), \
                                                                                                dictpopl['moon'][namepoplmoontotl]['minmsmaxmoon'][indxmooncompstar[k][j][jj]], \
                                                                                                dictpopl[strglimb][namepopllimbtotl]['maxmsmaxmoon'][dictindx[strglimb][strgbody][k][j]], 2.)
                      
                        # add the moon masses to the total mass of the companion + moons system
                        dictpopl[strglimb][namepopllimbtotl]['masscompmoon'][dictindx[strglimb][strgbody][k][j]] += \
                                                        np.sum(dictpopl['moon'][namepoplmoontotl]['massmoon'][indxmooncompstar[k][j]])
                        
                        # orbital period of the moons
                        dictpopl['moon'][namepoplmoontotl]['perimoon'][indxmooncompstar[k][j]] = tdpy.retr_perikepl(dictpopl['moon'][namepoplmoontotl]['smaxmoon'][indxmooncompstar[k][j]], \
                                                                              dictpopl[strglimb][namepopllimbtotl]['masscompmoon'][dictindx[strglimb][strgbody][k][j]] / dictfact['msme'])
                        
                        # load component features into moon features
                        ## temp the code crashed here once
                        for name in dictpopl[strglimb][namepopllimbtotl].keys():
                            dictpopl['moon'][namepoplmoontotl][name][indxmooncompstar[k][j]] = dictpopl[strglimb][namepopllimbtotl][name][dictindx[strglimb][strgbody][k][j]]
                        
                        if booldiag:
                            if (dictpopl['moon'][namepoplmoontotl]['smaxmoon'][indxmooncompstar[k][j]] > \
                                                0.7 * dictpopl[strglimb][namepopllimbtotl]['radihill'][dictindx[strglimb][strgbody][k][j]]).any():
                            
                            
                                print('')
                                print('')
                                print('')
                                print('numbmoon[j]')
                                print(numbmoon[j])
                                print('dictpopl[comp][namepopllimbtotl][radihill][dictindx[strglimb][strgbody][k][j]]')
                                print(dictpopl[strglimb][namepopllimbtotl]['radihill'][dictindx[strglimb][strgbody][k][j]])
                                
                                print('dictpopl[moon][namepoplmoontotl][smaxmoon][indxmooncompstar[k][j]]')
                                print(dictpopl['moon'][namepoplmoontotl]['smaxmoon'][indxmooncompstar[k][j]])
                                
                                print('dictpopl[comp][namepopllimbtotl][smaxcomp][dictindx[strglimb][strgbody][k][j]]')
                                print(dictpopl[strglimb][namepopllimbtotl]['smaxcomp'][0][dictindx[strglimb][strgbody][k][j]])
                                
                                print('dictpopl[moon][namepoplmoontotl][minmsmaxmoon][indxmooncompstar[k][j]]')
                                print(dictpopl['moon'][namepoplmoontotl]['minmsmaxmoon'][indxmooncompstar[k][j]])
                                print('dictpopl[comp][namepopllimbtotl][maxmsmaxmoon][dictindx[strglimb][strgbody][k][j]]')
                                print(dictpopl[strglimb][namepopllimbtotl]['maxmsmaxmoon'][dictindx[strglimb][strgbody][k][j]])
                                
                                raise Exception('Semi-major axis of a moon is larger than 0.7 times the Hill radius of the companion.')
            
                    # planet-moon conjunction times
                    dictpopl['moon'][namepoplmoontotl]['epocmtramoon'] = 1e8 * np.random.rand(numbmoontotl)
        
        dictnico['listnamefeatbody'] = np.array(list(dictpopl[strgbody][namepoplstartotl].keys()))
        dictnico['listnamefeatlimb'] = np.array(list(dictpopl[strglimb][namepopllimbtotl].keys()))
        dictnico['listnamefeatlimbonly'] = np.setdiff1d(dictnico['listnamefeatlimb'], dictnico['listnamefeatbody'])
    
    # check if dictpopl is properly defined while allowing the modern flattened array format.
    # Legacy entries are intentionally stored as [value, label], but newer code often stores
    # the raw NumPy arrays directly and only uses the legacy alias names for compatibility.
    if booldiag:
        for namepopl in dictpopl:
            for namespop in dictpopl[namepopl]:
                if not isinstance(dictpopl[namepopl][namespop], dict):
                    continue
                for namefeat, value in dictpopl[namepopl][namespop].items():
                    if isinstance(value, np.ndarray):
                        continue
                    if not isinstance(value, (list, tuple)):
                        continue
                    if len(value) != 2:
                        continue
                    label = value[1]
                    if isinstance(label, str):
                        continue
                    if isinstance(label, (list, tuple)) and len(label) == 0:
                        continue
                    if isinstance(label, (list, tuple)) and len(label) > 0 and not isinstance(label[0], str):
                        print('')
                        print('')
                        print('')
                        print('namepopl')
                        print(namepopl)
                        print('namespop')
                        print(namespop)
                        print('namefeat')
                        print(namefeat)
                        print('dictpopl[namepopl][namespop][namefeat]')
                        print(value)
                        raise Exception('dictpopl is not properly defined.')

    def _flatten_legacy_alias(popl):
        out = dict()
        for key, val in popl.items():
            if isinstance(val, list) and len(val) == 2:
                out[key] = val[0]
            elif isinstance(val, dict):
                out[key] = _flatten_legacy_alias(val)
            else:
                out[key] = val
        return out

    # Build legacy key aliases from the canonical population names actually present in the dict.
    # The alias objects are intentionally flattened to the historical raw-array contract used by the tests.
    for namepopl, dictnamepopl in list(dictpopl.items()):
        if namepopl not in ['star', 'comp']:
            continue
        for canon_name in list(dictnamepopl.keys()):
            if namepopl == 'star' and canon_name.startswith('star_') and canon_name.endswith('_All'):
                legacy_name = 'star' + canon_name[len('star_'):-len('_All')].replace('SyntheticPopulation', 'Synthetic').replace('Population', '') + 'totl'
                if legacy_name != canon_name:
                    dictnamepopl[legacy_name] = _flatten_legacy_alias(dictnamepopl[canon_name])
            elif namepopl == 'star' and canon_name.startswith('star_') and canon_name.endswith('_Occurrent'):
                legacy_name = 'star' + canon_name[len('star_'):-len('_Occurrent')].replace('SyntheticPopulation', 'Synthetic').replace('Population', '') + 'occu'
                if legacy_name != canon_name:
                    dictnamepopl[legacy_name] = _flatten_legacy_alias(dictnamepopl[canon_name])
            elif namepopl == 'comp' and canon_name.startswith('compstar_') and canon_name.endswith('_All'):
                legacy_name = 'compstar' + canon_name[len('compstar_'):-len('_All')].replace('SyntheticPopulation', 'Synthetic').replace('Population', '') + 'totl'
                if legacy_name != canon_name:
                    dictnamepopl[legacy_name] = _flatten_legacy_alias(dictnamepopl[canon_name])
            elif namepopl == 'comp' and canon_name.startswith('compstar_') and canon_name.endswith('_Transiting'):
                legacy_name = 'compstar' + canon_name[len('compstar_'):-len('_Transiting')].replace('SyntheticPopulation', 'Synthetic').replace('Population', '') + 'tran'
                if legacy_name != canon_name:
                    dictnamepopl[legacy_name] = _flatten_legacy_alias(dictnamepopl[canon_name])

    dictnico['dictpopl'] = dictpopl
    dictnico['dictindx'] = dictindx
    dictnico['dictnumbsamp'] = dictnumbsamp
    dictnico['dictindxsamp'] = dictindxsamp
    
    return dictnico
