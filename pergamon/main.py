import os
import sys

from tqdm import tqdm

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import ephesos
import chalcedon
import miletos
import tdpy
from tdpy import summgene 


def init( \
        # type of analysis
        typeanls, \
        
        # dictionary of features of the populations
        dictpopl=None, \

        # plotting
        typefileplot='png', \

        # verbosity level
        typeverb=1, \
        
        # Boolean flag to interpolate the scatter plots
        boolintpscat=False, \

        # Boolean flag to perform PCA decomposition
        boolprca=False, \

        # Boolean flag to plot scatter distributions of individual populations
        boolplotscatindi=False, \

        # Boolean flag to search for correlations
        boolsrchcorr=False, \
        
        # list of dictionaries of labels and colors for compared populations (one for each comparison)
        listdictlablcolrpopl=None, \
        
        # string holding the generic label of the samples
        lablsampgene=None, \

        # List of Booleans for each list of compared populations, indicating whether the populations consist of mutually-exclusive samples
        listboolcompexcl=None, \
        
        # label for the number of samples
        lablnumbsamp=None, \
        
        # list of feature names for which a cumulative histogram will be made
        listnamefeatcumu=None, \

        # list of titles for the comparison plots
        listtitlcomp=None, \
        
        # list of labels for the samples in the populations
        listlablsamp=None, \

        # type of language for labels and legends
        ## 'English'
        ## 'Turkish'
        typelang='English', \

        # Boolean flag to sort the populations according to their sizes
        boolsortpoplsize=True, \

        # name of the feature holding the labels of the samples
        namefeatlablsamp=None, \
        
        # base path
        pathbase=None, \

        # data path
        pathdata=None, \

        # Boolean flag to diagnose the code using potentially computationally-expensive sanity checks, which may slow down the execution
        booldiag=True, \

        # image path
        pathimag=None, \

        ):
    
    '''
    visualize the populations
    search for correlations in the populations
    search for clusters
    model the density and selection effects to determine the occurence rate
    '''
    
    # construct global object
    gdat = tdpy.gdatstrt()
    
    # copy unnamed inputs to the global object
    for attr, valu in locals().items():
        if '__' not in attr and attr != 'gdat':
            setattr(gdat, attr, valu)

    if gdat.listboolcompexcl is None:
        gdat.listboolcompexcl = []

    if gdat.listtitlcomp is None:
        gdat.listtitlcomp = []
    
    # paths
    if not (gdat.pathbase is not None and gdat.pathimag is None and gdat.pathdata is None or \
            gdat.pathbase is None and gdat.pathimag is not None and gdat.pathdata is not None or \
            gdat.pathbase is None and gdat.pathimag is None and gdat.pathdata is None \
           ):
        print('gdat.pathbase')
        print(gdat.pathbase)
        print('gdat.pathdata')
        print(gdat.pathdata)
        print('gdat.pathimag')
        print(gdat.pathimag)
        raise Exception('')
    
    print('pergamon initialized...')
    
    print('gdat.typeanls')
    print(gdat.typeanls)
        
    if gdat.pathbase is None:
        gdat.pathbase = os.environ['PERGAMON_DATA_PATH'] + '/'
    gdat.pathbase += '%s/' % gdat.typeanls

    if gdat.pathdata is None:
        gdat.pathdata = gdat.pathbase + 'data/'
    
    if gdat.pathimag is None:
        gdat.pathimag = gdat.pathbase + 'imag/'
    os.system('mkdir -p %s' % gdat.pathdata)
    os.system('mkdir -p %s' % gdat.pathimag)
    
    # settings
    ## plotting
    gdat.numbcyclcolrplot = 300
    gdat.alphraww = 0.2
    gdat.strgplotextn = 'pdf'
    ### percentile for zoom plots of relative flux
    gdat.pctlrflx = 95.
    gdat.typefileplot = 'pdf'
    gdat.figrsize = [6, 4]
    gdat.figrsizeydob = [8., 4.]
    gdat.figrsizeydobskin = [8., 2.5]
    
    if gdat.dictpopl is None:
        gdat.dictpopl = dict()
        booldictinpt = False
        
        if gdat.listdictlablcolrpopl is not None:
            raise Exception('')
    else:
        booldictinpt = True
        
        if gdat.listdictlablcolrpopl is None:
            gdat.listnamepopl = list(gdat.dictpopl.keys())
            gdat.numbpopl = len(gdat.listnamepopl)
            gdat.listdictlablcolrpopl = [[] for k in range(gdat.numbpopl)]
            for k in range(gdat.numbpopl):
                gdat.listdictlablcolrpopl[k] = dict()
                gdat.listdictlablcolrpopl[k][gdat.listnamepopl[k]] = [gdat.listnamepopl[k], 'k']
                gdat.listtitlcomp.append('')
                gdat.listboolcompexcl.append(False)
            gdat.typeanls = 'defa'
            gdat.lablsampgene = 'item'
        if gdat.booldiag:
            if gdat.listboolcompexcl is None:
                raise Exception('')
            if gdat.listtitlcomp is None:
                raise Exception('')

        if np.unique(np.array([len(gdat.listdictlablcolrpopl), len(gdat.listboolcompexcl), len(gdat.listtitlcomp)])).size != 1:
            print('')
            print('')
            print('')
            print('gdat.listdictlablcolrpopl')
            print(gdat.listdictlablcolrpopl)
            print('gdat.listboolcompexcl')
            print(gdat.listboolcompexcl)
            print('gdat.listtitlcomp')
            print(gdat.listtitlcomp)
            raise Exception('')

    print('booldictinpt')
    print(booldictinpt)

    if gdat.typeanls == 'micc':
        # merger-induced core collapse supernovae
        ## X-ray (2-10 keV) luminosity of J121001+495647 (VT 1210+4956)
        lumidete = 4e46 # [erg/s]
        distdete = 100. # [Mpc]
        # Mpc / AU
        factmpcsastu = 2.1e11
        
        # TESS magnitude
        magttess = -26.7 - 2.5 * np.log10(lumidete / (100. * factmpcsastu)**2 / 4e34)
        print('magttess')
        print(magttess)

    
    if gdat.typeanls == 'obsvjwstexop':
        path = gdat.pathdata + 'obsvjwstexop.csv'
        print('Reading from %s...' % path)
        
        if not booldictinpt:
            gdat.dictpopl['totl'] = pd.read_csv(path, skiprows=7).to_dict(orient='list')
            for namefeat in gdat.dictpopl['totl'].keys():
                gdat.dictpopl['totl'][namefeat] = np.array(gdat.dictpopl['totl'][namefeat])


    if gdat.typeanls == 'supntess':
        for a in [1, 2, 3]:
            strgcycl = 'cyc%d' % a
            if not booldictinpt:
                path = gdat.pathdata + 'Cycle%d-matched.csv' % a
                print('Reading from %s...' % path)
                gdat.dictpopl[strgcycl] = pd.read_csv(path).to_dict(orient='list')
            
                # change the name of some features
                gdat.dictpopl[strgcycl]['magtdisc'] = gdat.dictpopl[strgcycl].pop('Discovery Mag/Flux')
                gdat.dictpopl[strgcycl]['name'] = gdat.dictpopl[strgcycl].pop('Name')
                
                # only keep floats and integers and remove special characters from feature names
                listname = list(gdat.dictpopl[strgcycl].keys())
                for name in listname:
                    gdat.dictpopl[strgcycl][name] = np.array(gdat.dictpopl[strgcycl][name])
                    if not (gdat.dictpopl[strgcycl][name].dtype == np.float64 or gdat.dictpopl[strgcycl][name].dtype == np.int64 or name == 'name'):
                        del gdat.dictpopl[strgcycl][name]
                    else:
                        gdat.dictpopl[strgcycl][''.join(c for c in name if c not in ' ./:' )] = gdat.dictpopl[strgcycl].pop(name)
                
    gdat.dictindxsamp = dict()
    gdat.dictindxsamp['totl'] = dict()
    gdat.dictnumbsamp = dict()
        
    if gdat.typeanls == 'cosc' or gdat.typeanls == 'psys' or gdat.typeanls == 'plan':
            
        gdat.listcnfgpopl = [ \
                             #['gdr3m140', 's2nr', 'lsstwfdsnomi10yr'], \
                             #['gdr3m190', 's2nr', 'lsstwfdsnomi10yr'], \
                             #['gdr3m140', 's2nr', 'lsstwfdsroll10yr'], \
                             
                             ['ttarsc012min', 's2nr', 'tess2min'], \
                             
                             #['ttare1ms2min', 's2nr', 'tess2min'], \
                             #['ttare1ms20sc', 's2nr', 'tess20sc'], \
                             #['ttare1msffimm140', 's2nr', 'tess10mn'], \
                             #['ttare1msffimm100', 's2nr', 'tess10mn'], \
                             #
                             #['ttare2ms2min', 's2nr', 'tess2min'], \
                             #['ttare2ms20sc', 's2nr', 'tess20sc'], \
                             #['ttare2msffimm140', 's2nr', 'tess200s'], \
                             #
                             #['ttarprmsffim', 'inre', None], \
                             #['ttare1msffim', 'inre', None], \
                             #['ttare2msffimm140', 'inre', None], \
                            ]
        
        gdat.numbiterpopl = len(gdat.listcnfgpopl)
        gdat.indxiterpopl = np.arange(gdat.numbiterpopl)
            
    # dictionary for translation
    if gdat.typelang == 'Turkish':
        gdat.dictturk = tdpy.retr_dictturk()

    if not booldictinpt:
        if gdat.typeanls.startswith('toii'):
            gdat.lablsampgene = 'TOI'
        elif gdat.typeanls == 'autovett':
            gdat.lablsampgene = 'TCE'
        elif gdat.typeanls.startswith('exar'):
            gdat.lablsampgene = 'exoplanet'
        elif gdat.typeanls == 'plan':
            gdat.lablsampgene = 'planet'
        elif gdat.typeanls == 'psys':
            gdat.lablsampgene = 'planetary system'
        elif gdat.typeanls == 'isob':
            gdat.lablsampgene = 'SB'
        elif gdat.typeanls == 'cosc':
            gdat.lablsampgene = 'COSC'
        else:
            print('gdat.typeanls')
            print(gdat.typeanls)
            raise Exception('')
        
        if gdat.lablnumbsamp is None:
            gdat.lablnumbsamp = 'Number of %ss' % gdat.lablsampgene
        
        if gdat.typelang == 'Turkish':
            gdat.lablnumbsamp = gdat.dictturk[gdat.lablnumbsamp]
            gdat.lablsampgene = gdat.dictturk[gdat.lablsampgene]
        
        # get population features
        if gdat.typeanls.startswith('exar'):
            # features of confirmed exoplanets
            gdat.dictpopl['totl'] = chalcedon.retr_dictexar()
            gdat.dictpopl['totl']['noistess'] = chalcedon.retr_noistess(gdat.dictpopl['totl']['vmagsyst'])
            
        if gdat.typeanls == 'cosc' or gdat.typeanls == 'psys' or gdat.typeanls == 'plan':
            
            # type of system
            if gdat.typeanls == 'plan':
                typesyst = 'psys'
            else:
                typesyst = gdat.typeanls
            
            # setup LSST
            ## time stamps
            listtimelsst = np.random.rand(1000)
            
            ## number of years into the mission
            numbyear = 10
            numbsamp = int(numbyear * 100)
            
            minmtime = 0.
            maxmtime = 30.
            
            # cadence
            cade = 2. / 60. / 24. # [days]
            
            numbsamp = 1
            indxsamp = np.arange(numbsamp)
            
            print('gdat.listcnfgpopl')
            print(gdat.listcnfgpopl)
            for r in gdat.indxiterpopl:
                # type of target star population
                typepoplsyst = gdat.listcnfgpopl[r][0]
                print('typepoplsyst')
                print(typepoplsyst)
                
                # type of estimation
                ## 's2nr': based on signal-to-noise estimation
                typeesti = gdat.listcnfgpopl[r][1]
                print('typeesti')
                print(typeesti)

                # type of instrument, cadence, and temporal baseline
                typeinst = gdat.listcnfgpopl[r][2]
                print('typeinst')
                print(typeinst)
                
                # name of star population
                strgpoplstar = 'star' + typepoplsyst

                # name of companion population
                namepoplcomptotl = 'compstar' + typepoplsyst + 'totl'
                namepoplcomptran = 'compstar' + typepoplsyst + 'tran'

                # get a dictionary with features of stars and their companions
                gdat.dictpoplstar, gdat.dictpopl, gdat.dictpoplmoon, gdat.dictnumbsamp, gdat.dictindxsamp, indxcompstar, indxmooncompstar = \
                                                                                                                    ephesos.retr_dictpoplstarcomp(typesyst, typepoplsyst)
                
                # calculate photometric precision for the star population
                if typeinst.startswith('tess'):
                    gdat.dictpopl[namepoplcomptran]['nois'] = chalcedon.retr_noistess(gdat.dictpopl[namepoplcomptran]['tmag'])
                elif typeinst.startswith('lsst'):
                    gdat.dictpopl[namepoplcomptran]['nois'] = chalcedon.retr_noislsst(gdat.dictpopl[namepoplcomptran]['rmag'])
            
                # expected BLS signal detection efficiency
                if typeinst.startswith('lsst'):
                    numbvisi = 1000
                    gdat.dictpopl[namepoplcomptran]['sdee'] = gdat.dictpopl[namepoplcomptran]['depttrancomp'] / 5. / gdat.dictpopl[namepoplcomptran]['nois'] * \
                                                                                                         np.sqrt(gdat.dictpopl[namepoplcomptran]['dcyc'] * numbvisi)
                if typeinst.startswith('tess'):
                    if gdat.typeanls == 'plan':
                        print('namepoplcomptran')
                        print(namepoplcomptran)
                        gdat.dictpopl[namepoplcomptran]['sdee'] = np.sqrt(gdat.dictpopl[namepoplcomptran]['duratrantotl']) * \
                                                                            gdat.dictpopl[namepoplcomptran]['depttrancomp'] / gdat.dictpopl[namepoplcomptran]['nois']
                    if gdat.typeanls == 'cosc':
                        gdat.dictpopl[namepoplcomptran]['sdee'] = np.sqrt(gdat.dictpopl[namepoplcomptran]['duratrantotl']) * \
                                                                                                gdat.dictpopl[namepoplcomptran]['amplslen'] \
                                                                                                                       / gdat.dictpopl[namepoplcomptran]['nois']
                
                # expected detections
                #gdat.dictpopl[namepoplcomptran]['probdeteusam'] = np.exp(-(0.01 * gdat.dictpopl[namepoplcomptran]['pericomp'] * gdat.dictpopl[namepoplcomptran]['numbtsec']))
                #booldeteusam = np.random.rand(gdat.dictpopl[namepoplcomptran]['pericomp'].size) < gdat.dictpopl[namepoplcomptran]['probdeteusam']
                
                #indx = (gdat.dictpopl[namepoplcomptran]['sdee'] > 5) & booldeteusam
                indx = (gdat.dictpopl[namepoplcomptran]['sdee'] > 5)
                chalcedon.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, namepoplcomptran, 'compstar' + typepoplsyst + 'tranposi', indx)

                # expected non-detections
                #indx = (gdat.dictpopl[namepoplcomptran]['sdee'] < 5) | (~booldeteusam)
                indx = (gdat.dictpopl[namepoplcomptran]['sdee'] < 5)
                chalcedon.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, namepoplcomptran, 'compstar' + typepoplsyst + 'trannega', indx)

        if gdat.typeanls.startswith('toii'):
            # features of TOIs
            gdat.dictpopl['totl'] = chalcedon.retr_dicttoii()

        if gdat.typeanls.startswith('hosttoii'):
            # features of hosts of TOIs
            gdat.dictpopl['totl'] = chalcedon.retr_dicthostplan('toii')
        
        if gdat.typeanls.startswith('hostexar'):
            # features of hosts of exoplanets on NASA Exoplanet Archive
            gdat.dictpopl['totl'] = chalcedon.retr_dicthostplan('exar')
    else:
        if gdat.lablnumbsamp is None and gdat.lablsampgene is None:
            print('')
            print('')
            print('gdat.lablnumbsamp')
            print(gdat.lablnumbsamp)
            print('gdat.lablsampgene')
            print(gdat.lablsampgene)
            raise Exception('')

    # subpopulations
    if gdat.typeanls.startswith('exar'):
        
        # indices of samples of subpopulations
        dictindxexar = dict()
        
        # create subpopulations
        ## transiting planets
        dictindxexar['tran'] = np.where(gdat.dictpopl['totl']['booltran'])[0]
            
        ## with transit detections
        dictindxexar['detetran'] = np.where(gdat.dictpopl['totl']['methdisc'] == 'Transit')[0]
        
        ## with transit detections
        dictindxexar['hostsunl'] = np.where((gdat.dictpopl['totl']['tmptstar'] > 5500.) & (gdat.dictpopl['totl']['tmptstar'] < 6000.))[0]
        
        # planets with habitable irradiation
        dictindxexar['irrahabi'] = np.where((gdat.dictpopl['totl']['irra'] > 0.7) & (gdat.dictpopl['totl']['irra'] < 1.1))[0]
        
        # planets with good measured radii and masses
        for strg in ['radi', 'mass']:
            
            indxexargood = np.where( \
                                    (gdat.dictpopl['totl']['%splan' % strg] / gdat.dictpopl['totl']['stdv%splan' % strg] > 5.) & \
                                    (gdat.dictpopl['totl']['stdv%splan' % strg] > 0) \
                                    )[0]
            
            # temp -- this part may not be necessary
            #strgtemp = 'strgrefr%splan' % strg
            #indxexarmeas = []
            #for n  in range(gdat.dictpopl['totl'][strgtemp].size):
            #    if not 'Calculated Value' in gdat.dictpopl['totl'][strgtemp][n]:
            #        indxexarmeas.append(n)
            #indxexarmeas = np.array(indxexarmeas)
            #dictindxexar['prec%s' % strg] = np.setdiff1d(indxexarmeas, indxexargood)
            
            dictindxexar['prec%s' % strg] = indxexargood
            
        # young
        dictindxexar['yong'] = np.where(gdat.dictpopl['totl']['tagestar'] < 0.5)[0]

        # terrestrial mass
        dictindxexar['massterr'] = np.where(gdat.dictpopl['totl']['massplan'] < 10.)[0]

        # less massive than Neptune
        dictindxexar['masslnep'] = np.where((gdat.dictpopl['totl']['massplan'] > 10.) & (gdat.dictpopl['totl']['massplan'] < 18.))[0]

        # less massive than Jupiter
        dictindxexar['massljup'] = np.where((gdat.dictpopl['totl']['massplan'] > 18.) & (gdat.dictpopl['totl']['massplan'] < 300.))[0]

        # more massive than Jupiter
        dictindxexar['massmjup'] = np.where((gdat.dictpopl['totl']['massplan'] > 300.) & (gdat.dictpopl['totl']['massplan'] < 14 * 300.))[0]

        # more massive than Jupiter
        dictindxexar['massgianloww'] = np.where((gdat.dictpopl['totl']['radiplan'] > 10.) & \
                                                    (gdat.dictpopl['totl']['massplan'] > 300. * 0.4) & (gdat.dictpopl['totl']['massplan'] < 300. * 1.))[0]
        dictindxexar['massgianmedi'] = np.where((gdat.dictpopl['totl']['radiplan'] > 10.) & \
                                                    (gdat.dictpopl['totl']['massplan'] > 300. * 1.) & (gdat.dictpopl['totl']['massplan'] < 300. * 2.))[0]
        dictindxexar['massgianhigh'] = np.where((gdat.dictpopl['totl']['radiplan'] > 10.) & \
                                                    (gdat.dictpopl['totl']['massplan'] > 300. * 2.) & (gdat.dictpopl['totl']['massplan'] < 300. * 4.))[0]
        dictindxexar['massgianvhig'] = np.where((gdat.dictpopl['totl']['radiplan'] > 10.) & \
                                                    (gdat.dictpopl['totl']['massplan'] > 300. * 4.) & (gdat.dictpopl['totl']['massplan'] < 300 * 13.))[0]

        # low irradiation
        dictindxexar['irraloww'] = np.where(gdat.dictpopl['totl']['irra'] < 3.)[0]

        # medium irradiation
        dictindxexar['irramedi'] = np.where((gdat.dictpopl['totl']['irra'] > 3.) & (gdat.dictpopl['totl']['irra'] < 30.))[0]

        # high irradiation
        dictindxexar['irrahigh'] = np.where((gdat.dictpopl['totl']['irra'] > 30.) & (gdat.dictpopl['totl']['irra'] < 300.))[0]

        # very high irradiation
        dictindxexar['irravrhi'] = np.where((gdat.dictpopl['totl']['irra'] > 300.) & (gdat.dictpopl['totl']['irra'] < 3000.))[0]
        
        # extreme irradiation
        dictindxexar['irraexhi'] = np.where(gdat.dictpopl['totl']['irra'] > 3000.)[0]
        
        ## with RV detection
        dictindxexar['deteradv'] = np.where(gdat.dictpopl['totl']['methdisc'] == 'Radial Velocity')[0]
        
        ## with other detections (other than transit and RVs)
        dictindxexar['deteothr'] = np.where((gdat.dictpopl['totl']['methdisc'] != 'Radial Velocity') & \
                                                                            (gdat.dictpopl['totl']['methdisc'] != 'Transit'))[0]
        
        ## with imaging detection
        dictindxexar['deteimag'] = np.where(gdat.dictpopl['totl']['methdisc'] == 'Imaging')[0]
        
        ## with microlensing detections
        dictindxexar['detemicr'] = np.where(gdat.dictpopl['totl']['methdisc'] == 'Microlensing')[0]
        
        ## with transit timing variations
        dictindxexar['detetimetran'] = np.where(gdat.dictpopl['totl']['methdisc'] == 'Transit Timing Variations')[0]
        
        ## with eclipse timing variations
        dictindxexar['detetimeeclp'] = np.where(gdat.dictpopl['totl']['methdisc'] == 'Eclipse Timing Variations')[0]
        
        ## with phase variations
        dictindxexar['detephas'] = np.where(gdat.dictpopl['totl']['methdisc'] == 'Orbital Brightness Modulation')[0]
        
        ## with astrometry
        dictindxexar['deteastr'] = np.where(gdat.dictpopl['totl']['methdisc'] == 'Astrometry')[0]
        
        ## discovered by Kepler Telescope
        dictindxexar['kepl'] = np.where(gdat.dictpopl['totl']['facidisc'] == 'Kepler')[0]
        
        ## Kepler discoveries with precise masses 
        dictindxexar['precmasskepl'] = np.intersect1d(dictindxexar['precmass'], dictindxexar['kepl'])

        ## precise mass and radius (precise density)
        dictindxexar['precdens'] = np.intersect1d(dictindxexar['precmass'], dictindxexar['precradi'])

        
        ## can be observed from TUG
        #path = gdat.pathdata + 'trantugg_20220607.csv'
        #if os.path.exists(path):
        #    for p in range(len(indx)):
        #        pass
        #else:
        #    pass
        #indx = np.where(gdat.dictpopl['tran'][''] < 30.)[0]
        #chalcedon.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'tugg', 'tugg', indx)
        #chalcedon.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'tran', 'trantugg', indx)
        
        ## with weak mass
        sigm = gdat.dictpopl['totl']['massplan'] / gdat.dictpopl['totl']['stdvmassplan']
        dictindxexar['weakmass'] = np.where((sigm < 5.) & (sigm > 0.))[0]
        
        ## discovered by TESS
        dictindxexar['tess'] = np.where(gdat.dictpopl['totl']['facidisc'] == 'Transiting Exoplanet Survey Satellite (TESS)')[0]
        
        ## good atmospheric characterization potential
        ### high TSM
        dictindxexar['atmotsmm'] = np.where(gdat.dictpopl['totl']['tsmm'] > 50.)[0]
        ### high ESM
        dictindxexar['atmoesmm'] = np.where(gdat.dictpopl['totl']['esmm'] > 5.)[0]

    if gdat.typeanls.startswith('autovett'):
        
        pathbase = os.environ['DATA'] + '/external/AutomatedVetting/Metrics_Run1/'
        
        listname = ['inv', 'obs', 'inj']
        
        for name in listname:
            path = pathbase + '%s_metrics.csv' % name
            print('Reading from %s...' % path)
            gdat.dictpopl[name] = pd.read_csv(path).to_dict(orient='list')
            for namefeat in gdat.dictpopl[name].keys():
                gdat.dictpopl[name][namefeat] = np.array(gdat.dictpopl[name][namefeat])[:1000]
            #[:1000]
            

    if gdat.typeanls.startswith('toii'):
        # subpopulations
        # total number of TOIs
        numbtoii = len(gdat.dictpopl['totl']['strgcomm'])
        indxtoii = np.arange(numbtoii)

        # indices of samples of subpopulations
        dictindxtoii = dict()
        
        ## PC
        dictindxtoii['pcan'] = np.where(gdat.dictpopl['totl']['typedisptess'] == 'PC')[0]
                
        ## CP
        dictindxtoii['knwn'] = np.where(gdat.dictpopl['totl']['typedisptess'] == 'KP')[0]
        
        ## CP
        dictindxtoii['conp'] = np.where(gdat.dictpopl['totl']['typedisptess'] == 'CP')[0]
        
        ## FP
        dictindxtoii['fpos'] = np.where(gdat.dictpopl['totl']['typedisptess'] == 'EB')[0]
        
        ## sub-Neptunes (radii between two and four times the Earth)
        dictindxtoii['snep'] = np.where((gdat.dictpopl['totl']['radiplan'] > 2.) & (gdat.dictpopl['totl']['radiplan'] < 4.))[0]
        
        ## Faint-star search
        dictindxtoii['fstr'] = []
        for n in indxtoii:
            if isinstance(gdat.dictpopl['totl']['strgcomm'][n], str) and 'found in faint-star QLP search' in gdat.dictpopl['totl']['strgcomm'][n]:
                dictindxtoii['fstr'].append(n)
        dictindxtoii['fstr'] = np.array(dictindxtoii['fstr'])
        
        ## Faint-star search during PM
        dictindxtoii['fstrprms'] = np.intersect1d(dictindxtoii['fstr'], np.where(gdat.dictpopl['totl']['toii'] < 4500)[0])
        
        ## Faint-star search during EM1
        dictindxtoii['fstre1ms'] = np.setdiff1d(dictindxtoii['fstr'], dictindxtoii['fstrprms'])
    
        ## other TOIs
        dictindxtoii['othr'] = np.setdiff1d(indxtoii, dictindxtoii['fstr'])
        
        ## close to the ecliptic
        dictindxtoii['eclp'] = np.where(abs(gdat.dictpopl['totl']['laecstar']) < 10.)[0]
        
        ## distant
        dictindxtoii['dist'] = np.where(gdat.dictpopl['totl']['distsyst'] > 300.)[0]
        
        ## super-Neptunes
        dictindxtoii['larg'] = np.where(gdat.dictpopl['totl']['radiplan'] > 4.)[0]
        
        ## large predicted RV semi-amplitude
        dictindxtoii['rvelsm01'] = np.where(gdat.dictpopl['totl']['rvelsemapred'] > 1.)[0]
        
    
    if gdat.typeanls.startswith('exar') or gdat.typeanls.startswith('toii'):
        
        if gdat.typeanls.startswith('exar'):
            dictindxtemp = dictindxexar
        if gdat.typeanls.startswith('toii'):
            dictindxtemp = dictindxtoii
        
        ## host brighter than the given TESS magnitude
        dictindxtemp['brgttm11'] = np.where(gdat.dictpopl['totl']['tmagsyst'] < 11.)[0]
            
        ## sub-Neptunes (radii between two and four times the Earth)
        dictindxtemp['smal'] = np.where((gdat.dictpopl['totl']['radiplan'] < 4.) & (gdat.dictpopl['totl']['radiplan'] > 2.))[0]
        
        ## multiple transiting planets
        dictindxtemp['tranmult'] = np.where(gdat.dictpopl['totl']['numbplantranstar'] > 1)[0]
    
    ## intersections of subpopulations
    if gdat.typeanls.startswith('exar'):
        dictindxexar['trantess'] = np.intersect1d(dictindxexar['tran'], dictindxexar['tess'])
        dictindxexar['trantessbrgttm11'] = np.intersect1d(dictindxexar['trantess'], dictindxexar['brgttm11'])
        dictindxexar['trantessbrgttm11smal'] = np.intersect1d(dictindxexar['trantessbrgttm11'], dictindxexar['smal'])
        dictindxexar['trankepl'] = np.intersect1d(dictindxexar['tran'], dictindxexar['kepl'])
        dictindxexar['atmo'] = np.union1d(dictindxexar['atmotsmm'], dictindxexar['atmoesmm'])
        dictindxexar['atmoprecmass'] = np.intersect1d(dictindxexar['precmass'], dictindxexar['atmo'])
        dictindxexar['atmoweakmass'] = np.intersect1d(dictindxexar['weakmass'], dictindxexar['atmo'])
        dictindxexar['precmasstess'] = np.intersect1d(dictindxexar['precmass'], dictindxexar['tess'])
    if gdat.typeanls.startswith('toii'):
        for namefrst in ['eclp', 'dist', 'larg', 'pcan', 'knwn', 'conp', 'fpos']:
            for nameseco in ['fstr', 'fstrprms', 'fstre1ms', 'othr']:
                namepoplsubb = namefrst + nameseco
                indx = np.intersect1d(dictindxtoii[nameseco], dictindxtoii[namefrst])
                chalcedon.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'totl', namepoplsubb, indx)
            
    if gdat.typeanls.startswith('toii'):
        dictindxtoii['brgttm11rvelsm01'] = np.intersect1d(dictindxtemp['brgttm11'], dictindxtemp['rvelsm01'])
        namesubpthis = 'fstrbrgttm11rvelsm01'
        dictindxtoii[namesubpthis] = np.intersect1d(dictindxtemp['fstr'], np.intersect1d(dictindxtemp['brgttm11'], dictindxtemp['rvelsm01']))
    
        # Las Campanas Observatory (LCO), Chile
        latiobvt = -29.01418
        longobvt = -70.69239
        heigobvt = 2515.819
        offstimeobvt = 0. # give the times in UT

        strgtimeobvtyear = '2023-01-15 00:00:00'
        listdelttimeobvtyear = np.arange(0., 365.25 / 2., 1. / 24.)
        gdat.dictpopl['totl']['fractimeobsv'] = np.zeros_like(gdat.dictpopl['totl']['rascstar']) - 1.
        for nn in tqdm(range(len(gdat.dictpopl['totl']['declstar']))):
            if nn in dictindxtoii[namesubpthis]:
                dictmile = miletos.main.init( \
                          # provide the RA for the target
                          rasctarg=gdat.dictpopl['totl']['rascstar'][nn], \
                          # provide the Dec for the target
                          decltarg=gdat.dictpopl['totl']['declstar'][nn], \
                          # local time offset with respect to UTC
                          offstimeobvt=offstimeobvt, \
                          # latitude of the observatory
                          latiobvt=latiobvt, \
                          # longtiude of the observatory
                          longobvt=longobvt, \
                          # altitude of the observatory
                          heigobvt=heigobvt, \
                          # time samples
                          listdelttimeobvtyear=listdelttimeobvtyear, \
                          # a string indicating the midnight in the beginning of the observation year
                          strgtimeobvtyear=strgtimeobvtyear, \
                          # turn off time-domain data processing
                          booltserdata=False, \
                          # turn on visibity estimation
                          boolcalcvisi=True, \
                          # turn on visibity plotting
                          boolplotvisi=False, \
                         )
                gdat.dictpopl['totl']['fractimeobsv'][nn] = float(np.where((dictmile['massairr'] < 1.5) & (dictmile['massairr'] > 0.))[0].size) \
                                                                                                                                / dictmile['massairr'].size
        dictindxtemp[namesubpthis + 'lcoo'] = np.where(gdat.dictpopl['totl']['fractimeobsv'] > 0.1)[0]

    if gdat.typeanls.startswith('exar') or gdat.typeanls.startswith('toii'):
        for strg in dictindxtemp.keys():
            chalcedon.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'totl', strg, dictindxtemp[strg])
            
    if gdat.typeanls.startswith('toii') or gdat.typeanls.startswith('hosttoii'):
        
        ## candidate planets in multi systems
        indx = np.where(gdat.dictpopl['pcan']['numbplantranstar'] > 1)[0]
        chalcedon.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'pcan', 'pcanmult', indx)
        
        ## confirmed planets in multi systems
        indx = np.where(gdat.dictpopl['conp']['numbplantranstar'] > 1)[0]
        chalcedon.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'conp', 'conpmult', indx)
        
        ## singles
        #indx = np.where(gdat.dictpopl['pcan']['numbplantranstar'] == 1)[0]
        #chalcedon.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'totl', 'pcansing', indx)
    
    if not booldictinpt:
        if gdat.typeanls.startswith('toii'):
            
            # replace BJD with TESS-truncated BJD (TBJD)
            gdat.dictpopl['totl']['epocmtraplantess'] = gdat.dictpopl['totl']['epocmtraplan'] - 2457000
            
            listkeys = list(gdat.dictpopl['totl'].keys())

            # delete unnecessary features
            for name in listkeys:
                if name.startswith('stdv'):
                    del gdat.dictpopl['totl'][name]
                if name in ['tici', 'epocmtraplan', 'hmagsyst', 'kmagsyst', 'jmagsyst']:
                    del gdat.dictpopl['totl'][name]
        
        if gdat.typeanls == 'qtce':
            
            # base path for the faint-star search
            pathfstr = '/Users/tdaylan/Documents/work/data/external/FaintStars/'
            
            ## TESS sectors
            listtsecprms = range(1, 27)
            #listtsecprms = range(1, 3)
            numbtsecprms = len(listtsecprms)
            indxtsecprms = np.arange(numbtsecprms)
            
            # header
            listnamevarbfstr = 'tic per per_err epo epo_err rprs rprs_err b b_err ars ars_err rprs_odd rprs_err_odd rprs_even rprs_err_even dep dur din u1 '
            listnamevarbfstr += 'u2 SN SPN OOTmag sig_amp sig_pri sig_sec sig_ter sig_pos sig_oe dmm shape asym sig_fa1 sig_fa2 fred phs_pri phs_sec phs_ter phs_pos '
            listnamevarbfstr += 'dep_sec deperr_sec sig_oe_alt sig_12 sig_23 sig_13'
            listnamevarbfstr = listnamevarbfstr.split(' ')
            for n in range(len(listnamevarbfstr)):
                listnamevarbfstr[n] = ''.join(listnamevarbfstr[n].split('_'))
            
            print('listnamevarbfstr')
            print(listnamevarbfstr)
            
            gdat.dictpopl['totl'] = dict()
            for o in indxtsecprms:
                gdat.dictpopl['ts%02d' % o] = dict()
            for namevarbfstr in listnamevarbfstr:
                for o in indxtsecprms:
                    gdat.dictpopl['ts%02d' % o][namevarbfstr] = []
                
            print('Reading metrics...')
            for o in indxtsecprms:
                path = pathfstr + 'metr/metr_sc%02d.txt' % listtsecprms[o]
                print('Reading from %s...' % path)
                objt = open(path, 'r')
                for line in objt:
                    linesplt = line.split(',')

                    numbfiel = len(linesplt)
                    
                    if '' in linesplt:
                        continue
                    
                    peri = float(linesplt[1])
                    peri_err = float(linesplt[2])

                    if peri == -1 or peri_err == -1:# or peri_err > 1e10:
                        continue
                    
                    for n in range(numbfiel):
                        if listnamevarbfstr[n] == 'tic':
                            valu = int(linesplt[n])
                        else:
                            valu = float(linesplt[n])
                        gdat.dictpopl['ts%02d' % o][listnamevarbfstr[n]].append(valu)

                    # check
                    if numbfiel != 45:
                        raise Exception('')
                for namevarbfstr in listnamevarbfstr:
                    gdat.dictpopl['ts%02d' % o][namevarbfstr] = np.array(gdat.dictpopl['ts%02d' % o][namevarbfstr]) 
                    #if gdat.dictpopl['ts%02d' % o][namevarbfstr].size == 0:
                    #    raise Exception('')
                print('%d metrics...' % gdat.dictpopl['ts%02d' % o]['tic'].size)
            print('')
            
            # read TIC IDs and dispositions
            print('Reading dispositions...')
            pathbase = pathfstr + 'combined/'
            ticihvettsec = [[] for o in indxtsecprms]
            disphvettsec = [[] for o in indxtsecprms]
            for o in indxtsecprms:
                path = pathbase + 'sector-%d.tier2.csv' % listtsecprms[o]
                print('Reading from %s...' % path)
                dictquer = pd.read_csv(path).to_dict(orient='list')
                ticihvettsec[o] = np.array(dictquer['TIC'])
                disphvettsec[o] = dictquer['Final']
                print('%d dispositions...' % ticihvettsec[o].size)
            print('')
            
            print('Filtering out those targets for which there is a disposition, but not metric...')
            ticihvettsectemp = [[] for o in indxtsecprms]
            for o in indxtsecprms:
                for ticihvet in ticihvettsec[o]:
                    if ticihvet in gdat.dictpopl['ts%02d' % o]['tic']:
                        ticihvettsectemp[o].append(ticihvet)
                ticihvettsectemp[o] = np.array(ticihvettsectemp[o])
            for o in indxtsecprms:
                print('Sector %2d: %4d of %4d dispositions have metrics...' % (listtsecprms[o], ticihvettsectemp[o].size, ticihvettsec[o].size))
            ticihvettsec = ticihvettsectemp

            print('Merging lists of TIC IDs for targets with metrics...')
            ticiastrtsec = [[] for o in indxtsecprms]
            for o in indxtsecprms:
                ticiastrtsec[o] = gdat.dictpopl['ts%02d' % o]['tic']
            ticiastrconc = np.concatenate(ticiastrtsec)
            print('Total number of metrics: %d' % ticiastrconc.size)
            ticiastruniq, indxuniq, indxinve, numbtici = np.unique(ticiastrconc, return_index=True, return_inverse=True, return_counts=True)
            print('Total number of targets: %d' % ticiastruniq.size)

            print('Merging lists of TIC IDs for targets with dispositions...')
            ticihvetconc = np.concatenate(ticihvettsec)
            print('Total number of dispositions: %d' % ticihvetconc.size)
            ticihvetuniq, indxuniq, indxinve, numbtici = np.unique(ticihvetconc, return_index=True, return_inverse=True, return_counts=True)
            print('Total number of targets with dispositions: %d' % ticihvetuniq.size)

            #gdat.dictpopl['astr'][namevarbfstr] = []
            #gdat.dictpopl['hvet'][namevarbfstr] = []
            #gdat.dictpopl['fpos'][namevarbfstr] = []
            #gdat.dictpopl['pcan'][namevarbfstr] = []
            
            for a in range(2):
                if a == 0:
                    ticitsec = ticiastrtsec
                    ticiuniq = ticiastruniq
                if a == 1:
                    ticitsec = ticihvettsec
                    disptsec = disphvettsec
                    ticiuniq = ticihvetuniq
                
                for tici in ticiuniq:
                    # get dispositions and sector indices for this target
                    ## sector indices of this target
                    indxtsecthis = []
                    if a == 1:
                        ## dispositions of this target
                        dispthis = []
                    for o in indxtsecprms:
                        indx = np.where(tici == ticitsec[o])[0]
                        if indx.size > 0:
                            if a == 1:
                                dispthis.append(disptsec[o][indx[0]])
                            indxtsecthis.append(o)
                    
                    # index of the last sector of this target
                    indxtseclast = indxtsecthis[-1]

                    ## find the metric index of the target in the last sector
                    indx = np.where(gdat.dictpopl['ts%02d' % indxtseclast]['tic'] == tici)[0]
                    
                    if indx.size == 0:
                        raise Exception('')

                    ## collect metrics of the target in the last sector
                    if a == 0:
                        ### all human-vetted targets
                        for namevarbfstr in listnamevarbfstr:
                            gdat.dictpopl['totl'][namevarbfstr].append(gdat.dictpopl['astr'][indxtseclast][namevarbfstr][indx[0]])
                    if a == 1:
                        ### all human-vetted targets
                        for namevarbfstr in listnamevarbfstr:
                            gdat.dictpopl['hvet'][namevarbfstr].append(gdat.dictpopl['astr'][indxtseclast][namevarbfstr][indx[0]])
                        ### Ps
                        if dispthis[-1] == 'P':
                            for namevarbfstr in listnamevarbfstr:
                                gdat.dictpopl['pcan'][namevarbfstr].append(gdat.dictpopl['astr'][indxtseclast][namevarbfstr][indx[0]])
                        ### Fs
                        else:
                            for namevarbfstr in listnamevarbfstr:
                                gdat.dictpopl['fpos'][namevarbfstr].append(gdat.dictpopl['astr'][indxtseclast][namevarbfstr][indx[0]])
             
                
            for namevarbfstr in listnamevarbfstr:
                gdat.dictpopl['astr'][namevarbfstr] = np.array(gdat.dictpopl['astr'][namevarbfstr])
                if len(gdat.dictpopl['pcan'][namevarbfstr]) == 0:
                    raise Exception('')
                gdat.dictpopl['fpos'][namevarbfstr] = np.array(gdat.dictpopl['fpos'][namevarbfstr])
                gdat.dictpopl['pcan'][namevarbfstr] = np.array(gdat.dictpopl['rpcan'][namevarbfstr])
                gdat.dictpopl['hvet'][namevarbfstr] = np.array(gdat.dictpopl['hvet'][namevarbfstr])
            
            # crossmatch with TIC to get additional features
            #dictquer = ephesos.xmat_tici(ticifstrastr)

            #numbastr = len(ticifstrastr)
            #for namevarbfstr in listnamevarbfstr:
            #    gdat.dictpopl['astr'][namevarbfstr] = np.empty(numbastr)
            #for k, tici in enumerate(ticifstrastr):
            #    indx = np.where(ticifstr == tici)[0]
            #    if indx.size > 0:
            #        for namevarbfstr in listnamevarbfstr:
            #            gdat.dictpopl['astr'][namevarbfstr][k] = gdat.dictpopl['astr'][namevarbfstr][indx[0]]

            #figr, axis = plt.subplots(figsize=(4, 4))
            #axis.hist(numbtici)
            #axis.set_yscale('log')
            #axis.set_xlabel('Number of sectors')
            #axis.set_ylabel('N')
            #plt.tight_layout()
            #path = gdat.pathimag + 'histnumbtsec.%s' % (gdat.strgplotextn)
            #print('Writing to %s...' % path)
            #plt.savefig(path)
            #plt.close()
    
        if gdat.typeanls == 'hjuppcur':
            listlabl = ['log $g_{\mathrm{p}}$ [cgs]', r'$T_{\mathrm{eq}}$ [K]', '[Fe/H] [dex]', r'$T_{\mathrm{day}}$ [K]', '$A_g$']
            listname = ['logg', 'ptmp', 'meta', 'tday', 'albe']
            path = gdat.pathdata + 'catlpcurtess.csv'
            # get the dictionaries holding the population features
            gdat.dictpopl['hjuppcur'] = pd.read_csv(path)

        if gdat.typeanls == 'isob':
            gdat.dictpopl['totl']['nois'] = ephesos.samp_paraorbtsols()
        
    gdat.listnamepopl = np.array(list(gdat.dictpopl.keys()))
    
    # first item is becomes the y-axis, the second item becomes the x-axis
    gdat.listnameordrpair = [['radiplan', 'tmptplan']]

    # to be deleted?
    #if gdat.typeanls.startswith('exar'):
    #    gdat.lablpoplmast = 'Exoplanets'
    #if gdat.typeanls.startswith('hostexar'):
    #    gdat.lablpoplmast = 'Exoplanet hosts'
    #if gdat.typeanls.startswith('toii'):
    #    gdat.lablpoplmast = 'TOIs'
    #if gdat.typeanls.startswith('hosttoii'):
    #    gdat.lablpoplmast = 'TOI Hosts'
    
    if gdat.namefeatlablsamp is None and not booldictinpt:
        if gdat.typeanls.startswith('hosttoii') or gdat.typeanls.startswith('hostexar'):
            gdat.namefeatlablsamp = 'namestar'
        elif gdat.typeanls.startswith('cosc'):
            gdat.namefeatlablsamp = None
        elif gdat.typeanls.startswith('autovett'):
            gdat.namefeatlablsamp = 'tic'
        elif gdat.typeanls.startswith('toii'):
            gdat.namefeatlablsamp = 'nametoii'
        elif gdat.typeanls.startswith('exar'):
            gdat.namefeatlablsamp = 'nameplan'
        elif gdat.typeanls == 'psys':
            gdat.namefeatlablsamp = None
        elif gdat.typeanls == 'plan':
            gdat.namefeatlablsamp = None
        elif gdat.typeanls == 'supntess':
            gdat.namefeatlablsamp = 'name'
        elif gdat.typeanls == 'defa':
            gdat.namefeatlablsamp = None
        else:
            print('')
            print('')
            print('')
            print('gdat.typeanls')
            print(gdat.typeanls)
            raise Exception('')
        
    # fill the list of dictionaries with key carriying the name of the population and elements carrying the label and color
    if not booldictinpt:
        gdat.listdictlablcolrpopl = []

        if gdat.typeanls == 'targ_cosc':
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listdictlablcolrpopl[-1]['statu0ne'] = ['TS 0', 'blue']
            gdat.listdictlablcolrpopl[-1]['statu1ne'] = ['TS 1', 'g']
            gdat.listboolcompexcl.append(True)
        
        if gdat.typeanls == 'isob':
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listdictlablcolrpopl[-1]['isob'] = ['ISO', 'blue']
            gdat.listdictlablcolrpopl[-1]['ETNO'] = ['ETNO', 'gray']
            gdat.listboolcompexcl.append(True)
        
        if gdat.typeanls == 'supntess':
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listdictlablcolrpopl[-1]['cyc1'] = ['Cycle 1', 'gray']
            gdat.listdictlablcolrpopl[-1]['cyc2'] = ['Cycle 2', 'gray']
            gdat.listdictlablcolrpopl[-1]['cyc3'] = ['Cycle 3', 'gray']
        
        if gdat.typeanls == 'cosc' or gdat.typeanls == 'plan':
            
            print('gdat.listcnfgpopl')
            print(gdat.listcnfgpopl)
            for r in gdat.indxiterpopl:
                # type of target star population
                typepoplsyst = gdat.listcnfgpopl[r][0]
                strgpoplstar = 'star' + typepoplsyst
                
                print('typepoplsyst')
                print(typepoplsyst)

                gdat.listdictlablcolrpopl.append(dict())
                if gdat.typeanls == 'cosc':
                    gdat.listtitlcomp.append('Compact Objects with a stellar companion')
                if gdat.typeanls == 'plan':
                    gdat.listtitlcomp.append('Planets')
                gdat.listdictlablcolrpopl[-1]['comp'+strgpoplstar+'totl'] = ['All', 'black']
                gdat.listdictlablcolrpopl[-1]['comp'+strgpoplstar+'tran'] = ['Transiting', 'blue']
                gdat.listdictlablcolrpopl[-1]['comp'+strgpoplstar+'tranposi'] = ['Detected', 'green']
                gdat.listboolcompexcl.append(False)
                
                gdat.listdictlablcolrpopl.append(dict())
                if gdat.typeanls == 'cosc':
                    gdat.listtitlcomp.append('Compact Objects with a stellar companion')
                if gdat.typeanls == 'plan':
                    gdat.listtitlcomp.append('Planets')
                gdat.listdictlablcolrpopl[-1]['comp'+strgpoplstar+'tran'] = ['Transiting', 'blue']
                gdat.listdictlablcolrpopl[-1]['comp'+strgpoplstar+'tranposi'] = ['Detected', 'green']
                gdat.listboolcompexcl.append(False)
        
                gdat.listdictlablcolrpopl.append(dict())
                if gdat.typeanls == 'cosc':
                    gdat.listtitlcomp.append('Compact Objects with a stellar companion')
                if gdat.typeanls == 'plan':
                    gdat.listtitlcomp.append('Planets')
                gdat.listdictlablcolrpopl[-1]['comp'+strgpoplstar+'trannega'] = ['Not Detected', 'firebrick']
                gdat.listdictlablcolrpopl[-1]['comp'+strgpoplstar+'tranposi'] = ['Detected', 'green']
                gdat.listboolcompexcl.append(True)
        
        if gdat.typeanls.startswith('exar') or gdat.typeanls.startswith('hostexar'):
            
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('Exoplanet Detections')
            gdat.listdictlablcolrpopl[-1]['detetran'] = ['Transit', 'blue']
            gdat.listdictlablcolrpopl[-1]['deteradv'] = ['Radial Velocity', 'firebrick']
            gdat.listdictlablcolrpopl[-1]['deteothr'] = ['Other', 'orange']
            gdat.listboolcompexcl.append(True)
        
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('Exoplanet detections via direct imaging')
            gdat.listdictlablcolrpopl[-1]['deteimag'] = ['Direct imaging', 'deepskyblue']
            gdat.listboolcompexcl.append(True)
        
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('Exoplanet detections')
            gdat.listdictlablcolrpopl[-1]['deteimag'] = ['Imaging', 'gold']
            gdat.listdictlablcolrpopl[-1]['detemicr'] = ['Microlensing', 'darkgreen']
            gdat.listdictlablcolrpopl[-1]['deteastr'] = ['Astrometry', 'olive']
            gdat.listdictlablcolrpopl[-1]['detephas'] = ['Orbital Brightness Modulation', 'purple']
            gdat.listdictlablcolrpopl[-1]['detetimetran'] = ['Transit Timing Variations', 'magenta']
            gdat.listdictlablcolrpopl[-1]['detetimeeclp'] = ['Eclipse Timing Variations', 'deepskyblue']
            gdat.listboolcompexcl.append(True)
        
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('Exoplanets')
            gdat.listdictlablcolrpopl[-1]['hostsunl'] = ['Sun-like host', 'blue']
            gdat.listboolcompexcl.append(True)
        
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('Exoplanets')
            gdat.listdictlablcolrpopl[-1]['totl'] = ['All', 'black']
            gdat.listdictlablcolrpopl[-1]['tran'] = ['Transiting', 'blue']
            gdat.listboolcompexcl.append(False)
        
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('Transiting exoplanets discovered by TESS')
            gdat.listdictlablcolrpopl[-1]['trantess'] = ['All', 'black']
            gdat.listdictlablcolrpopl[-1]['trantessbrgttm11'] = ['Bright', 'purple']
            gdat.listdictlablcolrpopl[-1]['trantessbrgttm11smal'] = ['Bright and small', 'orange']
            gdat.listboolcompexcl.append(False)
        
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('Exoplanets')
            gdat.listdictlablcolrpopl[-1]['totl'] = ['All', 'gray']
            gdat.listdictlablcolrpopl[-1]['atmo'] = ['High TSM or ESM', 'blue']
            gdat.listdictlablcolrpopl[-1]['atmoprecmass'] = ['High TSM or ESM \& Precise mass', 'orange']
            gdat.listdictlablcolrpopl[-1]['atmoweakmass'] = ['High TSM or ESM \& Weak mass', 'firebrick']
            gdat.listboolcompexcl.append(False)
            
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('Exoplanets')
            gdat.listdictlablcolrpopl[-1]['totl'] = ['All', 'gray']
            gdat.listdictlablcolrpopl[-1]['tess'] = ['TESS discoveries', 'blue']
            gdat.listboolcompexcl.append(False)
            
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('Exoplanets')
            gdat.listdictlablcolrpopl[-1]['kepl'] = ['Kepler discoveries', 'firebrick']
            gdat.listdictlablcolrpopl[-1]['tess'] = ['TESS discoveries', 'blue']
            gdat.listboolcompexcl.append(True)
            
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('Exoplanets with precise masses')
            gdat.listdictlablcolrpopl[-1]['precmasskepl'] = ['Kepler discoveries', 'firebrick']
            gdat.listdictlablcolrpopl[-1]['precmasstess'] = ['TESS discoveries', 'blue']
            gdat.listboolcompexcl.append(True)
            
            #gdat.listdictlablcolrpopl.append(dict())
            #gdat.listtitlcomp.append('Exoplanets')
            #gdat.listdictlablcolrpopl[-1]['totl'] = ['All', 'gray']
            #gdat.listdictlablcolrpopl[-1]['yong'] = ['Young', 'blue']
            #gdat.listboolcompexcl.append(False)
            #
            #gdat.listdictlablcolrpopl.append(dict())
            #gdat.listtitlcomp.append('Exoplanets with precise masses')
            #gdat.listdictlablcolrpopl[-1]['atmoprecmass'] = ['Precise mass', 'blue']
            #gdat.listboolcompexcl.append(True)
            #
            #gdat.listdictlablcolrpopl.append(dict())
            #gdat.listtitlcomp.append('Exoplanets with weak mass')
            #gdat.listdictlablcolrpopl[-1]['atmoweakmass'] = ['Weak mass', 'blue']
            #gdat.listboolcompexcl.append(True)
            
            #gdat.listdictlablcolrpopl.append(dict())
            #gdat.listtitlcomp.append('Exoplanets')
            #gdat.listdictlablcolrpopl[-1]['totl'] = ['All', 'gray']
            #gdat.listdictlablcolrpopl[-1]['precdens'] = ['With precise mass and radius', 'blue']
            #gdat.listboolcompexcl.append(False)
            #
            #gdat.listdictlablcolrpopl.append(dict())
            #gdat.listtitlcomp.append('Exoplanets')
            #gdat.listdictlablcolrpopl[-1]['totl'] = ['All', 'gray']
            #gdat.listdictlablcolrpopl[-1]['irrahabi'] = ['Habitable zone', 'green']
            #gdat.listboolcompexcl.append(False)
            # 
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('Exoplanets')
            gdat.listdictlablcolrpopl[-1]['massterr'] = ['Terrestrial-mass', 'orange']
            gdat.listdictlablcolrpopl[-1]['masslnep'] = ['Less massive than Neptune', 'b']
            gdat.listdictlablcolrpopl[-1]['massljup'] = ['Less massive than Jupiter', 'r']
            gdat.listdictlablcolrpopl[-1]['massmjup'] = ['More massive than Jupiter', 'g']
            gdat.listboolcompexcl.append(True)
             
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('Exoplanets')
            gdat.listdictlablcolrpopl[-1]['massgianloww'] = ['$0.4 M_J < M < 1 M_J$', 'orange']
            gdat.listdictlablcolrpopl[-1]['massgianmedi'] = ['$1 M_J   < M < 2 M_J$', 'b']
            gdat.listdictlablcolrpopl[-1]['massgianhigh'] = ['$2 M_J   < M < 4 M_J$', 'r']
            gdat.listdictlablcolrpopl[-1]['massgianvhig'] = ['$4 M_J   < M < 13 M_J$', 'g']
            gdat.listboolcompexcl.append(True)
             
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('Exoplanets')
            #gdat.listdictlablcolrpopl[-1]['irramedi'] = ['Medium irradiation', 'g']
            gdat.listdictlablcolrpopl[-1]['irrahigh'] = ['High irradiation', 'b']
            gdat.listdictlablcolrpopl[-1]['irravrhi'] = ['Very high irradiation', 'r']
            gdat.listdictlablcolrpopl[-1]['irraexhi'] = ['Extreme irradiation', 'g']
            gdat.listboolcompexcl.append(True)
             
        if gdat.typeanls.startswith('autovett'):
            
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('TCEs')
            gdat.listdictlablcolrpopl[-1]['obs'] = ['Observed', 'g']
            gdat.listdictlablcolrpopl[-1]['inj'] = ['Injected', 'r']
            gdat.listdictlablcolrpopl[-1]['inv'] = ['Inverted', 'b']
            gdat.listboolcompexcl.append(True)
             
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('TCEs')
            gdat.listdictlablcolrpopl[-1]['obs'] = ['Observed', 'g']
            gdat.listboolcompexcl.append(True)
             
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('TCEs')
            gdat.listdictlablcolrpopl[-1]['inj'] = ['Injected', 'r']
            gdat.listboolcompexcl.append(True)
             
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('TCEs')
            gdat.listdictlablcolrpopl[-1]['inv'] = ['Inverted', 'b']
            gdat.listboolcompexcl.append(True)
             
            
        if gdat.typeanls.startswith('toii') or gdat.typeanls.startswith('hosttoii'):
            
            gdat.dicttoiistat = dict()
            gdat.dictlablcolrtoiistat = dict()
            gdat.dictlablcolrtoiistat['pcan'] = ['Planet\n Candidate', 'blue']
            gdat.dictlablcolrtoiistat['fpos'] = ['False\n Positive', 'firebrick']
            gdat.dictlablcolrtoiistat['knwn'] = ['Known\n Planet', 'orange']
            gdat.dictlablcolrtoiistat['conp'] = ['Confirmed\n Planet', 'green']
            gdat.listnametoiistat = list(gdat.dictlablcolrtoiistat.keys())
            
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('')
            gdat.listboolcompexcl.append(True)
            gdat.listdictlablcolrpopl[-1]['totl'] = ['', 'k']
            
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('TOIs')
            gdat.listboolcompexcl.append(True)
            gdat.listdictlablcolrpopl[-1]['totl'] = ['All', 'k']
            gdat.listdictlablcolrpopl[-1]['conp'] = ['Confirmed', 'green']
            
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('TOIs')
            gdat.listboolcompexcl.append(False)
            gdat.listdictlablcolrpopl[-1]['totl'] = ['All', 'k']
            gdat.listdictlablcolrpopl[-1]['fstr'] = ['Faint-star', 'green']
            
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('TOIs')
            for name in gdat.listnametoiistat:
                gdat.listdictlablcolrpopl[-1][name] = gdat.dictlablcolrtoiistat[name]
            gdat.listboolcompexcl.append(True)
        
            for a in range(2):
                
                if a == 0:
                    listnametypetoii = ['othr', 'fstr', 'fstrprms', 'fstre1ms']
                    listlabltypetoii = ['Other TOIs', 'Faint-star search', 'Faint-star search during PM', 'Faint-star search during EM1']
                    listcolrtypetoii = ['black', 'green', 'orange', 'blue']
                    gdat.listboolcompexcl.append(False)
                else:
                    listnametypetoii = ['othr', 'fstrprms', 'fstre1ms']
                    listlabltypetoii = ['Other TOIs', 'Faint-star search during PM', 'Faint-star search during EM1']
                    listcolrtypetoii = ['black', 'orange', 'blue']
                    gdat.listboolcompexcl.append(True)
                
                gdat.listdictlablcolrpopl.append(dict())
                gdat.listtitlcomp.append('TOIs')
                for k in range(len(listlabltypetoii)):
                    gdat.listdictlablcolrpopl[-1][listnametypetoii[k]] = [listlabltypetoii[k], listcolrtypetoii[k]]

                ## TOIs close to the ecliptic
                #gdat.listdictlablcolrpopl.append(dict())
                #gdat.listtitlcomp.append('TOIs close to the ecliptic')
                #for k in range(len(listlabltypetoii)):
                #    gdat.listdictlablcolrpopl[-1]['eclp' + listnametypetoii[k]] = [listlabltypetoii[k], listcolrtypetoii[k]]
                #
                ## distant TOIs
                #gdat.listdictlablcolrpopl.append(dict())
                #gdat.listtitlcomp.append('TOIs beyond 300 pc')
                #for k in range(len(listlabltypetoii)):
                #    gdat.listdictlablcolrpopl[-1]['dist' + listnametypetoii[k]] = [listlabltypetoii[k], listcolrtypetoii[k]]
                #
                ## large TOIs
                #gdat.listdictlablcolrpopl.append(dict())
                #gdat.listtitlcomp.append(r'TOIs larger than 4$R_{\oplus}$')
                #for k in range(len(listlabltypetoii)):
                #    gdat.listdictlablcolrpopl[-1]['larg' + listnametypetoii[k]] = [listlabltypetoii[k], listcolrtypetoii[k]]
            
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('Faint-Star TOIs')
            gdat.listdictlablcolrpopl[-1]['pcanfstr'] = gdat.dictlablcolrtoiistat['pcan']
            gdat.listdictlablcolrpopl[-1]['fposfstr'] = gdat.dictlablcolrtoiistat['fpos']
            gdat.listdictlablcolrpopl[-1]['knwnfstr'] = gdat.dictlablcolrtoiistat['knwn']
            gdat.listdictlablcolrpopl[-1]['conpfstr'] = gdat.dictlablcolrtoiistat['conp']
            gdat.listboolcompexcl.append(True)
        
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('Faint-Star TOIs')
            gdat.listdictlablcolrpopl[-1]['fstr'] = ['', 'green']
            gdat.listboolcompexcl.append(True)
        
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('Confirmed Faint-Star TOIs')
            gdat.listdictlablcolrpopl[-1]['conpfstr'] = ['', gdat.dictlablcolrtoiistat['conp'][1]]
            gdat.listboolcompexcl.append(True)
        
            print('gdat.listdictlablcolrpopl')
            print(gdat.listdictlablcolrpopl)
            #gdat.listdictlablcolrpopl.append(dict())
            #gdat.listtitlcomp.append('TOIs')
            #gdat.listdictlablcolrpopl[-1]['totl'] = ['All', 'gray']
            #gdat.listdictlablcolrpopl[-1]['fstrbrgttm11rvelsm01lcoo'] = \
            #                                    ['Found by the faint-star search, Host $T>10$, large predicted RV semi-amplitude, and visible from LCO', 'green']
            #gdat.listboolcompexcl.append(False)
        
            #gdat.listdictlablcolrpopl.append(dict())
            #gdat.listtitlcomp.append('Faint-star search during EM1')
            #gdat.listdictlablcolrpopl[-1]['fstre1ms'] = ['', 'green']
            #gdat.listboolcompexcl.append(True)
        
        if gdat.typeanls == 'hosttoii':
            pass

        if gdat.typeanls == 'hostexar':
            pass

        if gdat.typeanls == 'qtce':
            gdat.listdictlablcolrpopl.append(dict())
            gdat.listtitlcomp.append('TCEs')
            gdat.listdictlablcolrpopl[-1]['totl'] = ['All', 'black']
            gdat.listdictlablcolrpopl[-1]['tran'] = ['Transiting', 'blue']
            gdat.listboolcompexcl.append(False)
    

    if len(gdat.listdictlablcolrpopl) == 0:
        print('')
        print('')
        print('')
        print('List of dictionaries with labels and colors for the populations (listdictlablcolrpopl) should be provided as input.')
        raise Exception('')
    
    gdat.numbpopl = len(gdat.listnamepopl)
    if gdat.numbpopl == 0:
        raise Exception('')

    gdat.indxpopl = np.arange(gdat.numbpopl)
    
    listnamefeat = [[] for k in gdat.indxpopl]
    gdat.indxfeat = [[] for k in gdat.indxpopl]
    
    numbfeat = np.empty(gdat.numbpopl, dtype=int)
    for k in gdat.indxpopl:
        listnamefeat[k] = list(gdat.dictpopl[gdat.listnamepopl[k]].keys())
        
        numbfeat[k] = len(listnamefeat[k])
        if numbfeat[k] == 0:
            print('gdat.listnamepopl[k]')
            print(gdat.listnamepopl[k])
            print('gdat.dictpopl[gdat.listnamepopl[k]]')
            print(gdat.dictpopl[gdat.listnamepopl[k]])
            print('listnamefeat[k]')
            print(listnamefeat[k])
            raise Exception('')
        gdat.indxfeat[k] = np.arange(numbfeat[k])
    
    for k in gdat.indxpopl:
        
        namepopl = gdat.listnamepopl[k]

        # check for finiteness
        #indx = np.where(np.isfinite(gdat.dictpopl['massstar']))[0]
        #for name in listnamefeattemp:
        #    gdat.dictpopl[name] = gdat.dictpopl[name][indx]

        #indx = np.where(np.isfinite(gdat.dictpopl['duratran']))[0]
        #for name in listnamefeattemp:
        #    gdat.dictpopl[name] = gdat.dictpopl[name][indx]
        
        # check number of targets for each feature
        #numbkeys = len(listnamefeat[k])
        #numbtarg = np.empty(numbkeys, dtype=int)
        #for n in gdat.indxfeat[k]:
        #    print('gdat.dictpopl')
        #    print(gdat.dictpopl)
        #    numbtarg[n] = gdat.dictpopl[namepopl][listnamefeat[k][n]].size
        #    if numbtarg[n] == 0:
        #        print('Feature %s of population %s is empty!' % (listnamefeat[k][n], gdat.listnamepopl[k]))
        #if np.unique(numbtarg).size > 1:
        #    print('Number of targets is not the same for every feature.')
        #    print('k')
        #    print(k)
        #    print('namepopl')
        #    print(namepopl)
        #    print('listnamefeat[k]')
        #    print(listnamefeat[k])
        #    for n in gdat.indxfeat[k]:
        #        print('gdat.dictpopl[namepopl][listnamefeat[k][n]]')
        #        summgene(gdat.dictpopl[namepopl][listnamefeat[k][n]])
        #    raise Exception('')
        
    # filter some features of the population
    listsampfilt = [[] for k in gdat.indxpopl]
    listnamefeatfilt = [[] for k in gdat.indxpopl]
    gdat.dictpoplfilt = dict()

    if gdat.listlablsamp is None and gdat.namefeatlablsamp is not None:
        gdat.listlablsamp = [[] for k in gdat.indxpopl]

    for k in gdat.indxpopl:
        gdat.dictpoplfilt[gdat.listnamepopl[k]] = dict()
        
        for n in gdat.indxfeat[k]:
            
            boolgood = False
            if (gdat.typeanls.startswith('toii') or gdat.typeanls.startswith('hosttoii')) and listnamefeat[k][n] in [ \
                                                                           'radistar', 'massstar', 'metastar', 'loggstar', 'tagetar', 'distsyst', 'numbplantranstar', \
                                                                               'tmagsyst', 'vmagsyst']:
                boolgood = True
            
            if gdat.typeanls.startswith('toii') and listnamefeat[k][n] in [ \
                                                                                 'declstar', 'rascstar', 'laecstar', 'loecstar', \
                                                                                 'radistar', 'massstar', 'metastar', 'loggstar', \
                                                                                 'tagetar', 'distsyst', 'numbplantranstar', 'tagestar', \
                                                                                 'tmagsyst', 'vmagsyst', \
                                                                                 
                                                                                 'radiplan', 'tmptplan', 'rvelsemapred', \
                                                                                 'pericomp', 'periplan', 'duratran', 'dcyc', 'depttrancomp', \
                                                                                 'irra', \
                                                                                 'tsmm', 'esmm', \
                                                                                 
                                                                                 'yearaler', 'toii', 's2nr', \
                                                                                 'esmmacwg', 'tsmmacwg', \
                                                                                 'numbobsvtime', 'numbobsvimag', 'numbobsvspec', \
                                                                                 
                                                                                 ]:
                boolgood = True
        
            #if True or 
            if (gdat.typeanls.startswith('exar')) and listnamefeat[k][n] in [ \
                                                                                 'declstar', 'rascstar', 'laecstar', 'loecstar', \
                                                                                 'radistar', 'massstar', 'metastar', 'loggstar', \
                                                                                 'tagetar', 'distsyst', 'numbplantranstar', 'tagestar', \
                                                                                 'vmagsyst', 'periplan', 'densplan', 'massplan', \
                                                                                 'yeardisc', \
                                                                                 'radiplan', 'tmptplan', \
                                                                                 'pericomp', 'duratran', 'dcyc', 'depttrancomp', \
                                                                                 'irra', \
                                                                                 'tsmm', 'esmm', \
                                                                                 ]:
                boolgood = True
        
            if gdat.typeanls == 'autovett':
                if not listnamefeat[k][n] in ['tic']:
                    boolgood = True
            
            if gdat.typeanls == 'plan' or gdat.typeanls == 'cosc':
                if listnamefeat[k][n] not in ['idenstar', 'booltran']:
                    boolgood = True
            
            #if gdat.typeanls == 'plan':
            #and not listnamefeat[k][n] in ['cosi', 'densstar', 'massstar', 'nois', 'pericomp', 'numbcompstarmean', 'radiplan', \
            #                                                                                                            'radistar', 'rmag', 'rsma', 'smax', 'distsyst']:
            #    continue
            #if gdat.typeanls == 'cosc' and not listnamefeat[k][n] in [ \
            #                                                                     'radicomp', 'idenstar', 'rsum', 'booloccu', 'numbcompstar', 'numbcompstarmean', \
            #                                                                  ]:
            #    boolgood = True
            
            if booldictinpt:
                boolgood = True

            if not boolgood:
                continue
        
            gdat.dictpoplfilt[gdat.listnamepopl[k]][listnamefeat[k][n]] = gdat.dictpopl[gdat.listnamepopl[k]][listnamefeat[k][n]]
            
            # exclude features with string value
            if listnamefeat[k][n] in ['typedisptess', 'strgcomm', 'namestar', 'namesyst', 'nametoii', 'nameplan', \
                                                                'facidisc', 'strgprovmass', 'strgrefrradiplan', 'strgrefrmassplan', 'name']:
                continue
            
            # exclude features with unuseful IDs
            if listnamefeat[k][n] in ['tici']:
                continue
            
            samptemp = np.array(gdat.dictpopl[gdat.listnamepopl[k]][listnamefeat[k][n]])
            if np.isfinite(samptemp).size > 0 and not isinstance(samptemp[0], str):
                listsampfilt[k].append(samptemp.astype(float))
                listnamefeatfilt[k].append(listnamefeat[k][n])
        if gdat.namefeatlablsamp is not None:
            gdat.listlablsamp[k] = gdat.dictpopl[gdat.listnamepopl[k]][gdat.namefeatlablsamp]
        
        if len(listsampfilt[k]) == 0:
            print('')
            print('')
            print('')
            print('Warning!')
            print('gdat.listnamepopl')
            print(gdat.listnamepopl)
            print('listnamefeat')
            print(listnamefeat)
            print('k')
            print(k)
            print('gdat.indxfeat[k]')
            summgene(gdat.indxfeat[k])
            print('gdat.listnamepopl[k]')
            print(gdat.listnamepopl[k])
            print('listnamefeat[k]')
            print(listnamefeat[k])
            #raise Exception('')
        else:
            listsampfilt[k] = np.vstack(listsampfilt[k]).T
    listsamp = listsampfilt
    listnamefeat = listnamefeatfilt
    
    # list of pairs of feature names to be skipped
    if 'tmptplan' in listnamefeat and 'irra' in listnamefeat:
        gdat.listnamefeatskip = [['tmptplan', 'irra']]
    else:
        gdat.listnamefeatskip = None
    
    listlablfeat = [[] for k in gdat.indxpopl]
    listscalfeat = [[] for k in gdat.indxpopl]
    for k in gdat.indxpopl:
        # get the labels and scalings for the features in the population
        listlablfeat[k], listscalfeat[k], _, _, _ = tdpy.retr_listlablscalpara(listnamefeat[k])
        
    for k in gdat.indxpopl:
        numbfeat[k] = len(listnamefeat[k])
        gdat.indxfeat[k] = np.arange(numbfeat[k])
        listnamefeat[k] = np.array(listnamefeat[k])

    # store features on disc
    print('temp: Skipping storage on the disc')
    #for k in gdat.indxpopl:
    #    path = gdat.pathdata + 'dict_%s.csv' % gdat.listnamepopl[k] 
    #    print('Writing to %s...' % path)
    #    print('gdat.dictpoplfilt[gdat.listnamepopl[k]].keys()')
    #    print(gdat.dictpoplfilt[gdat.listnamepopl[k]].keys())
    #    pd.DataFrame.from_dict(gdat.dictpoplfilt[gdat.listnamepopl[k]]).to_csv(path, header=listlablfeattotl[k], index=False, float_format='%.8g')
    
    if gdat.typelang == 'Turkish':
        for k in gdat.indxpopl:
            for m in gdat.indxfeat[k]:
                if listlablfeat[k][m][0] in gdat.dictturk:
                    listlablfeat[k][m][0] = gdat.dictturk[listlablfeat[k][m][0]]
                else:
                    print('%s is not available in gdat.dictturk...' % listlablfeat[k][m][0])
                if listlablfeat[k][m][1] in gdat.dictturk:
                    listlablfeat[k][m][1] = gdat.dictturk[listlablfeat[k][m][1]]
                else:
                    print('%s is not available in gdat.dictturk...' % listlablfeat[k][m][1])
    
    typeplottdim = 'scat'
    
    # plot populations together
    numbplotcomm = len(gdat.listdictlablcolrpopl)
    indxplotcomm = np.arange(numbplotcomm)
    gdat.listnamepoplcomm = [[] for e in indxplotcomm]
    gdat.listcolrpoplcomm = [[] for e in indxplotcomm]
    gdat.listlablpoplcomm = [[] for e in indxplotcomm]
    
    print('Comparing populations...')
    for e in indxplotcomm:
        gdat.listnamepoplcomm[e] = list(gdat.listdictlablcolrpopl[e].keys())
        
        print('e')
        print(e)
        print('gdat.listnamepoplcomm')
        print(gdat.listnamepoplcomm)
        print('gdat.listdictlablcolrpopl')
        print(gdat.listdictlablcolrpopl)
        if not gdat.listnamepoplcomm[e][0] in gdat.listnamepopl:
            print('')
            print('One of the populations in gdat.listnamepoplcomm is not available in gdat.listnamepopl')
            print('e')
            print(e)
            print('gdat.listnamepoplcomm[e][0]')
            print(gdat.listnamepoplcomm[e][0])
            print('gdat.listnamepopl')
            print(gdat.listnamepopl)
            raise Exception('')
        indxfrst = np.where(gdat.listnamepoplcomm[e][0] == gdat.listnamepopl)[0][0]
        
        if gdat.booldiag:
            for name in gdat.listnamepoplcomm[e]:
                if len(gdat.listdictlablcolrpopl[e][name]) != 2:
                    print('')
                    print('')
                    print('name')
                    print(name)
                    print('gdat.listdictlablcolrpopl[e][name]')
                    print(gdat.listdictlablcolrpopl[e][name])
                    raise Exception('')

        # translate
        if gdat.typelang == 'Turkish':
            gdat.listtitlcomp[e] = gdat.dictturk[gdat.listtitlcomp[e]]
            for name in gdat.listnamepoplcomm[e]:
                
                if gdat.listdictlablcolrpopl[e][name][0] is not None:
                    gdat.listdictlablcolrpopl[e][name][0] = gdat.dictturk[gdat.listdictlablcolrpopl[e][name][0]]

        for n in gdat.indxfeat[indxfrst]:
            for k in range(len(gdat.listnamepoplcomm[e])):
                if not gdat.listnamepoplcomm[e][k] in gdat.listnamepopl:
                    print('')
                    print('')
                    print('')
                    print('One of the populations in gdat.listnamepoplcomm is not available in gdat.listnamepopl')
                    print('gdat.listdictlablcolrpopl[e]')
                    print(gdat.listdictlablcolrpopl[e])
                    print('gdat.listnamepoplcomm[e][k]')
                    print(gdat.listnamepoplcomm[e][k])
                    print('gdat.listnamepopl')
                    print(gdat.listnamepopl)
                    raise Exception('')
        
        #for n in gdat.indxfeat[indxfrst]:
        #    print('Feature: %s' % listnamefeat[indxfrst][n])
        #    for k in range(len(gdat.listnamepoplcomm[e])):
        #        indxseco = np.where(gdat.listnamepoplcomm[e][k] == gdat.listnamepopl)[0][0]
        #        summgene(gdat.dictpopl[gdat.listnamepopl[indxseco]][listnamefeat[indxseco][n]], boolslin=True)
        #    print('')
        #print('')
        
        numbpoplcomm = len(gdat.listnamepoplcomm[e])
        if numbpoplcomm > 0:
        
            boolplotpair = gdat.listboolcompexcl[e]
            boolpoplexcl = gdat.listboolcompexcl[e]

            ## find the indices of the populations to be plotted together
            print('gdat.listnamepopl')
            print(gdat.listnamepopl)
            print('gdat.listnamepoplcomm[e]')
            print(gdat.listnamepoplcomm[e])
            indxpoplcomm = []
            for namepoplcomm in gdat.listnamepoplcomm[e]:
                indxpoplcomm.append(np.where(gdat.listnamepopl == namepoplcomm)[0][0])
            
            indxpoplcomm = np.array(indxpoplcomm)

            ## find the list of feature names common across all populations
            gdat.listnamefeatcommtemp = []
            
            gdat.listnamefeatcomm = []
            gdat.listscalfeatcomm = []
            gdat.listlablfeatcomm = []
            
            for u in indxpoplcomm:
                gdat.listnamefeatcommtemp.append(listnamefeat[u])
            gdat.listnamefeatcommtemp = np.concatenate(gdat.listnamefeatcommtemp)
            
            gdat.listnamefeatcommtemptemp = []
            for name in gdat.listnamefeatcommtemp:
                if not name in gdat.listnamefeatcommtemptemp:
                    gdat.listnamefeatcommtemptemp.append(name)
            gdat.listnamefeatcommtemp = gdat.listnamefeatcommtemptemp
            gdat.listnamefeatcommtemp = np.array(gdat.listnamefeatcommtemp)
            
            for namefeat in gdat.listnamefeatcommtemp:
                booltemp = True
                for u in indxpoplcomm:
                    if not namefeat in listnamefeat[u]:
                        booltemp = False
                if booltemp:
                    gdat.listnamefeatcomm.append(namefeat)
            
            indxfeatcomm = [[] for u in indxpoplcomm]
            for uu, u in enumerate(indxpoplcomm):
                # list of feature indcices for each population, which will be plotted together
                for namefeat in gdat.listnamefeatcomm:
                    indxfeatcomm[uu].append(np.where(listnamefeat[u] == namefeat)[0][0])
                indxfeatcomm[uu] = np.array(indxfeatcomm[uu])
            
            for indx in indxfeatcomm[0]:
                gdat.listscalfeatcomm.append(listscalfeat[indxpoplcomm[0]][indx])
                gdat.listlablfeatcomm.append(listlablfeat[indxpoplcomm[0]][indx])
                
            gdat.listsampcomm = [[] for u in indxpoplcomm]
            gdat.listlablsampcomm = None
            if gdat.listlablsamp is not None:
                gdat.listlablsampcomm = [[] for u in indxpoplcomm]
            for uu, u in enumerate(indxpoplcomm):
                print('indxfeatcomm')
                print(indxfeatcomm)
                print('uu, u')
                print(uu, u)
                print('indxfeatcomm[uu]')
                print(indxfeatcomm[uu])
                gdat.listsampcomm[uu] = listsamp[u][:, indxfeatcomm[uu]]
                if gdat.listlablsamp is not None:
                    gdat.listlablsampcomm[uu] = gdat.listlablsamp[u]

                if gdat.booldiag:
                    if gdat.listlablsamp is not None and len(gdat.listlablsampcomm[uu]) == 0:
                        print('')
                        print('')
                        print('')
                        print('')
                        print('uu, u')
                        print(uu, u)
                        print('gdat.listlablsampcomm[uu]')
                        print(gdat.listlablsampcomm[uu])
                        print('gdat.listlablsamp')
                        print(gdat.listlablsamp)
                        print('gdat.listlablsamp[u]')
                        print(gdat.listlablsamp[u])
                        raise Exception('')
            
            # number of samples in each population
            numbsamppopl = np.empty(numbpoplcomm, dtype=int)
            for uu in range(numbpoplcomm):
                numbsamppopl[uu] = gdat.listsampcomm[uu].shape[0]
            
            # determine the list of labels for compared populations
            for uu, u in enumerate(indxpoplcomm):
                gdat.listlablpoplcomm[e].append(gdat.listdictlablcolrpopl[e][gdat.listnamepopl[u]][0])
            
            # determine the list of colors for compared populations
            for uu, u in enumerate(indxpoplcomm):
                gdat.listcolrpoplcomm[e].append(gdat.listdictlablcolrpopl[e][gdat.listnamepopl[u]][1])
            
            for uu in range(numbpoplcomm):
                print('%d samples in %s.' % (numbsamppopl[uu], gdat.listnamepoplcomm[e][uu]))
            print('')

            strgplot = ''
            for uu, u in enumerate(indxpoplcomm):
                strgplot += gdat.listnamepoplcomm[e][uu]
            
            print('e')
            print(e)
            print('gdat.listnamepoplcomm[e]')
            print(gdat.listnamepoplcomm[e])
            print('gdat.listcolrpoplcomm[e]')
            print(gdat.listcolrpoplcomm[e])
            print('gdat.listtitlcomp')
            print(gdat.listtitlcomp)
            print('gdat.listtitlcomp[e]')
            print(gdat.listtitlcomp[e])
            print('')

            for m in range(2):
                
                if m == 0:
                    pathbase = gdat.pathimag
                    boolmakelegd = True
                else:
                    pathbase = gdat.pathimag + 'without_legend/'
                    boolmakelegd = False
                
                tdpy.plot_grid( \
                               gdat.listlablfeatcomm, \
                               listpara=gdat.listsampcomm, \
                               strgplot=strgplot, \
                               pathbase=pathbase, \
                               boolplothistodim=True, boolplotpair=boolplotpair, \
                               boolpoplexcl=boolpoplexcl, \
                               boolplottria=False, \
                               listnamefeatcumu=listnamefeatcumu, \
                               typeplottdim=typeplottdim, \
                               typefileplot=typefileplot, \
                               lablnumbsamp=gdat.lablnumbsamp, \
                               listlablsamp=gdat.listlablsampcomm, \
                               lablsampgene=gdat.lablsampgene, \
                               listnamefeatskip=gdat.listnamefeatskip, \
                               boolmakelegd=boolmakelegd, \
                               listnameordrpair=gdat.listnameordrpair, \
                               listnamepara=gdat.listnamefeatcomm, \
                               listscalpara=gdat.listscalfeatcomm, \
                               listlablpopl=gdat.listlablpoplcomm[e], \
                               listcolrpopl=gdat.listcolrpoplcomm[e], \
                               titl=gdat.listtitlcomp[e], \
                              )
            
        
        if gdat.typeanls == 'cosc' or gdat.typeanls == 'plan':
            # names of the population to be used as relevant and positive population
            dictnamepopl = dict()
            dictnamepopl['rele'] = 'comp' + strgpoplstar + 'tran'
            dictnamepopl['posi'] = 'comp' + strgpoplstar + 'tranposi'
            
            #boolposirele = np.zeros(gdat.dictnumbsamp[namepoplcomptotl], dtype=bool)
            #boolposirele[gdat.dictindxsamp[namepoplcomptotl][namepoplrele]] = True
            #boolreleposi = np.ones(gdat.dictnumbsamp[namepoplrele], dtype=bool)
            
            # construct arrays of parameters for each type of population
            #dictlistvarb = dict()
            #for strgtypepopl in ['rele', 'anls']:
            #    listvarb = []
            #    for name in dictlistnamevarb['rele']:
            #        dictlistvarb[strgtypepopl] += [gdat.dictpopl[dictnamepopl[strgtypepopl]][name]]
            #    dictlistvarb[strgtypepopl] = np.vstack(dictlistvarb[strgtypepopl])
            
            tdpy.plot_recaprec(gdat.pathimag, dictnamepopl, gdat.dictpopl, \
                            #namepoplcomptotl, listvarbrele, listvarbanls, listnamevarbrele, listnamevarbanls, \
                                                         #listlablvarbrele, listlablvarbanls, \
                                                         #boolposirele=boolposirele, \
                                                         #boolreleposi=boolreleposi, \
                                                         
                                                         strgreca='Transit Fraction', \
                                                         )
            
            # ROC of detections
            #listnamevarbanls = ['sdee']
            #listlablvarbanls = [['SDE', '']]
            #boolposirele = np.zeros(gdat.dictnumbsamp[namepoplcomptran], dtype=bool)
            #boolposirele[gdat.dictindxsamp[namepoplcomptran][namepoplcomptranposi]] = True
            ##boolreleposi = np.ones(gdat.dictindxsamp[namepoplcomptran][namepoplcomptranposi].size, dtype=bool)
            #boolreleposi = np.ones(gdat.dictnumbsamp[namepoplcomptranposi], dtype=bool)
            #listvarbrele = np.vstack([gdat.dictpopl[namepoplcomptotl]['masscomp'], gdat.dictpopl[namepoplcomptotl]['pericomp'], gdat.dictpopl[namepoplcomptotl]['inclcomp']]).T
            #listvarbprec = np.vstack([gdat.dictpopl[namepoplcomptranposi]['sdee']]).T
            #tdpy.plot_recaprec(gdat.pathimag, namepoplcomptranposi, listvarbrele, listvarbprec, listnamevarbrele, listnamevarbprec, \
            #                                                                                listlablvarbrele, listlablvarbprec, boolposirele, boolreleposi)
            #
            ## overall selection
            #listnamevarbprec = ['sdee']
            #listlablvarbprec = [['SDE', '']]
            #boolposirele = np.zeros(gdat.dictnumbsamp[namepoplcomptranposi], dtype=bool)
            #boolposirele[indxtranoccu[gdat.dictindxsamp[namepoplcomptran][namepoplcomptranposi]]] = True
            #boolreleposi = np.ones(gdat.dictindxsamp[namepoplcomptran][namepoplcomptranposi].size, dtype=bool)
            #listvarbrele = np.vstack([gdat.dictpopl[namepoplcomptotl]['masscomp'], gdat.dictpopl[namepoplcomptotl]['pericomp'], gdat.dictpopl[namepoplcomptotl]['inclcomp']]).T
            #listvarbprec = np.vstack([gdat.dictpopl[namepoplcomptranposi]['sdee']]).T
            #tdpy.plot_recaprec(gdat.pathimag, 'totl', listvarbrele, listvarbprec, listnamevarbrele, listnamevarbprec, \
            #                                                             listlablvarbrele, listlablvarbprec, boolposirele, boolreleposi)
            #

            ## overall selection
            #listnamevarbprec = listnamevarbrele
            #listlablvarbprec = listlablvarbrele
            #boolposirele = np.zeros(gdat.dictnumbsamp['comp'], dtype=bool)
            #boolposirele[indxtranoccu[gdat.dictindxsamp[namepoplcomptran][namepoplcomptranposi]]] = True
            #boolreleposi = np.ones(gdat.dictindxsamp[namepoplcomptran][namepoplcomptranposi].size, dtype=bool)
            #listvarbrele = np.vstack([gdat.dictpopl[namepoplcomptotl]['masscomp'], gdat.dictpopl[namepoplcomptotl]['pericomp'], gdat.dictpopl[namepoplcomptotl]['inclcomp']]).T
            #listvarbprec = np.vstack([gdat.dictpopl[namepoplcomptranposi]['masscomp'], gdat.dictpopl[namepoplcomptranposi]['pericomp'], gdat.dictpopl[namepoplcomptranposi]['inclcomp']]).T
            #listvarbdete = []
            #tdpy.plot_recaprec(gdat.pathimag, 'comp', listvarbrele, listvarbprec, listnamevarbrele, listnamevarbprec, \
            #                                                             listlablvarbrele, listlablvarbprec, boolposirele, boolreleposi, listvarbdete=listvarbdete)
            

    # derived features
    # effective temperature of the star
    #tmptstar = arry[:, 0]
    #tmptstarstdv = arry[:, 1]
    ## metallicity of the star
    #metastar = arry[:, 2]
    #metastarstdv = arry[:, 3]
    ## radius of the star
    #radistar = arry[:, 4]
    #radistarstdv = arry[:, 5]
    ## mass of the planet
    #massplan = arry[:, 6]
    #massplanstdv = arry[:, 7]
    ## semi-major axis divided by the radius of the star
    #rsta = arry[:, 8]
    #rstastdv = arry[:, 9]
    ## radius of the planet divided by the radius of the star
    #rrat = arry[:, 10]
    #rratstdv = arry[:, 11]
    ## dayside temperature of the planet
    #tday = arry[:, 12]
    #tdaystdv = np.maximum(arry[:, 13], arry[:, 14])
    ## geometric albedo 
    #albe = arry[:, 15]
    #albestdv = np.maximum(arry[:, 16], arry[:, 17])
        
    ## semi-major axis
    #smax = radistar / rsta
    #smaxstdv = np.sqrt((radistarstdv / radistar)**2 + (rstastdv / rsta)**2)
    ## radius of the planet
    #radiplan = rrat * radistar
    #radiplanstdv = np.sqrt((rratstdv / rrat)**2 + (radistarstdv / radistar)**2)
    ## plenatery surface gravity
    #grav = 26.18 * massplan / radiplan**2
    #stdvgrav = grav * np.sqrt((massplanstdv / massplan) **2 + (2. * radiplanstdv / radiplan)**2)
    ## log of the plenatery surface gravity
    #logg = np.log10(grav)
    #loggstdv = 0.434 * stdvgrav / grav
    ## equilibrium temperature of the planet
    #tmpteqiu = tmptstar * np.sqrt(1. / 2. / smax)
    #tmpteqiustdv = tmpteqiu * np.sqrt((tmptstarstdv / tmptstar)**2 + (0.5 * smaxstdv / smax)**2) / np.sqrt(2.)
    #
    ## load variables into the array
    #arrytemp = np.empty((numbplan, numbfeat))
    #arrystdvtemp = np.empty((numbplan, numbfeat))
    ## logg 
    #arrytemp[:, 0] = logg
    #arrystdvtemp[:, 0] = loggstdv
    ## plenatery equilibrium temperature
    #arrytemp[:, 1] = tmpteqiu
    #arrystdvtemp[:, 1] = tmpteqiustdv
    ## stellar metallicity
    #arrytemp[:, 2] = metastar
    #arrystdvtemp[:, 2] = metastarstdv
    ## dayside temperature
    #arrytemp[:, 3] = tday
    #arrystdvtemp[:, 3] = tdaystdv
    ## dayside albedo
    #arrytemp[:, 4] = albe
    #arrystdvtemp[:, 4] = albestdv
        
    # search for correlations
    if boolsrchcorr:

        # sampler setup
        numbsampburnwalk = 1000
        numbsampwalk = 10000
        numbsampburnwalkseco = 2000
    
        for k in gdat.indxpopl:
            
            print('Searching for correlations...')    

            listlablpara = [[r'$\alpha$', ''], [r'$\rho$', '']]
            listscalpara = ['self', 'self']
            listmeangauspara = None
            liststdvgauspara = None
            listminmpara = np.array([0, -1e1])
            listmaxmpara = np.array([np.pi, 1e1])
            
            numbpara = len(listlablpara)
            numbwalk = max(20, 2 * numbpara)
            numbsamp = numbwalk * (numbsampwalk - numbsampburnwalkseco)
            indxsamp = np.arange(numbsamp)
            
            listcoef = np.zeros((numbfeat[k], numbfeat[k]))
            
            for n in gdat.indxfeat[k]: 
                gdat.tempfrst = gdat.dictpopl[gdat.listnamepopl[k]][listnamefeat[k][n]]
                gdat.tempfrststdv = arrystdvtemp[:, n]
                for m in gdat.indxfeat[k]: 
                    gdat.tempseco = gdat.dictpopl[gdat.listnamepopl[k]][listnamefeat[k][m]]
                    gdat.tempsecostdv = arrystdvtemp[:, m]
                    
                    minmxpos = np.amin(gdat.tempfrst) / 1.1
                    maxmxpos = np.amax(gdat.tempfrst) * 1.1
                    minmypos = np.amin(gdat.tempseco) / 1.1
                    maxmypos = np.amax(gdat.tempseco) * 1.1
                    
                    # calculate PCC
                    coef, pval = scipy.stats.pearsonr(gdat.tempfrst, gdat.tempseco)
                    listcoef[u, k] = coef
                    
                    # sample from a linear model
                    numbdata = 2 * numbplan
                    strgextn = '%d_' % b + listnamepara[k]
                    featpost = tdpy.samp(gdat, pathimag, numbsampwalk, numbsampburnwalk, numbsampburnwalkseco, retr_llik, \
                                                    listlablpara, listscalpara, listminmpara, listmaxmpara, listmeangauspara, liststdvgauspara, \
                                                        numbdata, strgextn=strgextn, strgplotextn=gdat.strgplotextn, verbtype=0)
                    
                    figr, axis = plt.subplots(figsize=(4, 4))
                    xerr = gdat.tempfrststdv
                    yerr = gdat.tempsecostdv
                    
                    if b == 0:
                        colr = 'b'
                    if b == 1:
                        colr = 'k'
                    if b == 0 or b == 1:
                        axis.errorbar(gdat.tempfrst, gdat.tempseco, yerr=yerr, xerr=xerr, fmt='o', color=colr)
                    if b == 2:
                        numbplantess = len(data[0])
                        axis.errorbar(gdat.tempfrst[:numbplantess], gdat.tempseco[:numbplantess], \
                                                                                yerr=yerr[:numbplantess], xerr=xerr[:numbplantess], fmt='o', color='b')
                        axis.errorbar(gdat.tempfrst[numbplantess:], gdat.tempseco[numbplantess:], \
                                                                                yerr=yerr[numbplantess:], xerr=xerr[numbplantess:], fmt='o', color='k')
                    
                    axis.set_xlim([minmxpos, maxmxpos])
                    axis.set_ylim([minmypos, maxmypos])
                    
                    postslop = -1. / np.tan(featpost[:, 0])
                    medislop = np.median(postslop)
                    lowrslop = np.median(postslop) - np.percentile(postslop, 16)
                    upprslop = np.percentile(postslop, 84) - np.median(postslop)
                    titl = 'PCC = %.3g, Slope: %.3g $\substack{+%.2g \\\\ -%.2g}$' % (listcoef[u, k], medislop, upprslop, lowrslop)
                    axis.set_xlabel(listlabl[k])
                    axis.set_ylabel(listlabl[u])
            
                    plt.tight_layout()
                    path = pathimag + 'scat_%s_%s_%d.%s' % (listnamefeat[k], listnamefeat[u], b, gdat.strgplotextn)
                    print('Writing to %s...' % path)
                    plt.savefig(path)
                    plt.close()

    # PCA
    if gdat.boolprca:
        import mergen
        mergen.exec_prca_skit(listvarb)

    # search for clusters

    if gdat.boolintpscat:
        ## interpolate distribution
        np.interp2d(grid)
    
    print('')
    print('')
    print('')
    print('')
    print('')

