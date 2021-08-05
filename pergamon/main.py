import os
import sys

import numpy as np
#import scipy.stat
import pandas as pd
import matplotlib.pyplot as plt

import ephesus
from tdpy import gdatstrt
import tdpy
from tdpy import summgene 
import miletos
import pcat

def retr_modl_corr(gdat, feat, inpt):
    
    angl = feat[0]
    gamm = feat[1]

    slop = -1. / np.tan(angl)
    intc = gamm / np.sin(angl)
    
    line = inpt * slop + intc
    
    return line, []


def retr_llik_corr(feat, gdat):
    
    angl = feat[0]
    gamm = feat[1]

    dist = np.cos(angl) * gdat.tempfrst + np.sin(angl) * gdat.tempseco - gamm
    vari = np.cos(angl)**2 * gdat.tempfrststdv**2 + np.sin(angl)**2 * gdat.tempsecostdv**2

    llik = -0.5 * np.sum(dist**2 / vari)

    return llik


def init( \
        # type of analysis
        typeanly, \

        # plotting
        typefileplot='pdf', \
        **args, \
        ):
    
    """
    visualize the population
    search for correlations
    search for clusters
    model the density and selection effects to determine the occurence rate
    """
    # construct global object
    gdat = tdpy.gdatstrt()
    
    # copy unnamed inputs to the global object
    for attr, valu in locals().items():
        if '__' not in attr and attr != 'gdat':
            setattr(gdat, attr, valu)

    # copy named arguments to the global object
    for strg, valu in args.items():
        setattr(gdat, strg, valu)

    # paths
    gdat.pathbase = os.environ['PERGAMON_DATA_PATH'] + '/'
    gdat.pathdata = gdat.pathbase + 'data/'
    gdat.pathimag = gdat.pathbase + 'imag/'
    os.system('mkdir -p %s' % gdat.pathdata)
    os.system('mkdir -p %s' % gdat.pathimag)
    
    # plotting
    strgplotextn = 'pdf'

    # list of population names
    ## 2-minute cadence target list of TESS during the nominal mission
    #if gdat.typeanly == 'exopmock2minnomi':
    #gdat.listnamepopl = ['2minnomi', 'exopmock2minnomi']
    
    ## 'bcanmock': mock BH candidate subpopulation of '2minnomi'
    ## '2minnomi': 2-min cadence target list during the nominal mission
    #if gdat.typeanly == 'bcanmock2minnomi':
    #gdat.listnamepopl = ['2minnomi', 'bcanmock2minnomi']
    ## 'exopmock': mock exoplanet subpopulation of '2minnomi'
    ## 'exoptoii': TOIs
    
    print('gdat.typeanly')
    print(gdat.typeanly)

    # correlation search
    numbsampfeww = 1000
    numbsampburnwalk = 1000
    numbsampwalk = 10000
    numbsampburnwalkseco = 2000
    
    ## conversion factors
    dictfact = ephesus.retr_factconv()
    
    dictpopl = dict()

    # preliminary setup for analyses
    if gdat.typeanly == 'mocklsst':
        # time stamps
        listtime = np.random.rand(1000)
    
        ## number of years into the mission
        numbyear = 10
        numbsamp = int(numbyear * 100)
    
        minmtime = 0.
        maxmtime = 30.
    
        # cadence
        cade = 2. / 60. / 24. # [days]
    
        numbsamp = 1
        indxsamp = np.arange(numbsamp)
    
    if gdat.typeanly == 'tess135nomiexopfstr':
        
        dictpopl['fstr'] = dict()
        dictpopl['fstrastr'] = dict()
        
        pathfstr = '/Users/tdaylan/Documents/work/data/external/FaintStars/'
        path = pathfstr + 'metr/metr_sc01.txt'
        print('Reading from %s...' % path)
        objt = open(path, 'r')
        k = 0
        ticifstr = []
        for line in objt:
            linesplt = line.split(' ')
            print('k')
            print(k)
            if k == 0:
                listnamefstr = linesplt
                #['peri', 'dura', 'dept', 'dcyctran', 'dcycinge', 'SN', 'SPN', 'OOTmag', 'sigamp']
                for namefstr in listnamefstr:
                    dictpopl['fstr'][namefstr] = []
            
                print('len(listnamefstr)')
                print(len(listnamefstr))
            else:

                #peri = float(linesplt[1])
                #if peri < 0 or peri > 1000:
                #    continue
                #
                #dept = float(linesplt[4])
                #if dept <= 0:
                #    continue
                #
                #dura = float(linesplt[3])
                #if dura <= 0:
                #    continue
                #
                #dcyctran = float(linesplt[5])
                #if dcyctran <= 0:
                #    continue
                #
                #boolcont = False
                #for p in range(11):
                #    if float(linesplt[p]) <= 0:
                #        boolcont = True
                #if boolcont:
                #    continue
                #    
                #dcycinge = float(linesplt[6])
                #if dcycinge <= 0:
                #    continue
                
                numbfiel = len(linesplt)
                ticifstr.append(float(linesplt[0]))
                print('numbfiel')
                print(numbfiel)
                for n in range(1, numbfiel):
                    print('n')
                    print(n)
                    print('linesplt[n]')
                    print(linesplt[n])
                    dictpopl['fstr'][listnamefstr[n]].append(float(linesplt[n]))
                    
                #dictpopl['fstr']['peri'].append(float(linesplt[1]))
                #dictpopl['fstr']['dura'].append(12. * float(linesplt[3]))
                #dictpopl['fstr']['dept'].append(float(linesplt[4]) * 1e3)
                #dictpopl['fstr']['dcyctran'].append(float(linesplt[5]))
                #dictpopl['fstr']['dcycinge'].append(float(linesplt[6]) * float(linesplt[5]))
                #dictpopl['fstr']['SN'].append(float(linesplt[7]))
                #dictpopl['fstr']['SPN'].append(float(linesplt[8]))
                #dictpopl['fstr']['OOTmag'].append(float(linesplt[9]))
                #dictpopl['fstr']['sigamp'].append(float(linesplt[10]))
                
            k += 1
        ticifstr = np.array(ticifstr) 
        for namefstr in listnamefstr:
            dictpopl['fstr'][namefstr] = np.array(dictpopl['fstr'][namefstr]) 
        #dictpopl['fstr'] = pd.read_csv(path, delimiter=' ').to_dict(orient='list')

        #dictquer = miletos.xmat_tici(ticifstrastr)
        
        # read TIC IDs and dispositions from each sector
        ## TESS sectors
        #listtsec = range(1, 27)
        listtsec = range(1, 2)
        numbtsec = len(listtsec)
        indxtsec = np.arange(numbtsec)
        ticifull = [[] for o in indxtsec]
        dispfull = [[] for o in indxtsec]
        pathbase = pathfstr + 'combined/'
        for o in indxtsec:
            path = pathbase + 'sector-%d.tier2.csv' % listtsec[o]
            print('Reading from %s...' % path)
            dictquer = pd.read_csv(path).to_dict(orient='list')
            ticifull[o] = dictquer['TIC']
            dispfull[o] = dictquer['Final']
        ## merge across sectors
        ticifullconc = np.concatenate(ticifull)
        dispfullconc = np.concatenate(dispfull)
        
        print('ticifullconc')
        summgene(ticifullconc)

        ticifstrastr, indxuniq, indxinve, numbtici = np.unique(ticifullconc, return_index=True, return_inverse=True, return_counts=True)
        dispfstrastr = ticifullconc[indxuniq]

        #ticifstrastr = ticifstrastr.astype(int)
        
        numbastr = len(ticifstrastr)
        for namefstr in listnamefstr:
            dictpopl['fstrastr'][namefstr] = np.empty(numbastr)
        for k, tici in enumerate(ticifstrastr):
            indx = np.where(ticifstr == tici)[0]
            if indx.size > 0:
                for namefstr in listnamefstr:
                    dictpopl['fstrastr'][namefstr][k] = dictpopl['fstr'][namefstr][indx[0]]
        
        figr, axis = plt.subplots(figsize=(4, 4))
        axis.hist(numbtici)
        axis.set_yscale('log')
        axis.set_xlabel('Number of sectors')
        axis.set_ylabel('N')
        plt.tight_layout()
        path = gdat.pathimag + 'histnumbtsec.%s' % (strgplotextn)
        print('Writing to %s...' % path)
        plt.savefig(path)
        plt.close()

    if gdat.typeanly == 'hjuppcur':
        listlabl = ['log $g_{\mathrm{p}}$ [cgs]', r'$T_{\mathrm{eq}}$ [K]', '[Fe/H] [dex]', r'$T_{\mathrm{day}}$ [K]', '$A_g$']
        listname = ['logg', 'ptmp', 'meta', 'tday', 'albe']
        path = gdat.pathdata + 'catlpcurtess.csv'
        # get the dictionaries holding the population properties
        dictpopl['hjuppcur'] = pd.read_csv(path)

    if gdat.typeanly == 'exoptoyy':
        listradistar = tdpy.icdf_powr(np.random.rand(numbsamp), 0.1, 10., -2.)
        listmassstar = tdpy.icdf_powr(np.random.rand(numbsamp), 0.1, 10., -2.)
        
    if gdat.typeanly == 'exopexar':
        dictpopl['exopexar'] = retr_dictexar()
        dictpopl['exopexar']['nois'] = ephesus.retr_noistess(dictpopl['exopexar']['vmagsyst'])
    
    if gdat.typeanly == 'pcanfstr':
        indxfullpcan = np.where(dispfull == 'PC')[0]
        dictpopl['tici'] = ticifull[indxfullpcan]
    
    #if gdat.typeanly == 'tessnomi2minexopmocktoyy':
    if gdat.typeanly[12:].startswith('bcanmock') or gdat.typeanly[12:].startswith('exopmock'):
        
        # name of the star population
        namepoplstar = gdat.typeanly[:12]
        print('namepoplstar')
        print(namepoplstar)
        
        # get the features of the star population
        dictpopl[namepoplstar] = miletos.retr_dictcatltic8(namepoplstar)
        dictpopl[namepoplstar]['densstar'] = 1.41 * dictpopl[namepoplstar]['massstar'] / dictpopl[namepoplstar]['radistar']**3
        numbtargtotl = dictpopl[namepoplstar]['tmag'].size
        print('numbtargtotl')
        print(numbtargtotl)
        
        # calculate TESS photometric precision for the star population
        dictpopl[namepoplstar]['nois'] = ephesus.retr_noistess(dictpopl[namepoplstar]['tmag'])

        # probability of occurence
        #dictpopl[namepoplstar]['proboccu'] = pcat.icdf_gaustrun(np.random.rand(dictpopl[namepoplstar]['radistar'].size), 0.02, 0.002, 0, np.inf)
        
        # Boolean flag of occurence
        #dictpopl[namepoplstar]['booloccu'] = np.random.rand(dictpopl[namepoplstar]['radistar'].size) < dictpopl[namepoplstar]['proboccu']
        
        # stars with occurence
        namepoploccu = gdat.typeanly[:20]
        
        minmmassbhol = 5.
        maxmmassbhol = 200.
        #binsmassbhol = np.logspace(minmmassbhol, maxmmassbhol, 11)
        #meanmassbhol = (binsmassbhol[1:] + binsmassbhol[:-1]) / 2.
        #binsperi = np.logspace(minmperi, maxmperi, 11)
        #meanperi = (binsperi[1:] + binsperi[:-1]) / 2.
                
        dictpopl[namepoploccu] = dict()
        #indx = np.where(dictpopl[namepoplstar]['booloccu'])[0]
        for name in dictpopl[namepoplstar].keys():
            #dictpopl[namepoploccu][name] = dictpopl[namepoplstar][name][indx]
            dictpopl[namepoploccu][name] = dictpopl[namepoplstar][name]
        numbtargoccu = dictpopl[namepoploccu]['radistar'].size
        print('numbtargoccu') 
        print(numbtargoccu)

        dictpopl[namepoploccu]['incl'] = np.random.rand(numbtargoccu) * 90.
        dictpopl[namepoploccu]['cosi'] = np.cos(np.pi / 180. * dictpopl[namepoploccu]['incl'])
        
        minmperi = 0.3
        maxmperi = 1000.
        dictpopl[namepoploccu]['peri'] = np.empty(numbtargoccu)
        print('dictpopl[namepoploccu][radistar][k]')
        summgene(dictpopl[namepoploccu]['radistar'])
        print('dictpopl[namepoploccu][massstar][k]')
        summgene(dictpopl[namepoploccu]['massstar'])
        print('dictpopl[namepoploccu][densstar][k]')
        summgene(dictpopl[namepoploccu]['densstar'])
        for k in range(numbtargoccu):
            if np.isfinite(dictpopl[namepoploccu]['densstar'][k]):
                dens = dictpopl[namepoploccu]['densstar'][k]
            else:
                dens = 1.
            minmperi = 0.43 / dens
            dictpopl[namepoploccu]['peri'][k] = tdpy.util.icdf_powr(np.random.rand(), minmperi, maxmperi, 2.)
        
        #indx = np.where(dictpopl[namepoploccu]['peri'] < 0.43 / dictpopl[namepoploccu]['densstar'])[0]
        #while indx.size > 0:
        #    dictpopl[namepoploccu]['peri'][indx] = tdpy.util.icdf_powr(np.random.rand(indx.size), minmperi, maxmperi, 2.)
        #    indx = np.where(dictpopl[namepoploccu]['peri'] < 0.43 / dictpopl[namepoploccu]['densstar'])[0]

        dictpopl[namepoploccu]['radiplan'] = tdpy.util.icdf_powr(np.random.rand(numbtargoccu), 1., 23., 2.)
        if gdat.typeanly[12:].startswith('bcanmock'):
            dictpopl[namepoploccu]['massbhol'] = tdpy.util.icdf_powr(np.random.rand(numbtargoccu), minmmassbhol, maxmmassbhol, 2.)
            dictpopl[namepoploccu]['masstotl'] = dictpopl[namepoploccu]['massbhol'] + dictpopl[namepoploccu]['massstar']
        if gdat.typeanly[12:].startswith('exopmock'):
            dictpopl[namepoploccu]['massplan'] = tdpy.util.icdf_powr(np.random.rand(numbtargoccu), 5., 200., 2.)
            dictpopl[namepoploccu]['masstotl'] = dictpopl[namepoploccu]['massplan'] / dictfact['msme'] + dictpopl[namepoploccu]['massstar']
        dictpopl[namepoploccu]['smax'] = ephesus.retr_smaxkepl(dictpopl[namepoploccu]['peri'], dictpopl[namepoploccu]['masstotl'])
        dictpopl[namepoploccu]['rsma'] = dictpopl[namepoploccu]['radistar'] / dictpopl[namepoploccu]['smax'] / dictfact['aurs']
        
        # stars with transiting occurence
        namepoploccutran = namepoploccu + 'tran'
        dictpopl[namepoploccutran] = dict()
        indxtranoccu = np.where(dictpopl[namepoploccu]['rsma'] > dictpopl[namepoploccu]['cosi'])[0]
        for name in dictpopl[namepoploccu].keys():
            dictpopl[namepoploccutran][name] = dictpopl[namepoploccu][name][indxtranoccu]
        numbtargoccutran = dictpopl[namepoploccutran]['radistar'].size
        print('numbtargoccutran') 
        print(numbtargoccutran)
        
        # transit duration
        dictpopl[namepoploccutran]['duratran'] = ephesus.retr_duratran(dictpopl[namepoploccutran]['peri'], \
                                                                       dictpopl[namepoploccutran]['rsma'], \
                                                                       dictpopl[namepoploccutran]['cosi'])
        dictpopl[namepoploccutran]['dcyc'] = dictpopl[namepoploccutran]['duratran'] / dictpopl[namepoploccutran]['peri'] / 24.
        if gdat.typeanly[12:].startswith('exopmock'):
            dictpopl[namepoploccutran]['rrat'] = dictpopl[namepoploccutran]['radiplan'] / dictpopl[namepoploccutran]['radistar'] / dictfact['rsre']
            dictpopl[namepoploccutran]['dept'] = 1e3 * dictpopl[namepoploccutran]['rrat']**2 # [ppt]
        if gdat.typeanly[12:].startswith('bcanmock'):
            dictpopl[namepoploccutran]['amplslen'] = ephesus.retr_amplslen(dictpopl[namepoploccutran]['peri'], dictpopl[namepoploccutran]['radistar'], \
                                                                                dictpopl[namepoploccutran]['massbhol'], dictpopl[namepoploccutran]['massstar'])
        
        # detection
        print('temp -- assuming all targets have one sector')
        dictpopl[namepoploccutran]['numbtsec'] = np.ones(numbtargoccutran)
        dictpopl[namepoploccutran]['numbtran'] = 27.3 * dictpopl[namepoploccutran]['numbtsec'] / dictpopl[namepoploccutran]['peri']
        ## SNR
        if gdat.typeanly[12:].startswith('exopmock'):
            dictpopl[namepoploccutran]['sdee'] = np.sqrt(dictpopl[namepoploccutran]['duratran']) * dictpopl[namepoploccutran]['dept'] / dictpopl[namepoploccutran]['nois']
        if gdat.typeanly[12:].startswith('bcanmock'):
            dictpopl[namepoploccutran]['sdee'] = np.sqrt(dictpopl[namepoploccutran]['duratran']) * dictpopl[namepoploccutran]['amplslen'] / dictpopl[namepoploccutran]['nois']
        dictpopl[namepoploccutran]['booldete'] = dictpopl[namepoploccutran]['sdee'] > 7
    
        # stars with transiting occurence that are detected
        namepoploccutrandete = namepoploccutran + 'dete'
        dictpopl[namepoploccutrandete] = dict()
        indxdetetran = np.where(dictpopl[namepoploccutran]['booldete'])[0]
        for name in dictpopl[namepoploccutran].keys():
            dictpopl[namepoploccutrandete][name] = dictpopl[namepoploccutran][name][indxdetetran]
        
        numbtargdete = indxdetetran.size
        print('numbtargdete') 
        print(numbtargdete)
        
        # stars with detections
        #namepopldete = namepoplstar + 'dete'
        #dictpopl[namepopldete] = dict()
        #indx = np.where(dictpopl[namepoplstar]['booldete'])[0]
        #for name in dictpopl[namepoplstar].keys():
        #    dictpopl[namepoploccutrandete][name] = dictpopl[namepoplstar][name][indx]
        
        # variables
        liststrgvarbreca = ['massbhol', 'peri', 'incl']
        listlablvarbreca = [['$M_c$', '$M_\odot$'], ['$P$', 'day'], ['$i$', '$^o$']]
        
        print('dictpopl[namepoploccu][peri]')
        summgene(dictpopl[namepoploccu]['peri'])
        
        # selection imposed by transit requirement
        liststrgvarbprec = None
        listlablvarbprec = None
        boolposirele = np.zeros(numbtargoccu, dtype=bool)
        boolposirele[indxtranoccu] = True
        boolreleposi = np.ones(numbtargoccutran, dtype=bool)
        listvarbreca = np.vstack([dictpopl[namepoploccu]['massbhol'], dictpopl[namepoploccu]['peri'], dictpopl[namepoploccu]['incl']]).T
        listvarbprec = None
        tdpy.plot_recaprec(gdat.pathimag, 'tran', listvarbreca, listvarbprec, liststrgvarbreca, liststrgvarbprec, \
                                                     listlablvarbreca, listlablvarbprec, boolposirele, boolreleposi, strgreca='Transit Fraction')
        
        #raise Exception('') 
        # ROC of detections
        liststrgvarbprec = ['sdee']
        listlablvarbprec = [['SDE', '']]
        boolposirele = np.zeros(numbtargoccutran, dtype=bool)
        boolposirele[indxdetetran] = True
        boolreleposi = np.ones(numbtargdete, dtype=bool)
        listvarbreca = np.vstack([dictpopl[namepoploccutran]['massbhol'], dictpopl[namepoploccutran]['peri'], dictpopl[namepoploccutran]['incl']]).T
        listvarbprec = np.vstack([dictpopl[namepoploccutrandete]['sdee']]).T
        tdpy.plot_recaprec(gdat.pathimag, 'dete', listvarbreca, listvarbprec, liststrgvarbreca, liststrgvarbprec, \
                                                                                        listlablvarbreca, listlablvarbprec, boolposirele, boolreleposi)
        
        # overall selection
        liststrgvarbprec = ['sdee']
        listlablvarbprec = [['SDE', '']]
        boolposirele = np.zeros(numbtargoccu, dtype=bool)
        boolposirele[indxtranoccu[indxdetetran]] = True
        boolreleposi = np.ones(numbtargdete, dtype=bool)
        listvarbreca = np.vstack([dictpopl[namepoploccu]['massbhol'], dictpopl[namepoploccu]['peri'], dictpopl[namepoploccu]['incl']]).T
        listvarbprec = np.vstack([dictpopl[namepoploccutrandete]['sdee']]).T
        tdpy.plot_recaprec(gdat.pathimag, 'totl', listvarbreca, listvarbprec, liststrgvarbreca, liststrgvarbprec, \
                                                                     listlablvarbreca, listlablvarbprec, boolposirele, boolreleposi)
        

        # overall selection
        liststrgvarbprec = liststrgvarbreca
        listlablvarbprec = listlablvarbreca
        boolposirele = np.zeros(numbtargoccu, dtype=bool)
        boolposirele[indxtranoccu[indxdetetran]] = True
        boolreleposi = np.ones(numbtargdete, dtype=bool)
        listvarbreca = np.vstack([dictpopl[namepoploccu]['massbhol'], dictpopl[namepoploccu]['peri'], dictpopl[namepoploccu]['incl']]).T
        listvarbprec = np.vstack([dictpopl[namepoploccutrandete]['massbhol'], dictpopl[namepoploccutrandete]['peri'], dictpopl[namepoploccutrandete]['incl']]).T
        listvarbdete = []
        tdpy.plot_recaprec(gdat.pathimag, 'occu', listvarbreca, listvarbprec, liststrgvarbreca, liststrgvarbprec, \
                                                                     listlablvarbreca, listlablvarbprec, boolposirele, boolreleposi, listvarbdete=[])
        

    gdat.listnamepopl = list(dictpopl.keys())
    print('gdat.listnamepopl')
    print(gdat.listnamepopl)

    gdat.numbpopl = len(gdat.listnamepopl)
    gdat.indxpopl = np.arange(gdat.numbpopl)

    listlablfeat = [[] for k in gdat.indxpopl]
    listscalfeat = [[] for k in gdat.indxpopl]
    listnamefeat = [[] for k in gdat.indxpopl]
    gdat.indxfeat = [[] for k in gdat.indxpopl]
    
    # make sure fields with strings are not included
    #if gdat.typepopl != '2minnomi':
    #    listnamefeat = []
    #    for name in listnamefeattemp:
    #        if not isinstance(dictpopl[name][0], str):
    #            listnamefeat.append(name)
    #else:
    #    listnamefeat = listnamefeattemp
    
    numbfeat = np.empty(gdat.numbpopl, dtype=int)
    for k in gdat.indxpopl:
        listnamefeat[k] = list(dictpopl[gdat.listnamepopl[k]].keys())
        listlablfeat[k], listscalfeat[k] = tdpy.retr_listlablscalpara(listnamefeat[k])
        numbfeat[k] = len(listnamefeat[k])
        gdat.indxfeat[k] = np.arange(numbfeat[k])
        print('gdat.listnamepopl[k]')
        print(gdat.listnamepopl[k])
        print('listnamefeat[k]')
        print(listnamefeat[k])
        print('')
    
    for k in gdat.indxpopl:
        
        namepopl = gdat.listnamepopl[k]

        # check for finiteness
        #indx = np.where(np.isfinite(dictpopl['massstar']))[0]
        #for name in listnamefeattemp:
        #    dictpopl[name] = dictpopl[name][indx]

        #indx = np.where(np.isfinite(dictpopl['duratran']))[0]
        #for name in listnamefeattemp:
        #    dictpopl[name] = dictpopl[name][indx]
        
        # check number of targets for each feature
        numbkeys = len(listnamefeat[k])
        numbtarg = np.empty(numbkeys, dtype=int)
        for n in gdat.indxfeat[k]:
            numbtarg[n] = dictpopl[namepopl][listnamefeat[k][n]].size
        
        if np.unique(numbtarg).size >  1:
            print('listnamefeat')
            print(listnamefeat)
            for n in gdat.indxfeat[k]:
                print('dictpopl[namepopl][listnamefeat[k][n]]')
                summgene(dictpopl[namepopl][listnamefeat[k][n]])
            raise Exception('')
        numbtarg = numbtarg[0]
        
        listsamp = np.empty((numbtarg, numbfeat[k]))
        for n in gdat.indxfeat[k]:
            listsamp[:, n] = dictpopl[namepopl][listnamefeat[k][n]]
            
        # visualize the population
        print('Visualizing the population...')
        boolscat = False
        tdpy.mcmc.plot_grid(gdat.pathimag, namepopl, listsamp, listlablfeat[k], boolscat=boolscat, boolplotpair=True, boolplottria=False, numbbinsplot=40, \
                                                        listscalpara=listscalfeat[k], typefileplot=typefileplot, liststrgvarb=listnamefeat[k])
        
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
    for k in gdat.indxpopl:
        
        continue

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
        
        for k in gdat.indxfeat[k]: 
            gdat.tempfrst = arrytemp[:, k]
            gdat.tempfrststdv = arrystdvtemp[:, k]
            for u in gdat.indxfeat[k]: 
                gdat.tempseco = arrytemp[:, u]
                gdat.tempsecostdv = arrystdvtemp[:, u]
                
                minmxpos = np.amin(gdat.tempfrst) / 1.1
                maxmxpos = np.amax(gdat.tempfrst) * 1.1
                minmypos = np.amin(gdat.tempseco) / 1.1
                maxmypos = np.amax(gdat.tempseco) * 1.1
                
                if not ( \
                        listnamefeat[k] == 'ptmp' and listnamefeat[u] == 'logg' or \
                        listnamefeat[k] == 'tday' and listnamefeat[u] == 'ptmp' or \
                        listnamefeat[k] == 'logg' and listnamefeat[u] == 'albe' or \
                        listnamefeat[k] == 'ptmp' and listnamefeat[u] == 'albe' or \
                        listnamefeat[k] == 'tday' and listnamefeat[u] == 'albe' or \
                        listnamefeat[k] == 'tday' and listnamefeat[u] == 'logg' \
                       ):
                    continue
                
                print('%s, %s:' % (listlabl[k], listlabl[u]))
                
                # calculate PCC
                coef, pval = scipy.stats.pearsonr(gdat.tempfrst, gdat.tempseco)
                listcoef[u, k] = coef
                
                # sample from a linear model
                numbdata = 2 * numbplan
                strgextn = '%d_' % b + listnamepara[k]
                featpost = tdpy.mcmc.samp(gdat, pathimag, numbsampwalk, numbsampburnwalk, numbsampburnwalkseco, retr_llik, \
                                                listlablpara, listscalpara, listminmpara, listmaxmpara, listmeangauspara, liststdvgauspara, \
                                                    numbdata, strgextn=strgextn, strgplotextn=strgplotextn, verbtype=0)
                
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
                path = pathimag + 'scat_%s_%s_%d.%s' % (listnamefeat[k], listnamefeat[u], b, strgplotextn)
                print('Writing to %s...' % path)
                plt.savefig(path)
                plt.close()

    # list of parameter labels and units
    listlablpara = [['$u_1$', ''], ['$u_2$', ''], ['$P$', 'days'], ['$i$', 'deg'], ['$\\rho$', ''], ['$C$', '']]
    # list of parameter scalings
    listscalpara = ['self', 'self', 'self', 'self', 'self', 'self']
    # list of parameter minima
    listminmpara = [-1., -1., 0.2,   0.,  0.,-1e-1]
    # list of parameter maxima
    listmaxmpara = [ 3.,  3., 0.4, 89.9, 0.6, 1e-1]
    
    for numbspottemp in range(gdat.numbspot):
        listlablpara += [['$\\theta_{%d}$' % numbspottemp, 'deg'], ['$\\phi_{%d}$' % numbspottemp, 'deg'], ['$R_{%d}$' % numbspottemp, '']]
        listscalpara += ['self', 'self', 'self']
        listminmpara += [-90.,   0.,  0.]
        listmaxmpara += [ 90., 360., 0.4]
        if gdat.boolevol:
            listlablpara += [['$T_{s;%d}$' % numbspottemp, 'day'], ['$\\sigma_{s;%d}$' % numbspottemp, '']]
            listscalpara += ['self', 'self']
            listminmpara += [gdat.minmtime, 0.1]
            listmaxmpara += [gdat.maxmtime, 20.]
            
    listminmpara = np.array(listminmpara)
    listmaxmpara = np.array(listmaxmpara)
    listmeangauspara = None
    liststdvgauspara = None
    
    # number of parameters
    numbpara = len(listlablpara)
    # number of walkers
    numbwalk = max(20, 2 * numbpara)
        
    numbdata = gdat.lcurdataused.size
    
    # number of degrees of freedom
    gdat.numbdoff = numbdata - numbpara
    
    indxpara = np.arange(numbpara)

    listpost = tdpy.mcmc.samp(gdat, gdat.pathimag, gdat.numbsampwalk, gdat.numbsampburnwalk, gdat.numbsampburnwalkseco, retr_llik, \
            listlablpara, listscalpara, listminmpara, listmaxmpara, listmeangauspara, liststdvgauspara, numbdata, boolpool=True, \
                    #retr_lpri=retr_lpri, \
                    strgextn=gdat.strgextn, samptype='emce')

    # search for clusters

    ## interpolate distribution
    np.interp2d(grid)


