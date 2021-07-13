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
        # type of the population
        ## '2minnomi': 2-minute cadence target list of TESS during the nominal mission
        typepopl=['2minnomi'], \
        
        listpopl=['exopmock'], \

        # type of measurement
        ## 'tess': type of measurement
        listtypemeas=['tess'], \
        # type of the data
        typedata = 'mock', \
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
    
    strgplotextn = 'pdf'

    # correlation search
    numbsampfeww = 1000
    numbsampburnwalk = 1000
    numbsampwalk = 10000
    numbsampburnwalkseco = 2000
    
    print('typepopl')
    print(typepopl)
    print('listpopl')
    print(listpopl)

    ## conversion factors
    dictfact = ephesus.retr_factconv()

    # plotting
    
    if gdat.typedata == 'mock':
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
    
    dictpopl = dict()
    gdat.numbpopl = len(gdat.typepopl)
    gdat.indxpopl = np.arange(gdat.numbpopl)
        
    if typepopl == 'fstr':
        
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
            if k == 0:
                listnamefstr = linesplt
                print('len(listnamefstr)')
                print(len(listnamefstr))
                #['peri', 'dura', 'dept', 'dcyctran', 'dcycinge', 'SN', 'SPN', 'OOTmag', 'sigamp']
                for namefstr in listnamefstr:
                    dictpopl['fstr'][namefstr] = []
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
                print('numbfiel')
                print(numbfiel)
                if numbfiel == 31:
                    ticifstr.append(float(linesplt[0]))
                    for n in range(1, numbfiel):
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

    # get the dictionaries holding the population properties
    if gdat.typepopl == 'hjuppcur':
        listlabl = ['log $g_{\mathrm{p}}$ [cgs]', r'$T_{\mathrm{eq}}$ [K]', '[Fe/H] [dex]', r'$T_{\mathrm{day}}$ [K]', '$A_g$']
        listname = ['logg', 'ptmp', 'meta', 'tday', 'albe']
        path = gdat.pathdata + 'catlpcurtess.csv'
        dictpopl['hjuppcur'] = pd.read_csv(path)

    if gdat.typepopl == 'ffimm135nomi' or gdat.typepopl == '2minnomi':
        dictpopl[gdat.typepopl] = miletos.retr_dictcatltic8(gdat.typepopl)
        dictpopl[gdat.typepopl]['nois'] = ephesus.retr_noistess(dictpopl[gdat.typepopl]['tmag'])

    if gdat.typepopl == 'toyystar':
        listradistar = tdpy.icdf_powr(np.random.rand(numbsamp), 0.1, 10., -2.)
        listmassstar = tdpy.icdf_powr(np.random.rand(numbsamp), 0.1, 10., -2.)
        
    listlablfeat = [[] for k in gdat.indxpopl]
    listscalfeat = [[] for k in gdat.indxpopl]
    listnamefeat = [[] for k in gdat.indxpopl]
    indxfeat = [[] for k in gdat.indxpopl]
    for k in gdat.indxpopl:
        if listpopl == 'exopexar':
            dictpopl['exopexar'] = retr_dictexar()
            dictpopl['exopexar']['nois'] = ephesus.retr_noistess(dictpopl['exopexar']['vmagsyst'])
        
        if listpopl == 'fstrpcan':
            indxfullpcan = np.where(dispfull == 'PC')[0]
            dictpopl['tici'] = ticifull[indxfullpcan]
        
        if listpopl == 'exopmock':
        
            dictpopl['2minnomi']['probexop'] = pcat.icdf_gaustrun(np.random.rand(dictpopl['2minnomi']['radistar'].size), 0.02, 0.002, 0, np.inf)
            dictpopl['2minnomi']['boolexop'] = np.random.rand(dictpopl['2minnomi']['radistar'].size) < dictpopl['2minnomi']['probexop']
            
            dictpopl['exopmock'] = dict()
            indx = np.where(dictpopl['2minnomi']['boolexop'])[0]
            for name in dictpopl['2minnomi'].keys():
                dictpopl['exopmock'][name] = dictpopl['2minnomi'][name][indx]

            numbtarg = dictpopl['exopmock']['radistar'].size
            
            dictpopl['exopmock']['incl'] = np.random.rand(numbtarg) * 90.
            dictpopl['exopmock']['cosi'] = np.cos(np.pi / 180. * dictpopl['exopmock']['incl'])
            dictpopl['exopmock']['peri'] = tdpy.util.icdf_powr(np.random.rand(numbtarg), 0.3, 20., 2.)
            dictpopl['exopmock']['radiplan'] = tdpy.util.icdf_powr(np.random.rand(numbtarg), 1., 23., 2.)
            if listpopl == 'slen':
                dictpopl['exopmock']['massfeat'] = tdpy.util.icdf_powr(np.random.rand(numbtarg), 5., 200., 2.)
                dictpopl['exopmock']['masstotl'] = dictpopl['exopmock']['massfeat'] + dictpopl['exopmock']['massstar']
            if listpopl == 'exopmock':
                dictpopl['exopmock']['massplan'] = tdpy.util.icdf_powr(np.random.rand(numbtarg), 5., 200., 2.)
                dictpopl['exopmock']['masstotl'] = dictpopl['exopmock']['massplan'] / dictfact['msme'] + dictpopl['exopmock']['massstar']
            dictpopl['exopmock']['smax'] = ephesus.retr_smaxkepl(dictpopl['exopmock']['peri'], dictpopl['exopmock']['masstotl'])
            dictpopl['exopmock']['rsma'] = dictpopl['exopmock']['radistar'] / dictpopl['exopmock']['smax'] / dictfact['aurs']
            
            dictpopl['exopmocktran'] = dict()
            indx = np.where(dictpopl['exopmock']['rsma'] > dictpopl['exopmock']['cosi'])[0]
            for name in dictpopl['exopmock'].keys():
                dictpopl['exopmocktran'][name] = dictpopl['exopmock'][name][indx]
            
            dictpopl['exopmocktran']['duratran'] = ephesus.retr_duratran(dictpopl['exopmock']['peri'][indx], \
                                                                               dictpopl['exopmock']['rsma'][indx], \
                                                                               dictpopl['exopmock']['cosi'][indx])
            if listpopl == 'slen':
                dictpopl['amplslen'] = retr_amplslen(dictpopl['peri'], dictpopl['radistar'], dictpopl['massfeat'], dictpopl['massstar'])
                dictpopl['s2nr'] = np.sqrt(dictpopl['duratran']) * dictpopl['amplslen'] / dictpopl['nois']
            if listpopl == 'exopmock':
                dictpopl['exopmocktran']['rrat'] = dictpopl['exopmocktran']['radiplan'] / dictpopl['exopmocktran']['radistar'] / dictfact['rsre']
                dictpopl['exopmocktran']['dept'] = 1e6 * dictpopl['exopmocktran']['rrat']**2 # [ppm]
                
                print('temp')
                dictpopl['exopmocktran']['numbtsec'] = np.ones_like(dictpopl['exopmocktran']['peri'])
                
                dictpopl['exopmocktran']['numbtran'] = 27.3 * dictpopl['exopmocktran']['numbtsec'] / dictpopl['exopmocktran']['peri']
                dictpopl['exopmocktran']['s2nrblss'] = np.sqrt(dictpopl['exopmocktran']['numbtran'] * dictpopl['exopmocktran']['duratran']) * \
                                                                    dictpopl['exopmocktran']['dept'] / dictpopl['exopmocktran']['nois']
                
                dictpopl['exopmocktran']['booldete'] = dictpopl['exopmocktran']['s2nrblss'] > 7
        
                # add the population of detections
                #dictpopl['exoptrandete'] = dict()
                #indx = np.random.choice(np.arange(dictpopl['2minnomi']['radistar'].size), size=2500)
                #for name in dictpopl['2minnomi'].keys():
                #    dictpopl['exoptrandete'][name] = dictpopl['2minnomi'][name][indx]
            
                # add the population of detections
                dictpopl['exoptranmockdete'] = dict()
                indx = np.where(dictpopl['exopmocktran']['booldete'])[0]
                for name in dictpopl['exopmocktran'].keys():
                    dictpopl['exoptranmockdete'][name] = dictpopl['exopmocktran'][name][indx]
                
    gdat.listnamepopl = list(dictpopl.keys())
    print('gdat.listnamepopl')
    print(gdat.listnamepopl)

    gdat.numbpopl = len(gdat.listnamepopl)
    gdat.indxpopl = np.arange(gdat.numbpopl)

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
        indxfeat[k] = np.arange(numbfeat[k])
    print('listnamefeat')
    print(listnamefeat)
    

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
        for n in indxfeat[k]:
            numbtarg[n] = dictpopl[namepopl][listnamefeat[k][n]].size
        
        if np.unique(numbtarg).size >  1:
            print('listnamefeat')
            print(listnamefeat)
            for n in indxfeat[k]:
                print('dictpopl[namepopl][listnamefeat[k][n]]')
                summgene(dictpopl[namepopl][listnamefeat[k][n]])
            raise Exception('')
        numbtarg = numbtarg[0]
        
        listsamp = np.empty((numbtarg, numbfeat[k]))
        for n in indxfeat[k]:
            listsamp[:, n] = dictpopl[namepopl][listnamefeat[k][n]]
            
            print(listnamefeat[k][n])
            summgene(listsamp[:, n], boolfull=True)
            print('')

        # visualize the population
        print('Visualizing the population...')
        boolscat = False
        tdpy.mcmc.plot_grid(gdat.pathimag, namepopl, listsamp, listlablfeat[k], boolscat=boolscat, \
                                                        listscalpara=listscalfeat[k], typefileplot=typefileplot)#, join=True)
        
        print('')
    
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
        
        for k in indxfeat[k]: 
            gdat.tempfrst = arrytemp[:, k]
            gdat.tempfrststdv = arrystdvtemp[:, k]
            for u in indxfeat[k]: 
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

    # search for clusters

    ## interpolate distribution
    np.interp2d(grid)

    ## selection effects
    miletos.retr_precreca(grid)
    
    ## occurence rate
    occu = dens / reca * prec


