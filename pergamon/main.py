import os
import sys

import numpy as np
#import scipy.stats
#import matplotlib.pyplot as plt

import ephesus
from tdpy import gdatstrt
import tdpy
from tdpy import summgene 
import miletos
import pcat

def retr_modl_corr(gdat, para, inpt):
    
    angl = para[0]
    gamm = para[1]

    slop = -1. / np.tan(angl)
    intc = gamm / np.sin(angl)
    
    line = inpt * slop + intc
    
    return line, []


def retr_llik_corr(para, gdat):
    
    angl = para[0]
    gamm = para[1]

    dist = np.cos(angl) * gdat.tempfrst + np.sin(angl) * gdat.tempseco - gamm
    vari = np.cos(angl)**2 * gdat.tempfrststdv**2 + np.sin(angl)**2 * gdat.tempsecostdv**2

    llik = -0.5 * np.sum(dist**2 / vari)

    return llik


def init( \
        # type of the population
        ## '2minnomi': 2-minute cadence target list of TESS during the nominal mission
        listtypepopl=['2minnomi'], \
        
        typeobjtdete='exopmock', \

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
    
    # correlation search
    numbsampfeww = 1000
    numbsampburnwalk = 1000
    numbsampwalk = 10000
    numbsampburnwalkseco = 2000
    
    ## conversion factors
    factrsrj, factrjre, factrsre, factmsmj, factmjme, factmsme, factaurs = ephesus.retr_factconv()

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
    gdat.numbpopl = len(gdat.listtypepopl)
    gdat.indxpopl = np.arange(gdat.numbpopl)
    for k in gdat.indxpopl:
        
        print('Population: %s' % gdat.listtypepopl[k])
        
        dictindxpopl = dict()
        # get the dictionaries holding the population properties
        if gdat.listtypepopl[k] == 'hjuppcur':
            listlablinpt = ['log $g_{\mathrm{p}}$ [cgs]', r'$T_{\mathrm{eq}}$ [K]', '[Fe/H] [dex]', r'$T_{\mathrm{day}}$ [K]', '$A_g$']
            listnameinpt = ['logg', 'ptmp', 'meta', 'tday', 'albe']
            numbcompinpt = len(listnameinpt)
            indxcompinpt = np.arange(numbcompinpt)
            path = gdat.pathdata + 'catlpcurtess.csv'
            dictpopl['hjuppcur'] = pd.read_csv(path)

        if gdat.listtypepopl[k] == '2minnomi':
            dictpopl['2minnomi'] = miletos.retr_dictcatltic8(gdat.listtypepopl[k])
            dictpopl['2minnomi']['nois'] = ephesus.retr_noistess(dictpopl['2minnomi']['tmag'])

        if gdat.listtypepopl[k] == 'toyystar':
            listradistar = tdpy.icdf_powr(np.random.rand(numbsamp), 0.1, 10., -2.)
            listmassstar = tdpy.icdf_powr(np.random.rand(numbsamp), 0.1, 10., -2.)
        
        if typeobjtdete == 'exopexar':
            dictpopl['exopexar'] = retr_dictexar()
            dictpopl['exopexar']['nois'] = ephesus.retr_noistess(dictpopl['exopexar']['vmagsyst'])
        
        if typeobjtdete == 'exopmock':
        
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
            if typeobjtdete == 'slen':
                dictpopl['exopmock']['masscomp'] = tdpy.util.icdf_powr(np.random.rand(numbtarg), 5., 200., 2.)
                dictpopl['exopmock']['masstotl'] = dictpopl['exopmock']['masscomp'] + dictpopl['exopmock']['massstar']
            if typeobjtdete == 'exopmock':
                dictpopl['exopmock']['massplan'] = tdpy.util.icdf_powr(np.random.rand(numbtarg), 5., 200., 2.)
                dictpopl['exopmock']['masstotl'] = dictpopl['exopmock']['massplan'] / factmsme + dictpopl['exopmock']['massstar']
            dictpopl['exopmock']['smax'] = ephesus.retr_smaxkepl(dictpopl['exopmock']['peri'], dictpopl['exopmock']['masstotl'])
            dictpopl['exopmock']['rsma'] = dictpopl['exopmock']['radistar'] / dictpopl['exopmock']['smax'] / factaurs
            
            dictpopl['exopmocktran'] = dict()
            indx = np.where(dictpopl['exopmock']['rsma'] > dictpopl['exopmock']['cosi'])[0]
            for name in dictpopl['exopmock'].keys():
                dictpopl['exopmocktran'][name] = dictpopl['exopmock'][name][indx]
            
            dictpopl['exopmocktran']['duratran'] = ephesus.retr_duratran(dictpopl['exopmock']['peri'][indx], \
                                                                               dictpopl['exopmock']['rsma'][indx], \
                                                                               dictpopl['exopmock']['cosi'][indx])
            if typeobjtdete == 'slen':
                dictpopl['amplslen'] = retr_amplslen(dictpopl['peri'], dictpopl['radistar'], dictpopl['masscomp'], dictpopl['massstar'])
                dictpopl['s2no'] = np.sqrt(dictpopl['duratran']) * dictpopl['amplslen'] / dictpopl['nois']
            if typeobjtdete == 'exopmock':
                dictpopl['exopmocktran']['rrat'] = dictpopl['exopmocktran']['radiplan'] / dictpopl['exopmocktran']['radistar'] / factrsre
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
            
    gdat.listnamepoplfinl = list(dictpopl.keys())
    gdat.numbpoplfinl = len(gdat.listnamepoplfinl)
    gdat.indxpoplfinl = np.arange(gdat.numbpoplfinl)

    # make sure fields with strings are not included
    #if gdat.listtypepopl[k] != '2minnomi':
    #    listnametotl = []
    #    for name in listnametotltemp:
    #        if not isinstance(dictpopl[name][0], str):
    #            listnametotl.append(name)
    #else:
    #    listnametotl = listnametotltemp
        
    for k in gdat.indxpoplfinl:
        
        namepopl = gdat.listnamepoplfinl[k]
        print(namepopl)
        listnametotl = list(dictpopl[namepopl].keys())

        listlablpara, listscalpara = tdpy.retr_listlablscalpara(listnametotl)
        
        numbnametotl = len(listnametotl)
        indxnametotl = np.arange(numbnametotl)
        
        # check for finiteness
        #indx = np.where(np.isfinite(dictpopl['massstar']))[0]
        #for name in listnametotltemp:
        #    dictpopl[name] = dictpopl[name][indx]

        #indx = np.where(np.isfinite(dictpopl['duratran']))[0]
        #for name in listnametotltemp:
        #    dictpopl[name] = dictpopl[name][indx]
        
        # check number of targets for each feature
        numbkeys = len(listnametotl)
        numbtarg = np.empty(numbkeys, dtype=int)
        for n in range(numbkeys):
            numbtarg[n] = dictpopl[namepopl][listnametotl[n]].size
        if np.unique(numbtarg).size != 1:
            print('listnametotl')
            print(listnametotl)
            print('numbtarg')
            print(numbtarg)
            raise Exception('')
        numbtarg = numbtarg[0]
        
        listsamp = np.empty((numbtarg, numbnametotl))
        for n in indxnametotl:
            listsamp[:, n] = dictpopl[namepopl][listnametotl[n]]
            
            print(listnametotl[n])
            summgene(listsamp[:, n], boolfull=True)
            print('')

        # visualize the population
        print('Visualizing the population...')
        boolscat = False
        tdpy.mcmc.plot_grid(gdat.pathimag, gdat.listnamepoplfinl[k], listsamp, listlablpara, boolscat=boolscat, \
                                                        listscalpara=listscalpara, typefileplot=typefileplot)#, join=True)
        
        print('')

    for k in gdat.indxpoplfinl:
        
        # search for correlations
        print('Searching for correlations...')    

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
            
        # derived quantities
        # semi-major axis
        smax = radistar / rsta
        smaxstdv = np.sqrt((radistarstdv / radistar)**2 + (rstastdv / rsta)**2)
        # radius of the planet
        radiplan = rrat * radistar
        radiplanstdv = np.sqrt((rratstdv / rrat)**2 + (radistarstdv / radistar)**2)
        # plenatery surface gravity
        grav = 26.18 * massplan / radiplan**2
        stdvgrav = grav * np.sqrt((massplanstdv / massplan) **2 + (2. * radiplanstdv / radiplan)**2)
        # log of the plenatery surface gravity
        logg = np.log10(grav)
        loggstdv = 0.434 * stdvgrav / grav
        # equilibrium temperature of the planet
        tmpteqiu = tmptstar * np.sqrt(1. / 2. / smax)
        tmpteqiustdv = tmpteqiu * np.sqrt((tmptstarstdv / tmptstar)**2 + (0.5 * smaxstdv / smax)**2) / np.sqrt(2.)
        
        # load variables into the array
        arryinpttemp = np.empty((numbplan, numbcompinpt))
        arryinptstdvtemp = np.empty((numbplan, numbcompinpt))
        # logg 
        arryinpttemp[:, 0] = logg
        arryinptstdvtemp[:, 0] = loggstdv
        # plenatery equilibrium temperature
        arryinpttemp[:, 1] = tmpteqiu
        arryinptstdvtemp[:, 1] = tmpteqiustdv
        # stellar metallicity
        arryinpttemp[:, 2] = metastar
        arryinptstdvtemp[:, 2] = metastarstdv
        # dayside temperature
        arryinpttemp[:, 3] = tday
        arryinptstdvtemp[:, 3] = tdaystdv
        # dayside albedo
        arryinpttemp[:, 4] = albe
        arryinptstdvtemp[:, 4] = albestdv
        
        listlablpara = [[r'$\alpha$', ''], [r'$\rho$', '']]
        listscalpara = ['self', 'self']
        listmeangauspara = None
        liststdvgauspara = None
        listminmpara = np.array([0, -1e1])
        listmaxmpara = np.array([np.pi, 1e1])
        numbpara = len(listscalpara)
        
        numbwalk = max(20, 2 * numbpara)
        numbsamp = numbwalk * (numbsampwalk - numbsampburnwalkseco)
        indxsamp = np.arange(numbsamp)
        
        listcoef = np.zeros((numbcompinpt, numbcompinpt))
        
        for k in indxcompinpt: 
            gdat.tempfrst = arryinpttemp[:, k]
            gdat.tempfrststdv = arryinptstdvtemp[:, k]
            for u in indxcompinpt: 
                gdat.tempseco = arryinpttemp[:, u]
                gdat.tempsecostdv = arryinptstdvtemp[:, u]
                
                minmxpos = np.amin(gdat.tempfrst) / 1.1
                maxmxpos = np.amax(gdat.tempfrst) * 1.1
                minmypos = np.amin(gdat.tempseco) / 1.1
                maxmypos = np.amax(gdat.tempseco) * 1.1
                
                if not ( \
                        listnameinpt[k] == 'ptmp' and listnameinpt[u] == 'logg' or \
                        listnameinpt[k] == 'tday' and listnameinpt[u] == 'ptmp' or \
                        listnameinpt[k] == 'logg' and listnameinpt[u] == 'albe' or \
                        listnameinpt[k] == 'ptmp' and listnameinpt[u] == 'albe' or \
                        listnameinpt[k] == 'tday' and listnameinpt[u] == 'albe' or \
                        listnameinpt[k] == 'tday' and listnameinpt[u] == 'logg' \
                       ):
                    continue
                
                print('%s, %s:' % (listlablinpt[k], listlablinpt[u]))
                
                # calculate PCC
                coef, pval = scipy.stats.pearsonr(gdat.tempfrst, gdat.tempseco)
                listcoef[u, k] = coef
                
                # sample from a linear model
                numbdata = 2 * numbplan
                strgextn = '%d_' % b + listnameinpt[k]
                parapost = tdpy.mcmc.samp(gdat, pathimag, numbsampwalk, numbsampburnwalk, numbsampburnwalkseco, retr_llik, \
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
                
                postslop = -1. / np.tan(parapost[:, 0])
                medislop = np.median(postslop)
                lowrslop = np.median(postslop) - np.percentile(postslop, 16)
                upprslop = np.percentile(postslop, 84) - np.median(postslop)
                titl = 'PCC = %.3g, Slope: %.3g $\substack{+%.2g \\\\ -%.2g}$' % (listcoef[u, k], medislop, upprslop, lowrslop)
                axis.set_xlabel(listlablinpt[k])
                axis.set_ylabel(listlablinpt[u])
        
                plt.tight_layout()
                path = pathimag + 'scat_%s_%s_%d.%s' % (listnameinpt[k], listnameinpt[u], b, strgplotextn)
                print('Writing to %s...' % path)
                plt.savefig(path)
                plt.close()

        # search for clusters

        ## interpolate distribution
        interp2d(grid)

        ## selection effects
        miletos.retr_precreca(grid)
        
        ## occurence rate
        occu = dens / reca * prec


