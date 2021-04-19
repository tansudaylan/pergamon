import os
import sys

import numpy as np
import matplotlib.pyplot as plt

import ephesus
import tdpy
from tdpy import summgene 
import miletos

def init():
   
    typedata = 'mock'

    # paths
    pathbase = os.environ['TDGU_DATA_PATH'] + '/lsst_tran/'
    pathdata = pathbase + 'data/'
    pathimag = pathbase + 'imag/'
    os.system('mkdir -p %s' % pathdata)
    os.system('mkdir -p %s' % pathimag)
    
    if typedata:
        # time stamps
        listtime = np.random.rand(1000)
    
        ## number of years into the mission
        numbyear = 10
        numbsamp = int(numbyear * 100)
        ## photometric precision
        stdvphot = 1e-3
    
    ephesus.expl_popl(typepopl='exoptran')
    
    # run BLS
    rsma = 0.1
    cosi = 0.
    dept = 0.01
    dura = ephesus.retr_dura(peri, rsma, cosi)
    dcyc = dura / peri
    meansampitra = dcyc * numbsamp
    meansampotra = numbsamp - meansampitra
    stdvitra = stdvphot / np.sqrt(meansampotra - 1)
    stdvotra = stdvphot / np.sqrt(meansampitra - 1)
    stdv = np.sqrt(stdvitra**2 + stdvotra**2)
    s2nr = dept / stdv
    
    # selection effect

    plotfiletype = 'pdf'
    
    figr, axis = plt.subplots(figsize=(6, 3))
    axis.plot(peri, s2nr)
    axis.set_xlabel('Period [day]')
    axis.set_xscale('log')
    axis.set_yscale('log')
    plt.tight_layout()
    path = pathimag + 'peris2nr.%s' % plotfiletype
    plt.savefig(path)
    plt.close()

from tdpy.util import summgene
import tdpy.util
import tdpy.mcmc
import ephesus.util

import pandas as pd

import numpy as np
import scipy.stats

import matplotlib.pyplot as plt

import os

def retr_modl(gdat, para, inpt):
    
    angl = para[0]
    gamm = para[1]

    slop = -1. / np.tan(angl)
    intc = gamm / np.sin(angl)
    
    line = inpt * slop + intc
    
    return line, []


def retr_llik(gdat, para):
    
    angl = para[0]
    gamm = para[1]

    dist = np.cos(angl) * gdat.tempfrst + np.sin(angl) * gdat.tempseco - gamm
    vari = np.cos(angl)**2 * gdat.tempfrststdv**2 + np.sin(angl)**2 * gdat.tempsecostdv**2

    llik = -0.5 * np.sum(dist**2 / vari)

    return llik

# paths
pathbase = os.environ['TDGU_DATA_PATH'] + '/pcurcorr/'
pathdata = pathbase + 'data/'
pathimag = pathbase + 'imag/'

# plotting
strgplotextn = 'pdf'

# construct global object
gdat = tdpy.util.gdatstrt()

#listdatainpt = np.array([\
#                    ['WASP-5 b'   ,   0., 0.16, 1996., 97.],\
#                    ['WASP-36 b'  , 0.15, 0.15, 1420., 165.],\
#                    ['WASP-43 b'  , 0.14, 0.06, 1632., 48.],\
#                    ['WASP-46 b'  , 0.27, 0.30, 1830., 140.],\
#                    ['WASP-64 b'  , 0.39, 0.28, 1915., 99.],\
#                    ['WASP-77 A b', 0.06, 0.05, 1830., 36.],\
#                    ['WASP-78 b'  , 0.19, 0.25, 2490., 185.],\
#                    ['WASP-100 b' , 0.26, 0.09, 2315., 78.],\
#                    ['WASP-19 b'  , 0.20, 0.08, 2469., 55.],\
#                    ['WASP-121 b' , 0.32, 0.08, 2469., 55.],\
#                    ['WASP-18 b'  , 0.  , 0.02, 3068., 57.],\
#                    ])

# read data
liststrgplan = [[], []]
liststrgtarg = ['TESS_targs', 'other_targs']
data = [[], []]
for a, strgtarg in enumerate(liststrgtarg):
    path = pathdata + strgtarg + '.txt'
    print('Reading from %s...' % path)
    objtfile = open(path, 'r')
    
    for k, line in enumerate(objtfile):
        
        if k < 1:
            #print(line)
            continue
        linesplt = line.split('\t')
        linesplt = [linesplttemp for linesplttemp in linesplt if linesplttemp != '']
        liststrgplan[a].append(linesplt[0])
        data[a].append(np.array(linesplt[1:]).astype(float))
    data[a] = np.vstack(data[a])

listlablinpt = ['log $g_{\mathrm{p}}$ [cgs]', r'$T_{\mathrm{eq}}$ [K]', '[Fe/H] [dex]', r'$T_{\mathrm{day}}$ [K]', '$A_g$']
listnameinpt = ['logg', 'ptmp', 'meta', 'tday', 'albe']
numbcompinpt = len(listnameinpt)
indxcompinpt = np.arange(numbcompinpt)

for b in range(3):
    
    if b == 0:
        arry = data[0]
        print('TESS only')
    if b == 1:
        arry = data[1]
        print('Others only')
    if b == 2:
        arry = np.concatenate(data, 0)
        print('All')
    
    # read data array    
    numbplan = arry.shape[0]
    # effective temperature of the star
    tmptstar = arry[:, 0]
    tmptstarstdv = arry[:, 1]
    # metallicity of the star
    metastar = arry[:, 2]
    metastarstdv = arry[:, 3]
    # radius of the star
    radistar = arry[:, 4]
    radistarstdv = arry[:, 5]
    # mass of the planet
    massplan = arry[:, 6]
    massplanstdv = arry[:, 7]
    # semi-major axis divided by the radius of the star
    rsta = arry[:, 8]
    rstastdv = arry[:, 9]
    # radius of the planet divided by the radius of the star
    rrat = arry[:, 10]
    rratstdv = arry[:, 11]
    # dayside temperature of the planet
    tday = arry[:, 12]
    tdaystdv = np.maximum(arry[:, 13], arry[:, 14])
    # geometric albedo 
    albe = arry[:, 15]
    albestdv = np.maximum(arry[:, 16], arry[:, 17])
        
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
    
    numbsampfeww = 1000
    numbsampburnwalk = 1000
    numbsampwalk = 10000
    numbsampburnwalkseco = 2000
    
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
            
            print(titl)
            axis.set_xlabel(listlablinpt[k])
            axis.set_ylabel(listlablinpt[u])
    
            plt.tight_layout()
            path = pathimag + 'scat_%s_%s_%d.%s' % (listnameinpt[k], listnameinpt[u], b, strgplotextn)
            #print('Writing to %s...' % path)
            plt.savefig(path)
            plt.close()
            print('')


