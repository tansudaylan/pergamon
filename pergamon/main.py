import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import ephesus
import tdpy
import lygos
from tdpy import summgene 
import miletos
import mergen

def retr_modl_corr(gdat, feat, inpt):
    '''
    Linear model
    '''    
    angl = feat[0]
    gamm = feat[1]

    slop = -1. / np.tan(angl)
    intc = gamm / np.sin(angl)
    
    line = inpt * slop + intc
    
    return line, []


def retr_llik_corr(feat, gdat):
    '''
    Likelihood for linear model
    '''
    angl = feat[0]
    gamm = feat[1]

    dist = np.cos(angl) * gdat.tempfrst + np.sin(angl) * gdat.tempseco - gamm
    vari = np.cos(angl)**2 * gdat.tempfrststdv**2 + np.sin(angl)**2 * gdat.tempsecostdv**2

    llik = -0.5 * np.sum(dist**2 / vari)

    return llik


def init( \
        # type of analysis
        typeanls, \
        
        # dictionary of features of the populations
        dictpopl=dict(), \

        # plotting
        typefileplot='pdf', \

        # verbosity level
        typeverb=1, \
        
        # Boolean flag to interpolate the scatter plots
        boolintpscat=False, \

        # Boolean flag to perform PCA decomposition
        boolprca=False, \

        # Boolean flag to search for correlations
        boolsrchcorr=False, \

        # name of the feature holding the labels of the samples
        namelablsamp=None, \
        
        # base path
        pathbase=None, \

        # data path
        pathdata=None, \

        # image path
        pathimag=None, \

        **args, \
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

    # copy named arguments to the global object
    for strg, valu in args.items():
        setattr(gdat, strg, valu)

    if gdat.pathbase is None:
        gdat.pathbase = os.environ['PERGAMON_DATA_PATH'] + '/%s/' % gdat.typeanls
    
    # paths
    if not (gdat.pathbase is not None and gdat.pathimag is None and gdat.pathdata is None ) or \
           (gdat.pathbase is None and gdat.pathimag is not None and gdat.pathdata is None) or \
           (gdat.pathbase is None and gdat.pathimag is None and gdat.pathdata is not None):
        print('gdat.pathbase')
        print(gdat.pathbase)
        print('gdat.pathdata')
        print(gdat.pathdata)
        print('gdat.pathimag')
        print(gdat.pathimag)
        raise Exception('')
    
    if gdat.pathdata is None or gdat.pathimag is None:
        gdat.pathdata = gdat.pathbase + 'data/'
        gdat.pathimag = gdat.pathbase + 'imag/'
    os.system('mkdir -p %s' % gdat.pathdata)
    os.system('mkdir -p %s' % gdat.pathimag)
    
    # plotting
    gdat.strgplotextn = 'pdf'

    ## conversion factors
    dictfact = ephesus.retr_factconv()
    
    print('gdat.typeanls')
    print(gdat.typeanls)

    # list of population names
    ## 2-minute cadence target list of TESS during the nominal mission
    #if gdat.typeanls == 'exopmock2minnomi':
    #gdat.listnamepopl = ['2minnomi', 'exopmock2minnomi']
    
    ## 'bholmock': mock BH candidate subpopulation of '2minnomi'
    ## '2minnomi': 2-min cadence target list during the nominal mission
    #if gdat.typeanls == 'bholmock2minnomi':
    #gdat.listnamepopl = ['2minnomi', 'bholmock2minnomi']
    ## 'exopmock': mock exoplanet subpopulation of '2minnomi'
    ## 'exoptoii': TOIs
    
    # preliminary setup for analyses
    if gdat.typeanls == 'mocklsst':
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


    if gdat.typeanls == 'featsupntess':
        for a in [1, 2, 3]:
            path = gdat.pathdata + 'Cycle%d-matched.csv' % a
            strgcycl = 'cyc%d' % a
            print('Reading from %s...' % path)
            dictpopl[strgcycl] = pd.read_csv(path).to_dict(orient='list')
            dictpopl[strgcycl]['magtdisc'] = dictpopl[strgcycl].pop('Discovery Mag/Flux')
            listname = list(dictpopl[strgcycl].keys())
            for name in listname:
                dictpopl[strgcycl][name] = np.array(dictpopl[strgcycl][name])
                if not (dictpopl[strgcycl][name].dtype == np.float64 or dictpopl[strgcycl][name].dtype == np.int64):
                    del dictpopl[strgcycl][name]
                else:
                    dictpopl[strgcycl][''.join(c for c in name if c not in ' ./:' )] = dictpopl[strgcycl].pop(name)
            del dictpopl[strgcycl]['Unnamed0']
        listlistnameplotcomm = ['cyc1', 'cyc2', 'cyc3']

    if gdat.typeanls == 'featexarmassradi' or gdat.typeanls == 'featexaratmo' or gdat.typeanls == 'featexarweakmass' or gdat.typeanls == 'featexartran':
        # properties of confirmed exoplanets
        dictpopl['totl'] = ephesus.retr_dictexar()
        
    if gdat.typeanls == 'featexarweakmass':
        # properties of TOIs
        dictexof = miletos.retr_dictexof()
        indx = np.where((dictexof['nameplan'] == 'TOI 712.01') | (dictexof['nameplan'] == 'TOI 712.02') | (dictexof['nameplan'] == 'TOI 712.03'))[0]
        print(dictexof['nameplan'])
        print(indx)
        for namefeat in dictpopl['totl'].keys():
            if namefeat in ['strgprovmass', 'cosi', 'strgrefrradiplan', 'strgrefrmassplan', 'tagestar', 'stdvtagestar']:
                arry = np.zeros_like(indx) + np.nan
            else:
                arry = dictexof[namefeat][indx]
            print(namefeat)
            dictpopl['totl'][namefeat] = np.concatenate((dictpopl['totl'][namefeat], arry))

    if gdat.typeanls == 'featexarmassradi':
        # create a subpopulation with planets with good measured radii and masses
        indx = []
        for n  in range(dictpopl['totl']['strgrefrmassplan'].size):
            if not ('Calculated Value' in dictpopl['totl']['strgrefrmassplan'][n] or \
                    'Calculated Value' in dictpopl['totl']['strgrefrradiplan'][n]):
                indx.append(n)
        indxmeas = np.array(indx)
        indxgood = np.where(dictpopl['totl']['stdvmassplan'] / dictpopl['totl']['stdvmassplan'] > 5.)[0]
        indx = np.setdiff1d(indxmeas, indxgood)
        dictpopl['meas'] = dict()
        for namefeat in dictpopl['totl'].keys():
            dictpopl['meas'][namefeat] = dictpopl['totl'][namefeat][indx]

    # name of the feature that holds the labels for samples
    if gdat.typeanls == 'featexarmassradi' or gdat.typeanls == 'featexaratmo' or gdat.typeanls == 'feattoiiatmo' or gdat.typeanls == 'featexarweakmass' or gdat.typeanls == 'featexartran':
        namelablsamp = 'nameplan'
            
    if gdat.typeanls == 'featexaratmo' or gdat.typeanls == 'featexarweakmass' or gdat.typeanls == 'featexartran':
        # create a subpopulation with transiting planets
        indx = np.where(dictpopl['totl']['booltran'])[0]
        dictpopl['tran'] = dict()
        for namefeat in dictpopl['totl'].keys():
            print('namefeat')
            print(namefeat)
            dictpopl['tran'][namefeat] = dictpopl['totl'][namefeat][indx]
        
    if gdat.typeanls == 'featexaratmo' or gdat.typeanls == 'featexarweakmass':
        del dictpopl['totl']

        # create a subpopulation with targets discovered by TESS
        indx = np.where(dictpopl['tran']['facidisc'] == 'Transiting Exoplanet Survey Satellite (TESS)')[0]
        dictpopl['tess'] = dict()
        for namefeat in dictpopl['tran'].keys():
            dictpopl['tess'][namefeat] = dictpopl['tran'][namefeat][indx]
        
        # create a subpopulation with targets brighter than 12th TESS magnitude
        indx = np.where(dictpopl['tess']['tmagsyst'] < 12.)[0]
        dictpopl['brgthost'] = dict()
        for namefeat in dictpopl['tess'].keys():
            dictpopl['brgthost'][namefeat] = dictpopl['tess'][namefeat][indx]
        
        # create a subpopulation with planets smaller than four times the Earth
        indx = np.where((dictpopl['brgthost']['radiplan'] < 4.) & (dictpopl['brgthost']['radiplan'] > 2.))[0]
        dictpopl['smalplan'] = dict()
        for namefeat in dictpopl['brgthost'].keys():
            dictpopl['smalplan'][namefeat] = dictpopl['brgthost'][namefeat][indx]
        
        # create a subpopulation with planets with good visibility from LCO (dec = -30 deg)
        indx = np.where(dictpopl['smalplan']['declstar'] < 30.)[0]
        dictpopl['lcoo'] = dict()
        for namefeat in dictpopl['smalplan'].keys():
            dictpopl['lcoo'][namefeat] = dictpopl['smalplan'][namefeat][indx]
        
        # create a subpopulation with planets with good atmospheric characterization potential
        indx = np.where((dictpopl['lcoo']['esmm'] > 5.) | (dictpopl['lcoo']['tsmm'] > 50.))[0]
        dictpopl['goodatmo'] = dict()
        for namefeat in dictpopl['lcoo'].keys():
            dictpopl['goodatmo'][namefeat] = dictpopl['lcoo'][namefeat][indx]
       
    if gdat.typeanls == 'featexaratmo':
        # create a subpopulation with good mass measurements
        indx = np.where(dictpopl['goodatmo']['massplan'] / dictpopl['goodatmo']['stdvmassplan'] > 5.)[0]
        dictpopl['goodmass'] = dict()
        for namefeat in dictpopl['goodatmo'].keys():
            dictpopl['goodmass'][namefeat] = dictpopl['goodatmo'][namefeat][indx]
        
    if gdat.typeanls == 'featexarweakmass':
        # create a subpopulation with weak mass measurements
        indx = np.where((dictpopl['goodatmo']['massplan'] / dictpopl['goodatmo']['stdvmassplan'] < 5.) | ~np.isfinite(dictpopl['goodatmo']['stdvmassplan']))[0]
        dictpopl['weakmass'] = dict()
        for namefeat in dictpopl['goodatmo'].keys():
            dictpopl['weakmass'][namefeat] = dictpopl['goodatmo'][namefeat][indx]
        
    if gdat.typeanls == 'featexarmassradi':
        
        minmradi = np.nanmin(dictpopl['meas']['radiplan'])
        maxmradi = np.nanmax(dictpopl['meas']['radiplan'])
        binsradi = np.linspace(minmradi, 24., 15)
        meanradi = (binsradi[1:] + binsradi[:-1]) / 2.
        arry = np.empty((meanradi.size, 5))
        arry[:, 0] = meanradi
        for l in range(meanradi.size):
            indx = np.where((dictpopl['meas']['radiplan'] > binsradi[l]) & (dictpopl['meas']['radiplan'] < binsradi[l+1]) & \
                                                                                        (dictpopl['meas']['massplan'] / dictpopl['meas']['stdvmassplan'] > 5.))[0]
            arry[l, 1] = np.nanmedian(dictpopl['meas']['massplan'][indx])
            arry[l, 2] = np.nanstd(dictpopl['meas']['massplan'][indx])
            arry[l, 3] = np.nanmedian(dictpopl['meas']['densplan'][indx])
            arry[l, 4] = np.nanstd(dictpopl['meas']['densplan'][indx])
        
        #meanmass = 
        #indx = np.where(meanradi < 1.23)[0]
        #np.log10(meanradi)
        #mass = 
        
        path = os.environ['MILETOS_DATA_PATH'] + '/data/massfromradi.csv'
        print('Writing to %s...' % path)
        np.savetxt(path, arry, fmt='%8.4g')

    if gdat.typeanls == 'feattoiiatmo' or gdat.typeanls == 'feattoiifstr':
        # properties of TOIs
        dictpopl['totl'] = miletos.retr_dictexof()

    # Tom's list
    #listnameplan = [ \
    #                # ESM list
    #                'AU Mic c', 'TOI-216.02', 'HD 63433 c', 'TOI-837 b', 'DS Tuc A b', 'TOI-2076 b', 'HR 858 d', 'LHS 3844 b', 'GJ 1252 b', 'HR 858 b', 'HR 858 c', \
    #                'TOI-1807 b', 'TOI-561 b', \
    #                # TSM list
    #                'AU Mic c', 'TOI-178 g', 'HD 63433 c', 'LP 791-18 c', 'TOI-178 d', 'L 98-59 b', 'TOI-540 b', 'LTT 1445 A b', 'LP 791-18 b', 'TOI-2076 b', \
    #                'HD 191939 b', 'TOI-1130 b', 'HD 63433 b', 'HR 858 d', 'LHS 3844 b', 'GJ 1252 b' 
    #               ]
                
    if gdat.typeanls == 'feattoiiatmo' or gdat.typeanls == 'feattoiifstr':
        
        # turn zero TSM ACWG or ESM ACWG into NaN
        indx = np.where(dictpopl['totl']['tsmmacwg'] == 0)[0]
        dictpopl['totl']['tsmmacwg'][indx] = np.nan
        indx = np.where(dictpopl['totl']['esmmacwg'] == 0)[0]
        dictpopl['totl']['esmmacwg'][indx] = np.nan

        # total number of TOIs
        numbtoii = len(dictpopl['totl']['strgcomm'])
        indxtoii = np.arange(numbtoii)

        ## PC
        dictpopl['pcan'] = dict()
        indxpcan = []
        for n in indxtoii:
            if not dictpopl['totl']['boolfpos'][n]:
                indxpcan.append(n)
        indxpcan = np.array(indxpcan, dtype=int)
        for name in dictpopl['totl'].keys():
            dictpopl['pcan'][name] = dictpopl['totl'][name][indxpcan]
        
        ## FP
        dictpopl['fpos'] = dict()
        indxfpos = np.setdiff1d(indxtoii, indxpcan)
        for name in dictpopl['totl'].keys():
            dictpopl['fpos'][name] = dictpopl['totl'][name][indxfpos]
        
        
    if gdat.typeanls == 'feattoiifstr':
        
        # replace BJD with TESS-truncated BJD (TBJD)
        dictpopl['totl']['epoctess'] = dictpopl['totl']['epoc'] - 2457000
        
        # delete unnecessary features
        listkeys = list(dictpopl['totl'].keys())
        for name in listkeys:
            if name.startswith('stdv'):
                del dictpopl['totl'][name]
            if name in ['tici', 'epoc', 'hmagsyst', 'kmagsyst', 'jmagsyst']:
                del dictpopl['totl'][name]
        
        ## FaintStars
        dictpopl['fstr'] = dict()
        indxfstr = []
        for n in indxtoii:
            if isinstance(dictpopl['totl']['strgcomm'][n], str) and 'faint-star' in dictpopl['totl']['strgcomm'][n]:
                indxfstr.append(n)
        indxfstr = np.array(indxfstr, dtype=int)
        for name in dictpopl['totl'].keys():
            dictpopl['fstr'][name] = dictpopl['totl'][name][indxfstr]
        
        print('indxpcan')
        summgene(indxpcan)
        print('indxfpos')
        summgene(indxfpos)
        print('indxfstr')
        summgene(indxfstr)

    # conversion factors
    gdat.dictfact = ephesus.retr_factconv()

    # settings
    ## plotting
    gdat.numbcyclcolrplot = 300
    gdat.alphraww = 0.2
    ### percentile for zoom plots of relative flux
    gdat.pctlrflx = 95.
    gdat.typefileplot = 'pdf'
    gdat.figrsize = [6, 4]
    gdat.figrsizeydob = [8., 4.]
    gdat.figrsizeydobskin = [8., 2.5]
        
    if gdat.typeanls == 'feattceefstr':
        
        # base path for FaintStars
        pathfstr = '/Users/tdaylan/Documents/work/data/external/FaintStars/'
        
        # astronet output
        dictpopl['fstrastr'] = dict()
        # FaintStars Fs + Ps
        dictpopl['fstrhvet'] = dict()
        # FaintStars Fs
        dictpopl['fstrfpos'] = dict()
        # FaintStars PCs
        dictpopl['fstrpcan'] = dict()
        
        ## TESS sectors
        listtsecnomi = range(1, 27)
        #listtsecnomi = range(1, 3)
        numbtsecnomi = len(listtsecnomi)
        indxtsecnomi = np.arange(numbtsecnomi)
        
        # header
        listnamevarbfstr = 'tic per per_err epo epo_err rprs rprs_err b b_err ars ars_err rprs_odd rprs_err_odd rprs_even rprs_err_even dep dur din u1 '
        listnamevarbfstr += 'u2 SN SPN OOTmag sig_amp sig_pri sig_sec sig_ter sig_pos sig_oe dmm shape asym sig_fa1 sig_fa2 fred phs_pri phs_sec phs_ter phs_pos '
        listnamevarbfstr += 'dep_sec deperr_sec sig_oe_alt sig_12 sig_23 sig_13'
        listnamevarbfstr = listnamevarbfstr.split(' ')
        for n in range(len(listnamevarbfstr)):
            listnamevarbfstr[n] = ''.join(listnamevarbfstr[n].split('_'))
        print('listnamevarbfstr')
        print(listnamevarbfstr)
        
        listdictpopl = dict()
        listdictpopl['fstrastr'] = [dict() for o in indxtsecnomi]
        for namevarbfstr in listnamevarbfstr:
            for o in indxtsecnomi:
                listdictpopl['fstrastr'][o][namevarbfstr] = []
            dictpopl['fstrastr'][namevarbfstr] = []
            dictpopl['fstrpcan'][namevarbfstr] = []
        
        print('Reading metrics...')
        for o in indxtsecnomi:
            path = pathfstr + 'metr/metr_sc%02d.txt' % listtsecnomi[o]
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
                    listdictpopl['fstrastr'][o][listnamevarbfstr[n]].append(valu)

                # check
                if numbfiel != 45:
                    raise Exception('')
            for namevarbfstr in listnamevarbfstr:
                listdictpopl['fstrastr'][o][namevarbfstr] = np.array(listdictpopl['fstrastr'][o][namevarbfstr]) 
                #if listdictpopl['fstrastr'][o][namevarbfstr].size == 0:
                #    raise Exception('')
            print('%d metrics...' % listdictpopl['fstrastr'][o]['tic'].size)
        print('')
        
        # read TIC IDs and dispositions
        print('Reading dispositions...')
        pathbase = pathfstr + 'combined/'
        ticihvettsec = [[] for o in indxtsecnomi]
        disphvettsec = [[] for o in indxtsecnomi]
        for o in indxtsecnomi:
            path = pathbase + 'sector-%d.tier2.csv' % listtsecnomi[o]
            print('Reading from %s...' % path)
            dictquer = pd.read_csv(path).to_dict(orient='list')
            ticihvettsec[o] = np.array(dictquer['TIC'])
            disphvettsec[o] = dictquer['Final']
            print('%d dispositions...' % ticihvettsec[o].size)
        print('')
        
        print('Filtering out those targets for which there is a disposition, but not metric...')
        ticihvettsectemp = [[] for o in indxtsecnomi]
        for o in indxtsecnomi:
            for ticihvet in ticihvettsec[o]:
                if ticihvet in listdictpopl['fstrastr'][o]['tic']:
                    ticihvettsectemp[o].append(ticihvet)
            ticihvettsectemp[o] = np.array(ticihvettsectemp[o])
        for o in indxtsecnomi:
            print('Sector %2d: %4d of %4d dispositions have metrics...' % (listtsecnomi[o], ticihvettsectemp[o].size, ticihvettsec[o].size))
        ticihvettsec = ticihvettsectemp

        print('Merging lists of TIC IDs for targets with metrics across nominal mission...')
        ticiastrtsec = [[] for o in indxtsecnomi]
        for o in indxtsecnomi:
            ticiastrtsec[o] = listdictpopl['fstrastr'][o]['tic']
        ticiastrconc = np.concatenate(ticiastrtsec)
        print('Total number of metrics: %d' % ticiastrconc.size)
        ticiastruniq, indxuniq, indxinve, numbtici = np.unique(ticiastrconc, return_index=True, return_inverse=True, return_counts=True)
        print('Total number of targets: %d' % ticiastruniq.size)

        print('Merging lists of TIC IDs for targets with dispositions across nominal mission...')
        ticihvetconc = np.concatenate(ticihvettsec)
        print('Total number of dispositions: %d' % ticihvetconc.size)
        ticihvetuniq, indxuniq, indxinve, numbtici = np.unique(ticihvetconc, return_index=True, return_inverse=True, return_counts=True)
        print('Total number of targets with dispositions: %d' % ticihvetuniq.size)

        for namevarbfstr in listnamevarbfstr:
            dictpopl['fstrastr'][namevarbfstr] = []
            dictpopl['fstrhvet'][namevarbfstr] = []
            dictpopl['fstrfpos'][namevarbfstr] = []
            dictpopl['fstrpcan'][namevarbfstr] = []
        
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
                for o in indxtsecnomi:
                    indx = np.where(tici == ticitsec[o])[0]
                    if indx.size > 0:
                        if a == 1:
                            dispthis.append(disptsec[o][indx[0]])
                        indxtsecthis.append(o)
                
                # index of the last sector of this target
                indxtseclast = indxtsecthis[-1]

                ## find the metric index of the target in the last sector
                indx = np.where(listdictpopl['fstrastr'][indxtseclast]['tic'] == tici)[0]
                
                if indx.size == 0:
                    raise Exception('')

                ## collect metrics of the target in the last sector
                if a == 0:
                    ### all human-vetted targets
                    for namevarbfstr in listnamevarbfstr:
                        dictpopl['fstrastr'][namevarbfstr].append(listdictpopl['fstrastr'][indxtseclast][namevarbfstr][indx[0]])
                if a == 1:
                    ### all human-vetted targets
                    for namevarbfstr in listnamevarbfstr:
                        dictpopl['fstrhvet'][namevarbfstr].append(listdictpopl['fstrastr'][indxtseclast][namevarbfstr][indx[0]])
                    ### Ps
                    if dispthis[-1] == 'P':
                        for namevarbfstr in listnamevarbfstr:
                            dictpopl['fstrpcan'][namevarbfstr].append(listdictpopl['fstrastr'][indxtseclast][namevarbfstr][indx[0]])
                    ### Fs
                    else:
                        for namevarbfstr in listnamevarbfstr:
                            dictpopl['fstrfpos'][namevarbfstr].append(listdictpopl['fstrastr'][indxtseclast][namevarbfstr][indx[0]])
         
            
        for namevarbfstr in listnamevarbfstr:
            dictpopl['fstrastr'][namevarbfstr] = np.array(dictpopl['fstrastr'][namevarbfstr])
            if len(dictpopl['fstrpcan'][namevarbfstr]) == 0:
                raise Exception('')
            dictpopl['fstrfpos'][namevarbfstr] = np.array(dictpopl['fstrfpos'][namevarbfstr])
            dictpopl['fstrpcan'][namevarbfstr] = np.array(dictpopl['fstrpcan'][namevarbfstr])
            dictpopl['fstrhvet'][namevarbfstr] = np.array(dictpopl['fstrhvet'][namevarbfstr])
        
        # crossmatch with TIC to get additional features
        #dictquer = miletos.xmat_tici(ticifstrastr)

        #numbastr = len(ticifstrastr)
        #for namevarbfstr in listnamevarbfstr:
        #    dictpopl['fstrastr'][namevarbfstr] = np.empty(numbastr)
        #for k, tici in enumerate(ticifstrastr):
        #    indx = np.where(ticifstr == tici)[0]
        #    if indx.size > 0:
        #        for namevarbfstr in listnamevarbfstr:
        #            dictpopl['fstrastr'][namevarbfstr][k] = dictpopl['fstrastr'][namevarbfstr][indx[0]]

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
        # get the dictionaries holding the population properties
        dictpopl['hjuppcur'] = pd.read_csv(path)

    if gdat.typeanls == 'exoptoyy':
        listradistar = tdpy.icdf_powr(np.random.rand(numbsamp), 0.1, 10., -2.)
        listmassstar = tdpy.icdf_powr(np.random.rand(numbsamp), 0.1, 10., -2.)
        
    if gdat.typeanls == 'exopexar':
        dictpopl['exopexar']['nois'] = ephesus.retr_noistess(dictpopl['exopexar']['vmagsyst'])
        


    #if gdat.typeanls == 'tessnomi2minexopmocktoyy':
    if gdat.typeanls[12:].startswith('bholmock') or gdat.typeanls[12:].startswith('exopmock'):
        
        # name of the star population
        namepoplstar = gdat.typeanls[:12]
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
        #dictpopl[namepoplstar]['proboccu'] = tdpy.icdf_gaustrun(np.random.rand(dictpopl[namepoplstar]['radistar'].size), 0.02, 0.002, 0, np.inf)
        
        # Boolean flag of occurence
        #dictpopl[namepoplstar]['booloccu'] = np.random.rand(dictpopl[namepoplstar]['radistar'].size) < dictpopl[namepoplstar]['proboccu']
        
        # stars with occurence
        namepoploccu = gdat.typeanls[:20]
        
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
        if gdat.typeanls[12:].startswith('bholmock'):
            dictpopl[namepoploccu]['massbhol'] = tdpy.util.icdf_powr(np.random.rand(numbtargoccu), minmmassbhol, maxmmassbhol, 2.)
            dictpopl[namepoploccu]['masstotl'] = dictpopl[namepoploccu]['massbhol'] + dictpopl[namepoploccu]['massstar']
        if gdat.typeanls[12:].startswith('exopmock'):
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
        if gdat.typeanls[12:].startswith('exopmock'):
            dictpopl[namepoploccutran]['rrat'] = dictpopl[namepoploccutran]['radiplan'] / dictpopl[namepoploccutran]['radistar'] / dictfact['rsre']
            dictpopl[namepoploccutran]['dept'] = 1e3 * dictpopl[namepoploccutran]['rrat']**2 # [ppt]
        if gdat.typeanls[12:].startswith('bholmock'):
            dictpopl[namepoploccutran]['amplslen'] = ephesus.retr_amplslen(dictpopl[namepoploccutran]['peri'], dictpopl[namepoploccutran]['radistar'], \
                                                                                dictpopl[namepoploccutran]['massbhol'], dictpopl[namepoploccutran]['massstar'])
        
        # detection
        print('temp -- assuming all targets have one sector')
        dictpopl[namepoploccutran]['numbtsec'] = np.ones(numbtargoccutran)
        dictpopl[namepoploccutran]['numbtran'] = 27.3 * dictpopl[namepoploccutran]['numbtsec'] / dictpopl[namepoploccutran]['peri']
        ## SNR
        if gdat.typeanls[12:].startswith('exopmock'):
            dictpopl[namepoploccutran]['sdee'] = np.sqrt(dictpopl[namepoploccutran]['duratran']) * dictpopl[namepoploccutran]['dept'] / dictpopl[namepoploccutran]['nois']
        if gdat.typeanls[12:].startswith('bholmock'):
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
        

    gdat.listnamepopl = np.array(list(dictpopl.keys()))
    print('gdat.listnamepopl')
    print(gdat.listnamepopl)

    gdat.numbpopl = len(gdat.listnamepopl)
    if gdat.numbpopl == 0:
        raise Exception('')

    gdat.indxpopl = np.arange(gdat.numbpopl)

    listlablfeat = [[] for k in gdat.indxpopl]
    listlablfeattotl = [[] for k in gdat.indxpopl]
    listscalfeat = [[] for k in gdat.indxpopl]
    listnamefeat = [[] for k in gdat.indxpopl]
    gdat.indxfeat = [[] for k in gdat.indxpopl]
    
    numbfeat = np.empty(gdat.numbpopl, dtype=int)
    for k in gdat.indxpopl:
        listnamefeat[k] = list(dictpopl[gdat.listnamepopl[k]].keys())
        listlablfeat[k], listscalfeat[k] = tdpy.retr_listlablscalpara(listnamefeat[k])
        numbfeat[k] = len(listnamefeat[k])
        gdat.indxfeat[k] = np.arange(numbfeat[k])
        
        listlablfeattotl[k] = tdpy.retr_labltotl(listlablfeat[k])

    for k in gdat.indxpopl:
        
        namepopl = gdat.listnamepopl[k]

        print('namepopl')
        print(namepopl)

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
            if numbtarg[n] == 0:
                raise Exception('')
        if np.unique(numbtarg).size >  1:
            print('listnamefeat')
            print(listnamefeat)
            for n in gdat.indxfeat[k]:
                print('dictpopl[namepopl][listnamefeat[k][n]]')
                summgene(dictpopl[namepopl][listnamefeat[k][n]])
            raise Exception('')
        
    # filter some features of the population
    listsampfilt = [[] for k in gdat.indxpopl]
    listlablfeatfilt = [[] for k in gdat.indxpopl]
    listlablfeattotlfilt = [[] for k in gdat.indxpopl]
    listnamefeatfilt = [[] for k in gdat.indxpopl]
    listscalfeatfilt = [[] for k in gdat.indxpopl]
    dictpoplfilt = dict()

    for k in gdat.indxpopl:
        dictpoplfilt[gdat.listnamepopl[k]] = dict()
        for n in gdat.indxfeat[k]:
            
            if gdat.typeanls == 'feattoiifstr' and not listnamefeat[k][n] in ['radiplan', 'toii']:
                continue
        
            if gdat.typeanls == 'featexartran' and not listnamefeat[k][n] in ['dcyc', 'peri', 'nameplan']:
                continue
        
            if gdat.typeanls == 'featexarmassradi' and not listnamefeat[k][n] in ['radiplan', 'massplan', 'tmptplan', 'radistar', 'densplan']:
                continue
        
            if (gdat.typeanls == 'feattoiiatmo' or gdat.typeanls == 'featexaratmo' or gdat.typeanls == 'featexarweakmass') and \
                        not listnamefeat[k][n] in [ \
                                                   'radiplan', 'massplan', 'tmptplan', 'radistar', 'stdvmassplan', 'peri', 'duratran', \
                                                   'tsmm', 'esmm', 'vmagsyst', 'densplan', 'nameplan', 'loecstar', 'laecstar', \
                                                   'rascstar', 'declstar', 'esmmacwg', 'tsmmacwg', \
                                                  ]:
                continue
            
            dictpoplfilt[gdat.listnamepopl[k]][listnamefeat[k][n]] = dictpopl[gdat.listnamepopl[k]][listnamefeat[k][n]]
            
            listlablfeattotlfilt[k].append(listlablfeattotl[k][n])
            
            # exclude features with string value
            if listnamefeat[k][n] in ['strgcomm', 'namesyst', 'nameplan', 'facidisc', 'strgprovmass', 'strgrefrradiplan', 'strgrefrmassplan']:
                continue
            
            # exclude features with unuseful IDs
            if listnamefeat[k][n] in ['tici']:
                continue
            
            listsampfilt[k].append(np.array(dictpopl[gdat.listnamepopl[k]][listnamefeat[k][n]]).astype(float))
            listlablfeatfilt[k].append(listlablfeat[k][n])
            listnamefeatfilt[k].append(listnamefeat[k][n])
            listscalfeatfilt[k].append(listscalfeat[k][n])
        listsampfilt[k] = np.vstack(listsampfilt[k]).T
    listsamp = listsampfilt
    listlablfeat = listlablfeatfilt
    listlabltotlfeat = listlablfeattotlfilt
    listnamefeat = listnamefeatfilt
    listscalfeat = listscalfeatfilt
    
    for k in gdat.indxpopl:
        numbfeat[k] = len(listnamefeat[k])
        gdat.indxfeat[k] = np.arange(numbfeat[k])

    if gdat.typeanls == 'featexarweakmass':
        
        import astropy
        from astropy.visualization import astropy_mpl_style, quantity_support
        import astropy.units as u
        
        # location object for LCO
        objtlocalcoo = astropy.coordinates.EarthLocation(lat=-29.01418*u.deg, lon=-70.69239*u.deg, height=2515.819*u.m)
        
        # time object for the year
        objttimeyear = astropy.time.Time(astropy.time.Time('2022-01-11 00:00:00').jd + np.linspace(0., 365., 10000), format='jd', location=objtlocalcoo)
        timeside = objttimeyear.sidereal_time('mean')
        
        # alt-az coordinate object for the Sun
        objtcoorsunnalazyear = astropy.coordinates.get_sun(objttimeyear)
            
        # delt time arry for night
        timedelt = np.linspace(-12., 12., 1000)
        
        for n in range(len(dictpoplfilt['weakmass'][namelablsamp])):
            
            continue

            #if n == 1:
            #    raise Exception('')

            strgsamp = ''.join(dictpoplfilt['weakmass'][namelablsamp][n].split(' '))
            
            # choose the time sample where the local sidereal time is closest to the right ascension
            indx = np.argmin(abs(180. - abs(objtcoorsunnalazyear.ra.degree - dictpopl['weakmass']['rascstar'][n])))# + abs(timeside.degree - dictpopl['goodatmo']['rascstar'][n]))
            
            # time object for night at midnight
            objttimenighcent = astropy.time.Time(int(objttimeyear[indx].jd), format='jd', location=objtlocalcoo) - 12 * u.hour
            objttimenigh = objttimenighcent + timedelt * u.hour
            
            # frame object for LCO at night
            objtframlcoonigh = astropy.coordinates.AltAz(obstime=objttimenigh, location=objtlocalcoo)
            
            # alt-az coordinate object for the Sun
            objtcoorsunnalaznigh = astropy.coordinates.get_sun(objttimenigh).transform_to(objtframlcoonigh)
            # alt-az coordinate object for the Moon
            objtcoormoonalaznigh = astropy.coordinates.get_moon(objttimenigh).transform_to(objtframlcoonigh)
            
            # alt-az coordinate object for the planet
            objtcoorplanalaznigh = astropy.coordinates.SkyCoord(ra=dictpoplfilt['weakmass']['rascstar'][n], \
                                                            dec=dictpoplfilt['weakmass']['declstar'][n], frame='icrs', unit='deg').transform_to(objtframlcoonigh)
            
            strgtitl = '%s, %s UTC' % (dictpoplfilt['weakmass'][namelablsamp][n], objttimenighcent.iso)

            # plot air mass
            figr, axis = plt.subplots(figsize=(4, 4))
            massairr = objtcoorplanalaznigh.secz
            indx = np.where(np.isfinite(massairr) & (massairr > 0))[0]
            plt.plot(timedelt[indx], massairr[indx])
            axis.fill_between(timedelt, 0, 90, objtcoorsunnalaznigh.alt < -0*u.deg, color='0.5', zorder=0)
            axis.fill_between(timedelt, 0, 90, objtcoorsunnalaznigh.alt < -18*u.deg, color='k', zorder=0)
            axis.fill_between(timedelt, 0, 90, (massairr > 2.) | (massairr < 1.), color='r', alpha=0.3, zorder=0)
            axis.set_xlabel('$\Delta t$ [hour]')
            axis.set_ylabel('Airmass')
            limtxdat = [np.amin(timedelt), np.amax(timedelt)]
            axis.set_title(strgtitl)
            axis.set_xlim(limtxdat)
            axis.set_ylim([1., 2.])
            path = gdat.pathimag + 'airmass_%s.%s' % (strgsamp, gdat.strgplotextn)
            print('Writing to %s...' % path)
            plt.savefig(path)
            
            # plot altitude
            figr, axis = plt.subplots(figsize=(4, 4))
            axis.plot(timedelt, objtcoorsunnalaznigh.alt, color='orange', label='Sun')
            axis.plot(timedelt, objtcoormoonalaznigh.alt, color='gray', label='Moon')
            axis.plot(timedelt, objtcoorplanalaznigh.alt, color='blue', label=dictpoplfilt['weakmass'][namelablsamp][n])
            axis.fill_between(timedelt, 0, 90, objtcoorsunnalaznigh.alt < -0*u.deg, color='0.5', zorder=0)
            axis.fill_between(timedelt, 0, 90, objtcoorsunnalaznigh.alt < -18*u.deg, color='k', zorder=0)
            axis.fill_between(timedelt, 0, 90, (massairr > 2.) | (massairr < 1.), color='r', alpha=0.3, zorder=0)
            axis.legend(loc='upper left')
            plt.ylim([0, 90])
            axis.set_title(strgtitl)
            axis.set_xlim(limtxdat)
            axis.set_xlabel('Hours from EDT Midnight')
            axis.set_ylabel('Altitude [deg]')
            
            path = gdat.pathimag + 'altitude_%s.%s' % (strgsamp, gdat.strgplotextn)
            print('Writing to %s...' % path)
            plt.savefig(path)


    # store features on disc
    for k in gdat.indxpopl:
        path = gdat.pathdata + 'dict_%s.csv' % gdat.listnamepopl[k] 
        print('Writing to %s...' % path)
        pd.DataFrame.from_dict(dictpoplfilt[gdat.listnamepopl[k]]).to_csv(path, header=listlablfeattotlfilt[k], index=False, float_format='%.8g')
    
    # visualize the population
    for k in gdat.indxpopl:
        print('Visualizing the population...')
        boolscat = True
        if namelablsamp is None:
            listlablsamp = None
        else:
            listlablsamp = dictpoplfilt[gdat.listnamepopl[k]][namelablsamp]
        tdpy.plot_grid(gdat.pathimag, gdat.listnamepopl[k], listsamp[k], listlablfeat[k], boolscat=boolscat, boolplotindi=True, boolplotpair=True, boolplottria=False, \
                                                                listscalpara=listscalfeat[k], typefileplot=typefileplot, liststrgvarb=listnamefeat[k], listlablsamp=listlablsamp)
    
    # plot populations together
    numblistnameplotcomm = len(listlistnameplotcomm)
    indxlistnameplotcomm = np.arange(indxlistnameplotcomm)
    for u in indxlistnameplotcomm:
        
        indx = []
        for nameplotcomm in listnameplotcomm[u]:
            indx.append(np.where())
        indx = np.array(indx)
        for k in gdat.indxpopl:
            tdpy.plot_grid(gdat.pathimag, gdat.listnamepopl[k], listsamp[k], listlablfeat[k], boolscat=boolscat, boolplotindi=True, boolplotpair=True, boolplottria=False, \
                                                                listscalpara=listscalfeat[k], typefileplot=typefileplot, liststrgvarb=listnamefeat[k])
        
    
    # print out population summary
    indxpoplcomp = None
    if gdat.typeanls == 'toiifstr':
        indxpoplcomp = 1
    if gdat.typeanls == 'feattoiiatmo':
        indxpoplcomp = 0
    if gdat.typeanls == 'tesspcanfstr':
        indxpoplcomp = 3
    
    errr = np.empty((4, gdat.numbpopl))
    sigm = np.empty(gdat.numbpopl - 1)
    for n in gdat.indxfeat[0]:
        print(listnamefeat[0][n])
        for k in gdat.indxpopl:
            strgvarb = gdat.listnamepopl[k]
            if gdat.numbpopl == 1 or k == indxpoplcomp or indxpoplcomp is None:
                varbcomp = None
            else:
                varbcomp = dictpopl[gdat.listnamepopl[indxpoplcomp]][listnamefeat[indxpoplcomp][n]]
            summgene(dictpopl[gdat.listnamepopl[k]][listnamefeat[k][n]], boolslin=True, strgvarb=strgvarb, varbcomp=varbcomp)
    
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
                gdat.tempfrst = dictpopl[gdat.listnamepopl[k]][listnamefeat[k][n]]
                gdat.tempfrststdv = arrystdvtemp[:, n]
                for m in gdat.indxfeat[k]: 
                    gdat.tempseco = dictpopl[gdat.listnamepopl[k]][listnamefeat[k][m]]
                    gdat.tempsecostdv = arrystdvtemp[:, m]
                    
                    minmxpos = np.amin(gdat.tempfrst) / 1.1
                    maxmxpos = np.amax(gdat.tempfrst) * 1.1
                    minmypos = np.amin(gdat.tempseco) / 1.1
                    maxmypos = np.amax(gdat.tempseco) * 1.1
                    
                    print('%s, %s:' % (listlabl[n], listlabl[m]))
                    
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
        mergen.exec_prca_skit(listvarb)

    # search for clusters

    if gdat.boolintpscat:
        ## interpolate distribution
        np.interp2d(grid)


