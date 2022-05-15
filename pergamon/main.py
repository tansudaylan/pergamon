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
        
        # list of dictionaries of labels and colors for compared populations (one for each comparison)
        listdictlablcolrpopl=None, \

        # Boolean flag to sort the populations according to their sizes
        boolsortpoplsize=True, \

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

    if gdat.listdictlablcolrpopl is None:
        gdat.listdictlablcolrpopl = []

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

    if gdat.pathdata is None:
        gdat.pathbase = os.environ['PERGAMON_DATA_PATH'] + '/%s/' % gdat.typeanls
        gdat.pathdata = gdat.pathbase + 'data/'
        gdat.pathimag = gdat.pathbase + 'imag/'
    os.system('mkdir -p %s' % gdat.pathdata)
    os.system('mkdir -p %s' % gdat.pathimag)
    
    # plotting
    gdat.strgplotextn = 'pdf'

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
        gdat.dictpopl['totl'] = pd.read_csv(path, skiprows=7).to_dict(orient='list')
        for namefeat in gdat.dictpopl['totl'].keys():
            gdat.dictpopl['totl'][namefeat] = np.array(gdat.dictpopl['totl'][namefeat])
            print('namefeat')
            print(namefeat)
            summgene(gdat.dictpopl['totl'][namefeat])


    if gdat.typeanls == 'supntess':
        for a in [1, 2, 3]:
            path = gdat.pathdata + 'Cycle%d-matched.csv' % a
            strgcycl = 'cyc%d' % a
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
            
            # delete unnecessary features
            del gdat.dictpopl[strgcycl]['Unnamed0']
        
        namelablsamp = 'name'

    # get population features
    if gdat.typeanls == 'exar':
        # features of confirmed exoplanets
        gdat.dictpopl['totl'] = ephesus.retr_dictexar()
        
        for uu in range(len(gdat.dictpopl['totl']['rascstar'])):
            print('%s: %g' % (gdat.dictpopl['totl']['nameplan'][uu], gdat.dictpopl['totl']['declstar'][uu]))
        print('')
        print('')
        print('')

    if gdat.typeanls == 'exarweakmass':
        # add TOI-712 into the population
        ## features of TOIs
        dicttoii = ephesus.retr_dicttoii()
        indx = np.where((dicttoii['nameplan'] == 'TOI 712.01') | (dicttoii['nameplan'] == 'TOI 712.02') | (dicttoii['nameplan'] == 'TOI 712.03'))[0]
        print(dicttoii['nameplan'])
        print(indx)
        for namefeat in gdat.dictpopl['totl'].keys():
            if namefeat in ['strgprovmass', 'cosi', 'strgrefrradiplan', 'strgrefrmassplan', 'tagestar', 'stdvtagestar']:
                arry = np.zeros_like(indx) + np.nan
            else:
                arry = dicttoii[namefeat][indx]
            print(namefeat)
            gdat.dictpopl['totl'][namefeat] = np.concatenate((gdat.dictpopl['totl'][namefeat], arry))
    
    gdat.dictindxsamp = dict()
    gdat.dictnumbsamp = dict()
    gdat.dictindxsamp['totl'] = dict()
    if gdat.typeanls.startswith('cosc') or gdat.typeanls.startswith('psys'):
        
        gdat.dictpopl['comp'] = dict()
        gdat.dictindxsamp['comp'] = dict()
    
        #featcosc_s2nr_tess2min_tessnomi2min
        ## 'cosc': features of COSCs
        ## 's2nr': based on signal-to-noise estimation
        ## 'tess2min': 2-min cadence TESS data
        ## 'tessnomi2min': 2-min targets in the nominal mission of TESS
        
        # name of the star population
        strgsplt = gdat.typeanls.split('_')
        # type of estimation
        typeesti = strgsplt[1]
        # type of instrument, cadence, and temporal baseline
        typeinst = strgsplt[2]
        # type of target star population
        typepoplstar = strgsplt[3]
        strgpoplstar = 'star' + typepoplstar
        gdat.dictindxsamp[strgpoplstar] = dict()
        # type of companion
        typecomp = gdat.typeanls[4:8]
        print('typeinst')
        print(typeinst)
        
        # setup for experiment
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
    
        # get a dictionary with features of stars and their companions
        dictpoplstarcomp, _, _ = ephesus.retr_dictpoplstarcomp(typepoplstar, typecomp)
        for name in dictpoplstarcomp.keys():
            gdat.dictpopl[name] = dictpoplstarcomp[name]

        # calculate photometric precision for the star population
        if typeinst.startswith('tess'):
            listkeys = list(gdat.dictpopl.keys())
            print('listkeys')
            print(listkeys)
            dictpopl[strgpoplstar]['nois'] = ephesus.retr_noistess(dictpopl[strgpoplstar]['tmag'])
        elif typeinst.startswith('lsst'):
            dictpopl[strgpoplstar]['nois'] = ephesus.retr_noislsst(dictpopl[strgpoplstar]['rmag'])
        
        # detection
        if typeinst.startswith('lsst'):
            numbvisi = 1000
            gdat.dictpopl['comptran']['sdee'] = gdat.dictpopl['comptran']['dept'] / 5. / gdat.dictpopl['comptran']['nois'] * \
                                                                                            np.sqrt(gdat.dictpopl['comptran']['dcyc'] * numbvisi)
        if typeinst.startswith('tess'):
            print('temp -- assuming all targets have one sector')
            gdat.dictpopl['comptran']['numbtsec'] = np.ones(gdat.dictnumbsamp['comptran'])
            gdat.dictpopl['comptran']['numbtran'] = 27.3 * gdat.dictpopl['comptran']['numbtsec'] / gdat.dictpopl['comptran']['peri']
            ## SNR
            if gdat.typeanls.startswith('psys'):
                gdat.dictpopl['comptran']['sdee'] = np.sqrt(gdat.dictpopl['comptran']['duratran']) * gdat.dictpopl['comptran']['dept'] / gdat.dictpopl['comptran']['nois']
            if gdat.typeanls.startswith('cosc'):
                gdat.dictpopl['comptran']['sdee'] = np.sqrt(gdat.dictpopl['comptran']['duratran']) * gdat.dictpopl['comptran']['amplslen'] / gdat.dictpopl['comptran']['nois']
        
        # detected
        indx = gdat.dictpopl['comptran']['sdee'] > 7
        ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'comptran', 'comptrandete', indx)
        
    if gdat.typeanls == 'toii':
        # features of TOIs
        gdat.dictpopl['totl'] = ephesus.retr_dicttoii()

    if gdat.typeanls == 'hosttoii':
        # features of hosts of TOIs
        gdat.dictpopl['totl'] = ephesus.retr_dicthostplan('toii')
    
    if gdat.typeanls == 'hostexar':
        # features of hosts of exoplanets on NASA Exoplanet Archive
        gdat.dictpopl['totl'] = ephesus.retr_dicthostplan('exar')

    # Tom's list
    #listnameplan = [ \
    #                # ESM list
    #                'AU Mic c', 'TOI-216.02', 'HD 63433 c', 'TOI-837 b', 'DS Tuc A b', 'TOI-2076 b', 'HR 858 d', 'LHS 3844 b', 'GJ 1252 b', 'HR 858 b', 'HR 858 c', \
    #                'TOI-1807 b', 'TOI-561 b', \
    #                # TSM list
    #                'AU Mic c', 'TOI-178 g', 'HD 63433 c', 'LP 791-18 c', 'TOI-178 d', 'L 98-59 b', 'TOI-540 b', 'LTT 1445 A b', 'LP 791-18 b', 'TOI-2076 b', \
    #                'HD 191939 b', 'TOI-1130 b', 'HD 63433 b', 'HR 858 d', 'LHS 3844 b', 'GJ 1252 b' 
    #               ]
                
    # subpopulations
    if gdat.typeanls == 'exar':
        
        if gdat.typeanls == 'exar' or gdat.typeanls == 'exarhabi':
            # planets with habitable insolation
            indx = np.where((gdat.dictpopl['totl']['inso'] > 0.7) & (gdat.dictpopl['totl']['inso'] < 1.1))[0]
            ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'totl', 'gmasinsohabi', indx)
        
        if gdat.typeanls == 'exar' or gdat.typeanls == 'exarmassradi':
            # planets with good measured radii and masses
            indx = []
            for n  in range(gdat.dictpopl['totl']['strgrefrmassplan'].size):
                if not ('Calculated Value' in gdat.dictpopl['totl']['strgrefrmassplan'][n] or \
                        'Calculated Value' in gdat.dictpopl['totl']['strgrefrradiplan'][n]):
                    indx.append(n)
            indxmeas = np.array(indx)
            indxgood = np.where(gdat.dictpopl['totl']['stdvmassplan'] / gdat.dictpopl['totl']['stdvmassplan'] > 5.)[0]
            
            indx = np.setdiff1d(indxmeas, indxgood)
            ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'totl', 'gmas', indx)
            
            # young
            indx = np.where(gdat.dictpopl['totl']['tagestar'] < 100.)[0]
            ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'totl', 'yong', indx)

            # low insolation
            indx = np.where(gdat.dictpopl['totl']['inso'] < 10.)[0]
            ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'totl', 'gmasinsohabi', indx)

            # medium insolation
            indx = np.where((gdat.dictpopl['totl']['inso'] > 10.) & (gdat.dictpopl['totl']['inso'] < 1000.))[0]
            ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'totl', 'gmasinsomedi', indx)

            # high insolation
            indx = np.where(gdat.dictpopl['totl']['inso'] > 1000.)[0]
            ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'totl', 'gmasinsohigh', indx)

        if gdat.typeanls == 'exar':
            # subpopulation with transiting planets
            indx = np.where(gdat.dictpopl['totl']['booltran'])[0]
            ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'totl', 'tran', indx)
            
            # create a subpopulation with targets discovered by TESS
            indx = np.where(gdat.dictpopl['tran']['facidisc'] == 'Transiting Exoplanet Survey Satellite (TESS)')[0]
            ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'tran', 'trantess', indx)
            
            ## TESS-confirmed planets in multi systems
            indx = np.where(gdat.dictpopl['trantess']['numbplantranstar'] > 1)[0]
            ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'trantess', 'trantessmult', indx)
        
            if False and gdat.typeanls == 'exar':

                # create a subpopulation with targets brighter than 12th TESS magnitude
                indx = np.where(gdat.dictpopl['trantess']['tmagsyst'] < 12.)[0]
                ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'trantess', 'trantessbrgt', indx)
                
                # create a subpopulation with planets smaller than four times the Earth
                indx = np.where((gdat.dictpopl['trantessbrgt']['radiplan'] < 4.) & (gdat.dictpopl['trantessbrgt']['radiplan'] > 2.))[0]
                ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'trantessbrgt', 'trantessbrgtsmal', indx)
                
                # create a subpopulation with planets with good visibility from LCO (dec = -30 deg)
                indx = np.where(gdat.dictpopl['trantessbrgtsmal']['declstar'] < 30.)[0]
                ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'trantessbrgtsmal', 'trantessbrgtsmallcoo', indx)
                
                # create a subpopulation with planets with good atmospheric characterization potential
                indx = np.where((gdat.dictpopl['trantessbrgtsmallcoo']['esmm'] > 5.) | (gdat.dictpopl['trantessbrgtsmallcoo']['tsmm'] > 50.))[0]
                ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'trantessbrgtsmallcoo', 'trantessbrgtsmallcooatmo', indx)
               
                # create a subpopulation with good mass measurements
                indx = np.where(gdat.dictpopl['trantessbrgtsmallcooatmo']['massplan'] / gdat.dictpopl['trantessbrgtsmallcooatmo']['stdvmassplan'] > 5.)[0]
                ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'trantessbrgtsmallcooatmo', 'trantessbrgtsmallcooatmogmas', indx)
            
                # create a subpopulation with weak mass measurements
                indx = np.where((gdat.dictpopl['trantessbrgtsmallcooatmo']['massplan'] / gdat.dictpopl['trantessbrgtsmallcooatmo']['stdvmassplan'] < 5.) | \
                                                                                              ~np.isfinite(gdat.dictpopl['trantessbrgtsmallcooatmo']['stdvmassplan']))[0]
                ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'trantessbrgtsmallcooatmo', 'trantessbrgtsmallcooatmowmas', indx)
        
    if gdat.typeanls == 'toii':
        # subpopulations of PCs and FPs
        # total number of TOIs
        numbtoii = len(gdat.dictpopl['totl']['strgcomm'])
        indxtoii = np.arange(numbtoii)

        ### PC
        indx = np.where(gdat.dictpopl['totl']['typedisptess'] == 'PC')[0]
        ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'totl', 'pcan', indx)
                
        # CP
        indx = np.where(gdat.dictpopl['totl']['typedisptess'] == 'CP')[0]
        ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'totl', 'conp', indx)
        
        ## EB
        indx = np.where(gdat.dictpopl['totl']['typedisptess'] == 'EB')[0]
        ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'totl', 'fpos', indx)
        
        # indices of samples of subpopulations
        
        ## FaintStars PM
        indx = []
        for n in indxtoii:
            if isinstance(gdat.dictpopl['totl']['strgcomm'][n], str) and 'found in faint-star QLP search' in gdat.dictpopl['totl']['strgcomm'][n] \
                                                                                                                    and int(gdat.dictpopl['totl']['toii'][n]) < 4500:
                indx.append(n)
        indx = np.array(indx)
        ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'totl', 'fstrpmis', indx)
        
        ## FaintStars EM1
        indx = []
        for n in indxtoii:
            if isinstance(gdat.dictpopl['totl']['strgcomm'][n], str) and 'found in faint-star QLP search' in gdat.dictpopl['totl']['strgcomm'][n] \
                                                                                                                    and int(gdat.dictpopl['totl']['toii'][n]) > 4500:
                indx.append(n)
        indx = np.array(indx)
        ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'totl', 'fstremis', indx)
        
        ## Others
        indx = []
        for n in indxtoii:
            if isinstance(gdat.dictpopl['totl']['strgcomm'][n], str) and not 'found in faint-star QLP search' in gdat.dictpopl['totl']['strgcomm'][n]:
                indx.append(n)
        indx = np.array(indx)
        ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'totl', 'othr', indx)
        
    if gdat.typeanls == 'hosttoii' or gdat.typeanls == 'toii':
        
        ## candidate planets in multi systems
        indx = np.where(gdat.dictpopl['pcan']['numbplantranstar'] > 1)[0]
        ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'pcan', 'pcanmult', indx)
        
        ## confirmed planets in multi systems
        indx = np.where(gdat.dictpopl['conp']['numbplantranstar'] > 1)[0]
        ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'conp', 'conpmult', indx)
        
        ## singles
        #indx = np.where(gdat.dictpopl['pcan']['numbplantranstar'] == 1)[0]
        #ephesus.retr_subp(gdat.dictpopl, gdat.dictnumbsamp, gdat.dictindxsamp, 'totl', 'pcansing', indx)
    
    if gdat.typeanls == 'toii' or gdat.typeanls == 'exar':
        del gdat.dictpopl['totl']['esmm']
    
    if gdat.typeanls == 'toii':
        
        # replace BJD with TESS-truncated BJD (TBJD)
        gdat.dictpopl['totl']['epoctess'] = gdat.dictpopl['totl']['epoc'] - 2457000
        
        listkeys = list(gdat.dictpopl['totl'].keys())

        # delete unnecessary features
        for name in listkeys:
            if name.startswith('stdv'):
                del gdat.dictpopl['totl'][name]
            if name in ['tici', 'epoc', 'hmagsyst', 'kmagsyst', 'jmagsyst']:
                del gdat.dictpopl['totl'][name]
        
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
        
    if gdat.typeanls == 'tceefstr':
        
        # base path for FaintStars
        pathfstr = '/Users/tdaylan/Documents/work/data/external/FaintStars/'
        
        # astronet output
        gdat.dictpopl['fstrastr'] = dict()
        # FaintStars Fs + Ps
        gdat.dictpopl['fstrhvet'] = dict()
        # FaintStars Fs
        gdat.dictpopl['fstrfpos'] = dict()
        # FaintStars PCs
        gdat.dictpopl['fstrpcan'] = dict()
        
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
        
        listgdat = tdpy.gdatstrt()
        listgdat.dictpopl = dict()
        listgdat.dictpopl['fstrastr'] = [dict() for o in indxtsecnomi]
        for namevarbfstr in listnamevarbfstr:
            for o in indxtsecnomi:
                listgdat.dictpopl['fstrastr'][o][namevarbfstr] = []
            gdat.dictpopl['fstrastr'][namevarbfstr] = []
            gdat.dictpopl['fstrpcan'][namevarbfstr] = []
        
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
                    listgdat.dictpopl['fstrastr'][o][listnamevarbfstr[n]].append(valu)

                # check
                if numbfiel != 45:
                    raise Exception('')
            for namevarbfstr in listnamevarbfstr:
                listgdat.dictpopl['fstrastr'][o][namevarbfstr] = np.array(listgdat.dictpopl['fstrastr'][o][namevarbfstr]) 
                #if listgdat.dictpopl['fstrastr'][o][namevarbfstr].size == 0:
                #    raise Exception('')
            print('%d metrics...' % listgdat.dictpopl['fstrastr'][o]['tic'].size)
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
                if ticihvet in listgdat.dictpopl['fstrastr'][o]['tic']:
                    ticihvettsectemp[o].append(ticihvet)
            ticihvettsectemp[o] = np.array(ticihvettsectemp[o])
        for o in indxtsecnomi:
            print('Sector %2d: %4d of %4d dispositions have metrics...' % (listtsecnomi[o], ticihvettsectemp[o].size, ticihvettsec[o].size))
        ticihvettsec = ticihvettsectemp

        print('Merging lists of TIC IDs for targets with metrics across nominal mission...')
        ticiastrtsec = [[] for o in indxtsecnomi]
        for o in indxtsecnomi:
            ticiastrtsec[o] = listgdat.dictpopl['fstrastr'][o]['tic']
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
            gdat.dictpopl['fstrastr'][namevarbfstr] = []
            gdat.dictpopl['fstrhvet'][namevarbfstr] = []
            gdat.dictpopl['fstrfpos'][namevarbfstr] = []
            gdat.dictpopl['fstrpcan'][namevarbfstr] = []
        
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
                indx = np.where(listgdat.dictpopl['fstrastr'][indxtseclast]['tic'] == tici)[0]
                
                if indx.size == 0:
                    raise Exception('')

                ## collect metrics of the target in the last sector
                if a == 0:
                    ### all human-vetted targets
                    for namevarbfstr in listnamevarbfstr:
                        gdat.dictpopl['fstrastr'][namevarbfstr].append(listgdat.dictpopl['fstrastr'][indxtseclast][namevarbfstr][indx[0]])
                if a == 1:
                    ### all human-vetted targets
                    for namevarbfstr in listnamevarbfstr:
                        gdat.dictpopl['fstrhvet'][namevarbfstr].append(listgdat.dictpopl['fstrastr'][indxtseclast][namevarbfstr][indx[0]])
                    ### Ps
                    if dispthis[-1] == 'P':
                        for namevarbfstr in listnamevarbfstr:
                            gdat.dictpopl['fstrpcan'][namevarbfstr].append(listgdat.dictpopl['fstrastr'][indxtseclast][namevarbfstr][indx[0]])
                    ### Fs
                    else:
                        for namevarbfstr in listnamevarbfstr:
                            gdat.dictpopl['fstrfpos'][namevarbfstr].append(listgdat.dictpopl['fstrastr'][indxtseclast][namevarbfstr][indx[0]])
         
            
        for namevarbfstr in listnamevarbfstr:
            gdat.dictpopl['fstrastr'][namevarbfstr] = np.array(gdat.dictpopl['fstrastr'][namevarbfstr])
            if len(gdat.dictpopl['fstrpcan'][namevarbfstr]) == 0:
                raise Exception('')
            gdat.dictpopl['fstrfpos'][namevarbfstr] = np.array(gdat.dictpopl['fstrfpos'][namevarbfstr])
            gdat.dictpopl['fstrpcan'][namevarbfstr] = np.array(gdat.dictpopl['fstrpcan'][namevarbfstr])
            gdat.dictpopl['fstrhvet'][namevarbfstr] = np.array(gdat.dictpopl['fstrhvet'][namevarbfstr])
        
        # crossmatch with TIC to get additional features
        #dictquer = miletos.xmat_tici(ticifstrastr)

        #numbastr = len(ticifstrastr)
        #for namevarbfstr in listnamevarbfstr:
        #    gdat.dictpopl['fstrastr'][namevarbfstr] = np.empty(numbastr)
        #for k, tici in enumerate(ticifstrastr):
        #    indx = np.where(ticifstr == tici)[0]
        #    if indx.size > 0:
        #        for namevarbfstr in listnamevarbfstr:
        #            gdat.dictpopl['fstrastr'][namevarbfstr][k] = gdat.dictpopl['fstrastr'][namevarbfstr][indx[0]]

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

    if gdat.typeanls == 'exopexar':
        gdat.dictpopl['exopexar']['nois'] = ephesus.retr_noistess(gdat.dictpopl['exopexar']['vmagsyst'])
        
    gdat.listnamepopl = np.array(list(gdat.dictpopl.keys()))
    
    if gdat.typeanls == 'exar':
        gdat.lablpoplmast = 'Exoplanets'
    if gdat.typeanls == 'hostexar':
        gdat.lablpoplmast = 'Exoplanet hosts'
    if gdat.typeanls == 'toii':
        gdat.lablpoplmast = 'TOIs'
    if gdat.typeanls == 'hosttoii':
        gdat.lablpoplmast = 'TOI Hosts'
    if gdat.typeanls == 'hosttoii' or gdat.typeanls == 'hostexar':
        namelablsamp = 'namestar'
    if gdat.typeanls == 'toii' or gdat.typeanls == 'exar':
        namelablsamp = 'nameplan'
        
    # fill the list of dictionaries with key carriying the name of the population and elements carrying the label and color
    if gdat.typeanls == 'supntess':
        gdat.listdictlablcolrpopl.append(dict())
        gdat.listdictlablcolrpopl[-1]['cyc1'] = ['Cycle 1', 'gray']
        gdat.listdictlablcolrpopl[-1]['cyc2'] = ['Cycle 2', 'gray']
        gdat.listdictlablcolrpopl[-1]['cyc3'] = ['Cycle 3', 'gray']
    
    if gdat.typeanls.startswith('cosc') or gdat.typeanls.startswith('psys'):
        if gdat.typeanls.startswith('cosc'):
            lablelem = 'BH'
        if gdat.typeanls.startswith('psys'):
            lablelem = 'planet'
        gdat.listdictlablcolrpopl.append(dict())
        gdat.listdictlablcolrpopl[-1][strgpoplstar] = ['All', 'gray']
        gdat.listdictlablcolrpopl[-1][strgpoplstar+'comp'] = ['w/ %s' % lablelem, 'gray']
        gdat.listdictlablcolrpopl[-1][strgpoplstar+'comptran'] = ['w/ transiting %s' % lablelem, 'gray']
        gdat.listdictlablcolrpopl[-1][strgpoplstar+'comptrandete'] = ['w/ detected transiting %s' % lablelem, 'gray']
    
    if gdat.typeanls == 'hosttoii' or gdat.typeanls == 'hostexar' or gdat.typeanls == 'toii' or gdat.typeanls == 'exar':
        gdat.listdictlablcolrpopl.append(dict())
        gdat.listdictlablcolrpopl[-1]['totl'] = ['All', 'gray']
        gdat.listdictlablcolrpopl[-1]['tran'] = ['Transiting', 'blue']
    
        gdat.listdictlablcolrpopl.append(dict())
        gdat.listdictlablcolrpopl[-1]['trantess'] = ['Transiting, TESS discoveries', 'blue']
        gdat.listdictlablcolrpopl[-1]['trantessbrgt'] = ['Bright and transiting TESS discoveries', 'blue']
        gdat.listdictlablcolrpopl[-1]['trantessbrgtsmal'] = ['Bright, small, and transiting TESS discoveries', 'blue']
        gdat.listdictlablcolrpopl[-1]['trantessbrgtsmallcoo'] = ['Bright, small, and transiting TESS discoveries visible from LCO', 'blue']
        gdat.listdictlablcolrpopl[-1]['trantessbrgtsmallcooatmo'] = ['Bright, small, and transiting TESS discoveries visible from LCO, favorable to AC', 'blue']
    
    if gdat.typeanls == 'exar' or gdat.typeanls == 'hostexar':
        gdat.listdictlablcolrpopl.append(dict())
        gdat.listdictlablcolrpopl[-1]['totl'] = ['All', 'gray']
        gdat.listdictlablcolrpopl[-1]['yong'] = ['Young', 'blue']
        
        gdat.listdictlablcolrpopl.append(dict())
        gdat.listdictlablcolrpopl[-1]['totl'] = ['All', 'gray']
        gdat.listdictlablcolrpopl[-1]['gmas'] = ['Precise mass', 'blue']
        
        gdat.listdictlablcolrpopl.append(dict())
        gdat.listdictlablcolrpopl[-1]['totl'] = ['All', 'gray']
        gdat.listdictlablcolrpopl[-1]['gmasinsohabi'] = ['HZ, precise mass', 'blue']
        gdat.listdictlablcolrpopl[-1]['gmasinsomedi'] = ['Medium-irradiated, precise mass', 'blue']
        gdat.listdictlablcolrpopl[-1]['gmasinsohigh'] = ['Highly-irradiated, precise mass', 'blue']
        
    if gdat.typeanls == 'toii' or gdat.typeanls == 'hosttoii':
        gdat.listdictlablcolrpopl.append(dict())
        gdat.listdictlablcolrpopl[-1]['totl'] = ['All', 'gray']
        gdat.listdictlablcolrpopl[-1]['othr'] = ['TOI Process', 'gray']
        gdat.listdictlablcolrpopl[-1]['fstrpmis'] = ['FaintStar PM', 'blue']
        gdat.listdictlablcolrpopl[-1]['fstremis'] = ['FaintStar EM1', 'orange']
    
        gdat.listdictlablcolrpopl.append(dict())
        gdat.listdictlablcolrpopl[-1]['othr'] = ['Other TOIs', 'gray']
        gdat.listdictlablcolrpopl[-1]['pcan'] = ['PC', 'blue']
        gdat.listdictlablcolrpopl[-1]['fpos'] = ['FP', 'green']
        gdat.listdictlablcolrpopl[-1]['knwn'] = ['KP', 'orange']
    
    if gdat.typeanls == 'toii':
        pass

    if gdat.typeanls == 'hosttoii':
        pass

    if gdat.typeanls == 'tceefstr':
        pass
    
    if len(gdat.listdictlablcolrpopl) == 0:
        raise Exception('')
    
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
        listnamefeat[k] = list(gdat.dictpopl[gdat.listnamepopl[k]].keys())
        listlablfeat[k], listscalfeat[k] = tdpy.retr_listlablscalpara(listnamefeat[k])
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
        
        listlablfeattotl[k] = tdpy.retr_labltotl(listlablfeat[k])

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
        numbkeys = len(listnamefeat[k])
        numbtarg = np.empty(numbkeys, dtype=int)
        for n in gdat.indxfeat[k]:
            numbtarg[n] = gdat.dictpopl[namepopl][listnamefeat[k][n]].size
            if numbtarg[n] == 0:
                print('Feature %s of population %s is empty!' % (listnamefeat[k][n], gdat.listnamepopl[k]))
        if np.unique(numbtarg).size > 1:
            print('Number of targets is not the same for every feature.')
            print('k')
            print(k)
            print('namepopl')
            print(namepopl)
            print('listnamefeat[k]')
            print(listnamefeat[k])
            for n in gdat.indxfeat[k]:
                print('gdat.dictpopl[namepopl][listnamefeat[k][n]]')
                summgene(gdat.dictpopl[namepopl][listnamefeat[k][n]])
            raise Exception('')
        
    # filter some features of the population
    listsampfilt = [[] for k in gdat.indxpopl]
    listlablfeatfilt = [[] for k in gdat.indxpopl]
    listlablfeattotlfilt = [[] for k in gdat.indxpopl]
    listnamefeatfilt = [[] for k in gdat.indxpopl]
    listscalfeatfilt = [[] for k in gdat.indxpopl]
    gdat.dictpoplfilt = dict()

    for k in gdat.indxpopl:
        gdat.dictpoplfilt[gdat.listnamepopl[k]] = dict()
        
        for n in gdat.indxfeat[k]:
            #if gdat.typeanls.startswith('psys') and not listnamefeat[k][n] in ['cosi', 'densstar', 'massstar', 'nois', 'peri', 'numbcompstarmean', 'radiplan', \
            #                                                                                                            'radistar', 'rmag', 'rsma', 'smax', 'distsyst']:
            #    continue
            
            boolgood = False
            if (gdat.typeanls == 'hostexar' or gdat.typeanls == 'hosttoii' or gdat.typeanls == 'exar' or gdat.typeanls == 'toii') and listnamefeat[k][n] in [ \
                                                                                 'rascstar', 'declstar', 'loecstar', 'laecstar', \
                                                                                 'radistar', 'massstar', 'metastar', 'loggstar', 'tagetar', 'distsyst', 'numbplantranstar', \
                                                                                 'tmagsyst', 'vmagsyst']:
                boolgood = True
            
            if (gdat.typeanls == 'exar' or gdat.typeanls == 'toii') and listnamefeat[k][n] in [ \
                                                                                 'radiplan', 'massplan', 'densplan', 'tmptplan', 'stdvmassplan', 'vsec', \
                                                                                 'peri', 'duratran', 'dcyc', 'dept', \
                                                                                 'inso', \
                                                                                 'tsmm', 'esmm', 'esmmacwg', 'tsmmacwg']:
                boolgood = True
        
            elif gdat.typeanls == 'toii' and listnamefeat[k][n] in [ \
                                                                      'toii', \
                                                                     ]:
                boolgood = True
        
            if not boolgood:
                continue

            gdat.dictpoplfilt[gdat.listnamepopl[k]][listnamefeat[k][n]] = gdat.dictpopl[gdat.listnamepopl[k]][listnamefeat[k][n]]
            
            listlablfeattotlfilt[k].append(listlablfeattotl[k][n])
            
            # exclude features with string value
            if listnamefeat[k][n] in ['typedisptess', 'strgcomm', 'namestar', 'namesyst', 'nameplan', 'facidisc', 'strgprovmass', 'strgrefrradiplan', 'strgrefrmassplan', 'name']:
                continue
            
            # exclude features with unuseful IDs
            if listnamefeat[k][n] in ['tici']:
                continue
            
            listsampfilt[k].append(np.array(gdat.dictpopl[gdat.listnamepopl[k]][listnamefeat[k][n]]).astype(float))
            listlablfeatfilt[k].append(listlablfeat[k][n])
            listnamefeatfilt[k].append(listnamefeat[k][n])
            listscalfeatfilt[k].append(listscalfeat[k][n])
        
        if len(listsampfilt[k]) == 0:
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
            raise Exception('')

        listsampfilt[k] = np.vstack(listsampfilt[k]).T
    listsamp = listsampfilt
    listlablfeat = listlablfeatfilt
    listlabltotlfeat = listlablfeattotlfilt
    listnamefeat = listnamefeatfilt
    listscalfeat = listscalfeatfilt
    
    for k in gdat.indxpopl:
        numbfeat[k] = len(listnamefeat[k])
        gdat.indxfeat[k] = np.arange(numbfeat[k])
        listnamefeat[k] = np.array(listnamefeat[k])

    # store features on disc
    for k in gdat.indxpopl:
        path = gdat.pathdata + 'dict_%s.csv' % gdat.listnamepopl[k] 
        print('Writing to %s...' % path)
        pd.DataFrame.from_dict(gdat.dictpoplfilt[gdat.listnamepopl[k]]).to_csv(path, header=listlablfeattotlfilt[k], index=False, float_format='%.8g')
    
    typeplottdim = 'best'
    
    # visualize each population separately
    #for k in gdat.indxpopl:
    #    if namelablsamp is None:
    #        listlablsamp = None
    #    else:
    #        listlablsamp = gdat.dictpopl[gdat.listnamepopl[k]][namelablsamp]
    #    path = gdat.pathimag + gdat.listnamepopl[k] + '/'
    #    tdpy.plot_grid(path, gdat.listnamepopl[k], listsamp[k], listlablfeat[k], typeplottdim=typeplottdim, boolplotindi=True, boolplotpair=True, \
    #                                                            # a bool parameter is giving error when this is True due to lower and upper limits being equal
    #                                                            boolplottria=False, \
    #                                                            listscalpara=listscalfeat[k], typefileplot=typefileplot, listnamepara=listnamefeat[k], listlablsamp=listlablsamp)
    
    print('gdat.listnamepopl')
    print(gdat.listnamepopl)
    print('gdat.listdictlablcolrpopl')
    print(gdat.listdictlablcolrpopl)
    
    # plot populations together
    numbplotcomm = len(gdat.listdictlablcolrpopl)
    indxplotcomm = np.arange(numbplotcomm)
    gdat.listnamepoplcomm = [[] for e in indxplotcomm]
    gdat.listcolrpoplcomm = [[] for e in indxplotcomm]
    gdat.listlablpoplcomm = [[] for e in indxplotcomm]
    
    print('Comparing populations...')
    for e in indxplotcomm:
        #gdat.listnamepoplcomm[e] = list(gdat.listdictlablcolrpopl[e].keys())
        #gdat.listnamepoplcomm[e] = np.array(gdat.listnamepoplcomm[e])
        
        gdat.listnamepoplcomm[e] = list(gdat.listdictlablcolrpopl[e].keys())
        indxfrst = np.where(gdat.listnamepoplcomm[e][0] == gdat.listnamepopl)[0][0]
        
        print('gdat.listnamepoplcomm[e]')
        print(gdat.listnamepoplcomm[e])
        for n in gdat.indxfeat[indxfrst]:
            print('Feature: %s' % listnamefeat[indxfrst][n])
            for k in range(len(gdat.listnamepoplcomm[e])):
                print('gdat.listnamepoplcomm[e][k]')
                print(gdat.listnamepoplcomm[e][k])
                print('gdat.listnamepopl')
                print(gdat.listnamepopl)
                indxseco = np.where(gdat.listnamepoplcomm[e][k] == gdat.listnamepopl)[0][0]
                summgene(gdat.dictpopl[gdat.listnamepopl[indxseco]][listnamefeat[indxseco][n]], boolslin=True)
            print('')
        print('')
        print('')
        numbpoplcomm = len(gdat.listnamepoplcomm[e])
        if numbpoplcomm > 0:
            ## find the indices of the populations to be plotted together
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
            gdat.listnamefeatcommtemp = np.unique(gdat.listnamefeatcommtemp)
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
            for uu, u in enumerate(indxpoplcomm):
                gdat.listsampcomm[uu] = listsamp[u][:, indxfeatcomm[uu]]
            
            # number of samples in each population
            numbsamppopl = np.empty(numbpoplcomm, dtype=int)
            for uu in range(numbpoplcomm):
                numbsamppopl[uu] = gdat.listsampcomm[uu].shape[0]
            
            # determine the list of labels for compared populations
            for uu, u in enumerate(indxpoplcomm):
                gdat.listlablpoplcomm[e].append(gdat.listdictlablcolrpopl[e][gdat.listnamepopl[u]][0])
            
            print('gdat.listdictlablcolrpopl')
            print(gdat.listdictlablcolrpopl)
            print('gdat.listlablpoplcomm')
            print(gdat.listlablpoplcomm)

            # determine the list of colors for compared populations
            for uu, u in enumerate(indxpoplcomm):
                gdat.listcolrpoplcomm[e].append(gdat.listdictlablcolrpopl[e][gdat.listnamepopl[u]][1])
            
            print('gdat.listnamepoplcomm[e]')
            print(gdat.listnamepoplcomm[e])
            print('gdat.listlablpoplcomm[e]')
            print(gdat.listlablpoplcomm[e])
            print('gdat.listcolrpoplcomm[e]')
            print(gdat.listcolrpoplcomm[e])
            
            strgplot = ''
            for uu, u in enumerate(indxpoplcomm):
                strgplot += gdat.listnamepoplcomm[e][uu]
            
            tdpy.plot_grid(gdat.pathimag, strgplot, gdat.listsampcomm, gdat.listlablfeatcomm, boolplotindi=True, boolplotpair=True, boolplottria=False, typeplottdim=typeplottdim, \
                           typefileplot=typefileplot, \
                           listnamepara=gdat.listnamefeatcomm, \
                           listscalpara=gdat.listscalfeatcomm, \
                           listlablpopl=gdat.listlablpoplcomm[e], \
                           listcolrpopl=gdat.listcolrpoplcomm[e], \
                           )
            
        
        if gdat.typeanls.startswith('cosc') or gdat.typeanls[12:].startswith('exopmock'):
            # variables
            listnamevarbreca = ['massbhol', 'peri', 'incl']
            listlablvarbreca = [['$M_c$', '$M_\odot$'], ['$P$', 'day'], ['$i$', '$^o$']]
            
            # selection imposed by transit requirement
            listnamevarbprec = None
            listlablvarbprec = None
            boolposirele = np.zeros(gdat.dictnumbsamp['comp'], dtype=bool)
            boolposirele[gdat.dictindxsamp['comp']['comptran']] = True
            boolreleposi = np.ones(gdat.dictnumbsamp['comptran'], dtype=bool)
            listvarbreca = np.vstack([gdat.dictpopl['comp']['massbhol'], gdat.dictpopl['comp']['peri'], gdat.dictpopl['comp']['incl']]).T
            listvarbprec = None
            tdpy.plot_recaprec(gdat.pathimag, 'tran', listvarbreca, listvarbprec, listnamevarbreca, listnamevarbprec, \
                                                         listlablvarbreca, listlablvarbprec, boolposirele, boolreleposi, strgreca='Transit Fraction')
            
            # ROC of detections
            listnamevarbprec = ['sdee']
            listlablvarbprec = [['SDE', '']]
            boolposirele = np.zeros(gdat.dictnumbsamp['comptran'], dtype=bool)
            boolposirele[gdat.dictindxsamp['comptran']['comptrandete']] = True
            boolreleposi = np.ones(gdat.dictindxsamp['comptran']['comptrandete'].size, dtype=bool)
            listvarbreca = np.vstack([gdat.dictpopl['comp']['massbhol'], gdat.dictpopl['comp']['peri'], gdat.dictpopl['comp']['incl']]).T
            listvarbprec = np.vstack([gdat.dictpopl['comptrandete']['sdee']]).T
            tdpy.plot_recaprec(gdat.pathimag, 'comptrandete', listvarbreca, listvarbprec, listnamevarbreca, listnamevarbprec, \
                                                                                            listlablvarbreca, listlablvarbprec, boolposirele, boolreleposi)
            
            # overall selection
            listnamevarbprec = ['sdee']
            listlablvarbprec = [['SDE', '']]
            boolposirele = np.zeros(gdat.dictnumbsamp['comp'], dtype=bool)
            boolposirele[indxtranoccu[gdat.dictindxsamp['comptran']['comptrandete']]] = True
            boolreleposi = np.ones(gdat.dictindxsamp['comptran']['comptrandete'].size, dtype=bool)
            listvarbreca = np.vstack([gdat.dictpopl['comp']['massbhol'], gdat.dictpopl['comp']['peri'], gdat.dictpopl['comp']['incl']]).T
            listvarbprec = np.vstack([gdat.dictpopl['comptrandete']['sdee']]).T
            tdpy.plot_recaprec(gdat.pathimag, 'totl', listvarbreca, listvarbprec, listnamevarbreca, listnamevarbprec, \
                                                                         listlablvarbreca, listlablvarbprec, boolposirele, boolreleposi)
            

            # overall selection
            listnamevarbprec = listnamevarbreca
            listlablvarbprec = listlablvarbreca
            boolposirele = np.zeros(gdat.dictnumbsamp['comp'], dtype=bool)
            boolposirele[indxtranoccu[gdat.dictindxsamp['comptran']['comptrandete']]] = True
            boolreleposi = np.ones(gdat.dictindxsamp['comptran']['comptrandete'].size, dtype=bool)
            listvarbreca = np.vstack([gdat.dictpopl['comp']['massbhol'], gdat.dictpopl['comp']['peri'], gdat.dictpopl['comp']['incl']]).T
            listvarbprec = np.vstack([gdat.dictpopl['comptrandete']['massbhol'], gdat.dictpopl['comptrandete']['peri'], gdat.dictpopl['comptrandete']['incl']]).T
            listvarbdete = []
            tdpy.plot_recaprec(gdat.pathimag, 'comp', listvarbreca, listvarbprec, listnamevarbreca, listnamevarbprec, \
                                                                         listlablvarbreca, listlablvarbprec, boolposirele, boolreleposi, listvarbdete=[])
            

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

