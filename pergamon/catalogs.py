"""Dictionary catalog adapters transferred from Nicomedia to Pergamon."""

import os
import fnmatch
import json
from pathlib import Path
import numpy as np
import pandas as pd
import astropy.coordinates
import astropy.time
import astropy.units
import tdpy
from tdpy import summgene
from tdpy.verbosity import print, tqdm
from tdpy.paths import RepositoryPaths
from tdpy.mast import quer_mast, xmat_tici
from .emission_spectroscopy import calc_tsmmesmm

def _repository_data_path(repository, *parts):
    if repository == 'EPHESOS':
        root = Path(os.environ.get('PERGAMON_PATH', Path(__file__).parents[1]))
        path = root / 'data'
        path = path.joinpath(*parts)
        legacy_root = os.environ.get('EPHESOS_PATH')
        if not path.exists() and legacy_root is not None:
            legacy_path = Path(legacy_root) / 'data'
            legacy_path = legacy_path.joinpath(*parts)
            if legacy_path.is_file():
                return str(legacy_path)
        return str(path)
    return os.path.join(str(RepositoryPaths(repository + '_PATH').get_data_path()), *parts)


def retr_toiitici(tici, typeverb=1, dicttoii=None):
    
    if dicttoii is None:
        dicttoii = retr_dicttoii()
    
    toii = None
    indx = np.where(dicttoii['TICID'][0] == tici)[0]
    if indx.size > 0:
        toii = int(str(dicttoii['TOIID'][0][indx[0]]).split('.')[0])
        if typeverb > 0:
            print('Matched the input TIC ID with TOI-%d.' % toii)
    
    return toii


def retr_dictpopltic8( \
                      typepopl, \
                      numbsyst=None, \
                      # type of verbosity
                      ## -1: absolutely no text
                      ##  0: no text output except critical warnings
                      ##  1: minimal description of the execution
                      ##  2: detailed description of the execution
                      typeverb=1, \
                      # Boolean flag to turn on diagnostic mode
                      booldiag=True, \
                     ):
    '''
    Get a dictionary of the sources in the TIC8 with the fields in the TIC8.
    
    Keyword arguments   
        typepopl: type of the population
            'TIC_prmshcon': TIC targets with contamination larger than
            'TIC_m060': TIC targets brighter than TESS magnitude 6.0
            'TIC_m100': TIC targets brighter than TESS magnitude 10.0
            'TIC_m140': TIC targets brighter than TESS magnitude 14.0
            'targtess_prms_ffimm060': TESS targets observed during PM on FFIs brighter than mag 6.0
            'targtess_prms_2min': 2-minute TESS targets obtained by merging the SPOC 2-min bulk downloads

    Returns a dictionary with keys:
        rasc: RA
        decl: declination
        magtsystTESS: TESS magnitude
        radistar: radius of the star
        massstar: mass of the star
    '''
    
    if typeverb > 0:
        print('Retrieving a dictionary of TIC8 for population %s...' % typepopl)
    
    if typepopl.startswith('CTL'):
        strgtypepoplsplt = typepopl.split('_')
        
        if booldiag:
            if len(strgtypepoplsplt) < 3:
                print('')
                print('')
                print('')
                print('typepopl')
                print(typepopl)
                print('strgtypepoplsplt')
                print(strgtypepoplsplt)
                raise Exception('len(strgtypepoplsplt) < 2')
        
        strgtimetess = strgtypepoplsplt[1]
        if strgtimetess == 'yr01':
            listtsec = np.arange(1, 14) # [1-13]
        elif strgtimetess == 'S1':
            listtsec = np.arange(1, 2) # [1]
        elif strgtimetess == 'yr02':
            listtsec = np.arange(13, 27) # [13-26]
        elif strgtimetess == 'yr03':
            listtsec = np.arange(27, 40) # [27-39]
        elif strgtimetess == 'yr04':
            listtsec = np.arange(40, 56) # [40-55]
        elif strgtimetess == 'yr05':
            listtsec = np.arange(56, 70) # [56-69]
        elif strgtimetess == 'sc01':
            listtsec = np.arange(1, 2)
        elif strgtimetess == 'prms':
            listtsec = np.arange(1, 27)
        elif strgtimetess == 'e1ms':
            listtsec = np.arange(27, 56)
        elif strgtimetess == 'e2ms':
            listtsec = np.arange(56, 70)
        else:
            print('typepopl')
            print(typepopl)
            raise Exception('')
        numbtsec = len(listtsec)
        indxtsec = np.arange(numbtsec)

    path = _repository_data_path('EPHESOS', 'listticidata', 'listticidata_%s.csv' % typepopl)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if not os.path.exists(path):
        
        # dictionary of strings that will be keys of the output dictionary
        dictnamefeatmast = dict()
        dictnamefeatmast['ID'] = 'TICID'
        dictnamefeatmast['ra'] = 'rascstar'
        dictnamefeatmast['dec'] = 'declstar'
        dictnamefeatmast['Tmag'] = 'magtsystTESS'
        dictnamefeatmast['rad'] = 'radistar'
        dictnamefeatmast['mass'] = 'massstar'
        dictnamefeatmast['Teff'] = 'tmptstar'
        dictnamefeatmast['logg'] = 'loggstar'
        dictnamefeatmast['MH'] = 'metastar'
        liststrg = list(dictnamefeatmast.keys())
        
        print('typepopl')
        print(typepopl)
        if typepopl.startswith('CTL'):
            
            if strgtypepoplsplt[2] == '20sc':
                strgurll = '_20s_'
                labltemp = '20-second'
            elif strgtypepoplsplt[2] == '2min':
                strgurll = '_'
                labltemp = '2-minute'
            else:
                print('')
                print('')
                print('')
                print('typepopl')
                print(typepopl)
                raise Exception('')

            dictquer = dict()
            listtici = []
            for o in indxtsec:
                if typepopl.endswith('bulk'):
                    pathtess = _repository_data_path('TESS', 'lcur', 'sector-%02d' % listtsec[o])
                    listnamefile = fnmatch.filter(os.listdir(pathtess), '*.fits')
                    listticitsec = []
                    for namefile in listnamefile:
                        listticitsec.append(str(int(namefile.split('-')[2])))
                    listticitsec = np.array(listticitsec)
                else:
                    urlo = 'https://tess.mit.edu/wp-content/uploads/all_targets%sS%03d_v1.csv' % (strgurll, listtsec[o])
                    c = pd.read_csv(urlo, header=5)
                    listticitsec = c['TICID'].values
                    listticitsec = listticitsec.astype(str)
                numbtargtsec = listticitsec.size
                
                if typeverb > 0:
                    print('%d observed %s targets in Sector %d...' % (numbtargtsec, labltemp, listtsec[o]))
                
                if numbtargtsec > 0:
                    dictquertemp = xmat_tici(listticitsec, dictnamefeatmast)
                else:
                    dictquertemp = {name: [np.array([]), ''] for name in dictnamefeatmast.values()}
                
                if o == 0:
                    dictquerinte = dict()
                    for name in dictnamefeatmast.keys():
                        dictquerinte[dictnamefeatmast[name]] = [[] for o in indxtsec]
                
                for name in dictnamefeatmast.keys():
                    print('dictquertemp')
                    print(dictquertemp)
                    dictquerinte[dictnamefeatmast[name]][o].extend(dictquertemp[dictnamefeatmast[name]][0])
                    print('name')
                    print(name)
                    print('dictquerinte[dictnamefeatmast[name]]')
                    summgene(dictquerinte[dictnamefeatmast[name]])
                    print('')

            print('Concatenating arrays from different sectors...')
            for name in dictnamefeatmast.keys():
                dictquer[dictnamefeatmast[name]] = [np.concatenate(dictquerinte[dictnamefeatmast[name]]), '']
            
            u, indxuniq, cnts = np.unique(dictquer['TICID'][0], return_index=True, return_counts=True)
            for name in dictnamefeatmast.keys():
                dictquer[dictnamefeatmast[name]][0] = dictquer[dictnamefeatmast[name]][0][indxuniq]
            dictquer['numbtsec'] = [cnts, '']

        elif typepopl.startswith('TIC'):
            if typepopl.endswith('hcon'):
                request = {'service':'Mast.Catalogs.Filtered.Tic.Rows', 'format':'json', 'params':{ \
                'columns':'ID, ra, dec, Tmag, rad, mass, contratio', \
                                                             'filters':[{'paramName':'contratio', 'values':[{"min":10., "max":1e3}]}]}}
            elif typepopl.endswith('m060'):
                request = {'service':'Mast.Catalogs.Filtered.Tic.Rows', 'format':'json', 'params':{ \
                'columns':'ID, ra, dec, Tmag, rad, mass', \
                                                             'filters':[{'paramName':'Tmag', 'values':[{"min":-100., "max":6.0}]}]}}
            elif typepopl.endswith('m090'):
                request = {'service':'Mast.Catalogs.Filtered.Tic.Rows', 'format':'json', 'params':{ \
                'columns':'ID, ra, dec, Tmag, rad, mass', \
                                                             'filters':[{'paramName':'Tmag', 'values':[{"min":-100., "max":9.0}]}]}}
            elif typepopl.endswith('m100'):
                request = {'service':'Mast.Catalogs.Filtered.Tic.Rows', 'format':'json', 'params':{ \
                'columns':'ID, ra, dec, Tmag, rad, mass', \
                                                             'filters':[{'paramName':'Tmag', 'values':[{"min":-100., "max":10.0}]}]}}
            elif typepopl.endswith('m140'):
                request = {'service':'Mast.Catalogs.Filtered.Tic.Rows', 'format':'json', 'params':{ \
                'columns':'ID, ra, dec, Tmag, rad, mass', \
                                                             'filters':[{'paramName':'Tmag', 'values':[{"min":-100., "max":14.0}]}]}}
            else:
                raise Exception('')
    
            # temp
            ## this can be alternatively done as
            #catalog_data = Catalogs.query_criteria(catalog='Tic', Vmag=[0., 5], objtype='STAR')
            #print(catalog_data.keys())
            #x_value = catalog_data['ID']
            #y_value = catalog_data['HIP']

            headers, outString = quer_mast(request)
            listdictquer = json.loads(outString)['data']
            if typeverb > 0:
                print('%d matches...' % len(listdictquer))
            dictquer = dict()
            features = {'ID': ('TICID', ''), 'Tmag': ('magtsystTESS', 'mag'),
                        'ra': ('rascstar', 'degree'), 'dec': ('declstar', 'degree'),
                        'rad': ('radistar', r'$R_{\odot}$'), 'mass': ('massstar', r'$M_{\odot}$'),
                        'contratio': ('contratio', '')}
            for name, (namefeat, lablunit) in features.items():
                if not listdictquer or name in listdictquer[0]:
                    dictquer[namefeat] = [np.asarray([row[name] for row in listdictquer]), lablunit]
        else:
            raise Exception('Unrecognized population name: %s' % typepopl)
        
        numbtarg = dictquer['radistar'][0].size
            
        if typeverb > 0:
            print('%d targets...' % numbtarg)
            print('Writing to %s...' % path)
        #columns = ['TICID', 'radi', 'mass']
        dictquertemp = dict()
        listnamefeat = list(dictquer.keys())
        listhead = []
        for namefeat in listnamefeat:
            strghead = namefeat
            if dictquer[namefeat][1] != '':
                strghead += ' [%s]' % dictquer[namefeat][1]
            listhead.append(strghead)
            dictquertemp[namefeat] = dictquer[namefeat][0]
        print('dictquertemp')
        print(dictquertemp)
        print('listhead')
        print(listhead)
        pd.DataFrame.from_dict(dictquertemp).to_csv(path, header=listhead, index=False)#, columns=columns)
    else:
        if typeverb > 0:
            print('Reading from %s...' % path)
        dictquertemp = pd.read_csv(path, nrows=numbsyst).to_dict(orient='list')
        dictquer = dict()
        for heading, values in dictquertemp.items():
            name, separator, label = heading.partition(' [')
            dictquer[name] = [np.asarray(values), label[:-1] if separator else '']

    if numbsyst is not None:
        for feature in dictquer.values():
            feature[0] = feature[0][:numbsyst]

    #if gdat.typedata == 'simuinje':
    #    indx = np.where((~np.isfinite(gdat.dictfeat['true']['ssys']['massstar'])) | (~np.isfinite(gdat.dictfeat['true']['ssys']['radistar'])))[0]
    #    gdat.dictfeat['true']['ssys']['radistar'][indx] = 1.
    #    gdat.dictfeat['true']['ssys']['massstar'][indx] = 1.
    #    gdat.dictfeat['true']['totl']['magtsystTESS'] = dicttic8['magtsystTESS']
        
    # check if dictquer is properly defined, whose leaves should be a list of two items (of values and labels, respectively)
    if booldiag:
        for namefeat in dictquer:
            if len(dictquer[namefeat]) != 2 or not isinstance(dictquer[namefeat][1], str):
                print('')
                print('')
                print('')
                print('namefeat')
                print(namefeat)
                
                print('dictquer[namefeat]')
                print(dictquer[namefeat])
                print('typepopl')
                print(typepopl)
                raise Exception('dictquer is not properly defined.')

    return dictquer


def retr_dictpoplrvel(typeverb=1, path=None):
    
    if typeverb > 0:
        print('Reading Sauls Gaia high RV catalog...')
    if path is None:
        path = _repository_data_path('TROIA', 'Gaia_high_RV_errors.txt')
    for line in open(path):
        listnamesaul = line[:-1].split('\t')
        break
    if typeverb > 0:
        print('Reading from %s...' % path)
    data = np.loadtxt(path, skiprows=1)
    dictcatl = dict()
    dictcatl['rascstar'] = data[:, 0]
    dictcatl['declstar'] = data[:, 1]
    dictcatl['stdvrvel'] = data[:, -4]
    
    return dictcatl


def retr_dicthostplan(namepopl, \
                      # type of verbosity
                      ## -1: absolutely no text
                      ##  0: no text output except critical warnings
                      ##  1: minimal description of the execution
                      ##  2: detailed description of the execution
                      typeverb=1, \
                      ):
    
    path = _repository_data_path('EPHESOS', 'dicthost%s.csv' % namepopl)
    if os.path.exists(path):
        if typeverb > 0:
            print('Reading from %s...' % path)
        dicthost = pd.read_csv(path).to_dict(orient='list')
        for name in dicthost.keys():
            if name in ('epocmtracomp', 'pericomp', 'duratrantotl', 'radicomp', 'masscomp'):
                dicthost[name] = [np.asarray(json.loads(value)) for value in dicthost[name]]
            else:
                dicthost[name] = np.array(dicthost[name])
        
    else:
        dicthost = dict()
        if namepopl in ('TOIID', 'toii'):
            dictplan = retr_dicttoii()
        else:
            dictplan = retr_dictexar()
        dictplan = {name: value[0] if isinstance(value, list) and len(value) == 2 else value
                for name, value in dictplan.items()}
        listnamestar = np.unique(dictplan['namestar'])
        dicthost['namestar'] = listnamestar
        numbstar = listnamestar.size
        indxstar = np.arange(numbstar)
        
        listnamefeatstar = ['numbplanstar', 'numbplantranstar', 'radistar', 'massstar']
        listnamefeatcomp = ['epocmtracomp', 'pericomp', 'duratrantotl', 'radicomp', 'masscomp']
        for namefeat in listnamefeatstar:
            dicthost[namefeat] = np.empty(numbstar)
        for namefeat in listnamefeatcomp:
            dicthost[namefeat] = [[] for k in indxstar]
        for k in indxstar:
            indx = np.where(dictplan['namestar'] == listnamestar[k])[0]
            for namefeat in listnamefeatstar:
                if namefeat == 'numbplanstar':
                    dicthost[namefeat][k] = indx.size
                elif namefeat == 'numbplantranstar':
                    dicthost[namefeat][k] = np.sum(dictplan['booltran'][indx])
                else:
                    dicthost[namefeat][k] = dictplan[namefeat][indx[0]]
            for namefeat in listnamefeatcomp:
                dicthost[namefeat][k] = dictplan[namefeat][indx]
                
        print('Writing to %s...' % path)
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        cached = {name: [json.dumps(values.tolist()) for values in data] if name in listnamefeatcomp else data
              for name, data in dicthost.items()}
        pd.DataFrame.from_dict(cached).to_csv(path, index=False)

    return dicthost


def retr_dicttoii(toiitarg=None, boolreplexar=False, \
                  
                  strgelem='comp', \
                  
                  # type of verbosity
                  ## -1: absolutely no text
                  ##  0: no text output except critical warnings
                  ##  1: minimal description of the execution
                  ##  2: detailed description of the execution
                  typeverb=1, \
                  
                  ):
    
    dictfact = tdpy.retr_factconv()
    
    pathephe = str(Path(os.environ.get('PERGAMON_PATH', Path(__file__).parents[1]))) + os.sep
    os.makedirs(pathephe + 'data', exist_ok=True)
    pathexof = _repository_data_path('EPHESOS', 'exofop_toilists.csv')
    if typeverb > 0:
        print('Reading from %s...' % pathexof)
    objtexof = pd.read_csv(pathexof, skiprows=0)
    
    strgradielem = 'radi' + strgelem
    strgstdvradi = 'stdv' + strgradielem
    strgmasselem = 'mass' + strgelem
    strgstdvmass = 'stdv' + strgmasselem
    
    strgperielem = 'peri' + strgelem
    
    dicttoii = {}
    tdpy.setp_dict(dicttoii, 'TOIID', objtexof['TOI'].values)
    numbcomp = dicttoii['TOIID'][0].size
    indxcomp = np.arange(numbcomp)
    toiitargexof = np.empty(numbcomp, dtype=object)
    for k in indxcomp:
        toiitargexof[k] = int(dicttoii['TOIID'][0][k])
        
    if toiitarg is None:
        indxcomp = np.arange(numbcomp)
    else:
        indxcomp = np.where(toiitargexof == toiitarg)[0]
    
    dicttoii['TOIID'][0] = dicttoii['TOIID'][0][indxcomp]
    
    numbcomp = indxcomp.size
    
    if indxcomp.size == 0:
        if typeverb > 0:
            print('The host name, %s, was not found in the ExoFOP TOI Catalog.' % toiitarg)
        return None
    else:
        tdpy.setp_dict(dicttoii, 'namestar', np.empty(numbcomp, dtype=object))
        tdpy.setp_dict(dicttoii, 'nametoii', np.empty(numbcomp, dtype=object))
        for kk, k in enumerate(indxcomp):
            dicttoii['nametoii'][0][kk] = 'TOI-' + str(dicttoii['TOIID'][0][kk])
            dicttoii['namestar'][0][kk] = 'TOI-' + str(dicttoii['TOIID'][0][kk])[:-3]
        
        tdpy.setp_dict(dicttoii, 'depttrancomp', objtexof['Depth (ppm)'].values[indxcomp] * 1e-3, 'ppt')
        tdpy.setp_dict(dicttoii, 'rratcomp', np.sqrt(dicttoii['depttrancomp'][0] * 1e-3))
        tdpy.setp_dict(dicttoii, strgradielem, objtexof['Planet Radius (R_Earth)'][indxcomp].values, r'R$_{\oplus}$')
        tdpy.setp_dict(dicttoii, strgstdvradi, objtexof['Planet Radius (R_Earth) err'][indxcomp].values, r'R$_{\oplus}$')
        
        rascstarstrg = objtexof['RA'][indxcomp].values
        declstarstrg = objtexof['Dec'][indxcomp].values
        tdpy.setp_dict(dicttoii, 'rascstar', np.empty(numbcomp))
        tdpy.setp_dict(dicttoii, 'declstar', np.empty(numbcomp))
        for k in range(dicttoii[strgradielem][0].size):
            objt = astropy.coordinates.SkyCoord('%s %s' % (rascstarstrg[k], declstarstrg[k]), unit=(astropy.units.hourangle, astropy.units.deg))
            dicttoii['rascstar'][0][k] = objt.ra.degree
            dicttoii['declstar'][0][k] = objt.dec.degree

        # a string holding the comments
        tdpy.setp_dict(dicttoii, 'strgcomm', objtexof['Comments'][indxcomp].values)
        
        # transit duration
        tdpy.setp_dict(dicttoii, 'duratrantotl', objtexof['Duration (hours)'].values[indxcomp], 'hours')
        
        # coordinates
        objticrs = astropy.coordinates.SkyCoord(ra=dicttoii['rascstar'][0], dec=dicttoii['declstar'][0], frame='icrs', unit='deg')
        
        ## galactic longitude
        tdpy.setp_dict(dicttoii, 'lgalstar', np.array([objticrs.galactic.l])[0, :], 'degrees')
        
        ## galactic latitude
        tdpy.setp_dict(dicttoii, 'bgalstar', np.array([objticrs.galactic.b])[0, :], 'degrees')
        
        ## ecliptic longitude
        tdpy.setp_dict(dicttoii, 'loecstar', np.array([objticrs.barycentricmeanecliptic.lon.degree])[0, :], 'degrees')
        
        ## ecliptic latitude
        tdpy.setp_dict(dicttoii, 'laecstar', np.array([objticrs.barycentricmeanecliptic.lat.degree])[0, :], 'degrees')

        # SNR
        tdpy.setp_dict(dicttoii, 's2nr', objtexof['Planet SNR'][indxcomp].values)
        
        # amount of follow-up data
        tdpy.setp_dict(dicttoii, 'numbobsvtime', objtexof['Time Series Observations'][indxcomp].values)
        tdpy.setp_dict(dicttoii, 'numbobsvspec', objtexof['Spectroscopy Observations'][indxcomp].values)
        tdpy.setp_dict(dicttoii, 'numbobsvimag', objtexof['Imaging Observations'][indxcomp].values)
        
        # alert year
        tdpy.setp_dict(dicttoii, 'yearaler', objtexof['Date TOI Alerted (UTC)'][indxcomp].values)
        for k in range(len(dicttoii['yearaler'][0])):
            dicttoii['yearaler'][0][k] = astropy.time.Time(dicttoii['yearaler'][0][k] + ' 00:00:00.000').decimalyear
        dicttoii['yearaler'][0] = dicttoii['yearaler'][0].astype(float)

        tdpy.setp_dict(dicttoii, 'tsmmacwg', objtexof['TSM'][indxcomp].values)
        tdpy.setp_dict(dicttoii, 'esmmacwg', objtexof['ESM'][indxcomp].values)
    
        tdpy.setp_dict(dicttoii, 'facidisc', np.empty(numbcomp, dtype=object))
        dicttoii['facidisc'][0][:] = 'Transiting Exoplanet Survey Satellite (TESS)'
        
        tdpy.setp_dict(dicttoii, strgperielem, objtexof['Period (days)'][indxcomp].values)
        dicttoii[strgperielem][0][np.where(dicttoii[strgperielem][0] == 0)] = np.nan

        tdpy.setp_dict(dicttoii, 'epocmtra' + strgelem, objtexof['Epoch (BJD)'][indxcomp].values)

        tdpy.setp_dict(dicttoii, 'magtsystTESS', objtexof['TESS Mag'][indxcomp].values)
        tdpy.setp_dict(dicttoii, 'stdvmagtsystTESS', objtexof['TESS Mag err'][indxcomp].values)

        # transit duty cycle
        tdpy.setp_dict(dicttoii, 'dcyc', dicttoii['duratrantotl'][0] / dicttoii[strgperielem][0] / 24.)
        
        liststrgfeatstartici = ['massstar', 'magtsystVbnd', 'magtsystJbnd', 'magtsystHbnd', 'magtsystKbnd', 'distsyst', 'metastar', 'radistar', 'tmptstar', 'loggstar']
        liststrgfeatstarticiinhe = ['mass', 'Vmag', 'Jmag', 'Hmag', 'Kmag', 'd', 'MH', 'rad', 'Teff', 'logg']
        
        numbstrgfeatstartici = len(liststrgfeatstartici)
        indxstrgfeatstartici = np.arange(numbstrgfeatstartici)

        for strgfeat in liststrgfeatstartici:
            tdpy.setp_dict(dicttoii, strgfeat, np.zeros(numbcomp))
            tdpy.setp_dict(dicttoii, 'stdv' + strgfeat, np.zeros(numbcomp))
        
        ## crossmatch with TIC
        print('Retrieving TIC columns of TOI hosts...')
        tdpy.setp_dict(dicttoii, 'TICID', objtexof['TIC ID'][indxcomp].values)
        listticiuniq = np.unique(dicttoii['TICID'][0].astype(str))
        request = {'service':'Mast.Catalogs.Filtered.Tic', 'format':'json', 'params':{'columns':"*", \
                                                              'filters':[{'paramName':'ID', 'values':list(listticiuniq)}]}}
        headers, outString = quer_mast(request)
        listdictquer = json.loads(outString)['data']
        for k in range(len(listdictquer)):
            indxtemp = np.where(dicttoii['TICID'][0] == listdictquer[k]['ID'])[0]
            if indxtemp.size == 0:
                raise Exception('')
            for n in indxstrgfeatstartici:
                dicttoii[liststrgfeatstartici[n]][0][indxtemp] = listdictquer[k][liststrgfeatstarticiinhe[n]]
                dicttoii['stdv' + liststrgfeatstartici[n]][0][indxtemp] = listdictquer[k]['e_' + liststrgfeatstarticiinhe[n]]
        
        tdpy.setp_dict(dicttoii, 'typedisptess', objtexof['TESS Disposition'][indxcomp].values)
        tdpy.setp_dict(dicttoii, 'boolfpos', objtexof['TFOPWG Disposition'][indxcomp].values == 'FP')
        
        # augment
        tdpy.setp_dict(dicttoii, 'numb%sstar' % strgelem, np.zeros(numbcomp))
        boolfrst = np.zeros(numbcomp, dtype=bool)
        for kk, k in enumerate(indxcomp):
            indxcompthis = np.where(dicttoii['namestar'][0][kk] == dicttoii['namestar'][0])[0]
            if kk == indxcompthis[0]:
                boolfrst[kk] = True
            dicttoii['numb%sstar' % strgelem][0][kk] = indxcompthis.size
        
        tdpy.setp_dict(dicttoii, 'numb%stranstar' % strgelem, dicttoii['numb%sstar' % strgelem][0])
        tdpy.setp_dict(dicttoii, 'lumistar', dicttoii['radistar'][0]**2 * (dicttoii['tmptstar'][0] / 5778.)**4)
        tdpy.setp_dict(dicttoii, 'stdvlumistar', dicttoii['lumistar'][0] * np.sqrt((2 * dicttoii['stdvradistar'][0] / dicttoii['radistar'][0])**2 + \
                                                                                            (4 * dicttoii['stdvtmptstar'][0] / dicttoii['tmptstar'][0])**2))
        
        # predicted mass from radii
        path = pathephe + 'data/exofop_toi_mass_saved.csv'
        if not os.path.exists(path):
            dictelemmass = dict()
            dictelemmass[strgmasselem] = np.ones_like(dicttoii[strgradielem][0]) + np.nan
            dictelemmass[strgstdvmass] = np.ones_like(dicttoii[strgradielem][0]) + np.nan
            
            numbsamppopl = 10
            indx = np.where(np.isfinite(dicttoii[strgradielem][0]))[0]
            for n in tqdm(range(indx.size)):
                k = indx[n]
                meanvarb = dicttoii[strgradielem][0][k]
                stdvvarb = dicttoii['stdvradi' + strgelem][0][k]
                
                # if radius uncertainty is not available, assume that it is small, so the mass uncertainty will be dominated by population uncertainty
                if not np.isfinite(stdvvarb):
                    stdvvarb = 1e-3 * dicttoii[strgradielem][0][k]
                else:
                    stdvvarb = dicttoii['stdvradi' + strgelem][0][k]
                
                # sample from a truncated Gaussian
                listradicomp = tdpy.samp_gaustrun(1000, dicttoii[strgradielem][0][k], stdvvarb, 0., np.inf)
                
                # estimate the mass from samples
                listmassplan = tdpy.retr_massfromradi(listradicomp)
                
                dictelemmass[strgmasselem][k] = np.mean(listmassplan)
                dictelemmass[strgstdvmass][k] = np.std(listmassplan)
                
            if typeverb > 0:
                print('Writing to %s...' % path)
            pd.DataFrame.from_dict(dictelemmass).to_csv(path, index=False)
        else:
            if typeverb > 0:
                print('Reading from %s...' % path)
            dictelemmass = pd.read_csv(path).to_dict(orient='list')
            
            for name in dictelemmass:
                dictelemmass[name] = np.array(dictelemmass[name])
                if toiitarg is not None:
                    dictelemmass[name] = dictelemmass[name][indxcomp]
            
            if dicttoii['radistar'][0].size != dictelemmass[strgmasselem].size:
                raise Exception('dictelemmass size is anomalous.')

        tdpy.setp_dict(dicttoii, strgmasselem, dictelemmass[strgmasselem])
        tdpy.setp_dict(dicttoii, strgstdvmass, dictelemmass[strgstdvmass])
        
        perielem = dicttoii[strgperielem][0]
        masselem = dicttoii[strgmasselem][0]

        rvelsemapred = tdpy.retr_rvelsema(perielem, dicttoii['massstar'][0], masselem, 90., 0.)
        tdpy.setp_dict(dicttoii, 'rvelsemapred', rvelsemapred, 'degrees')
        
        masstotl = dicttoii['massstar'][0] + dicttoii[strgmasselem][0] / dictfact['msme']
        tdpy.setp_dict(dicttoii, 'masstotl', masstotl)
        
        tdpy.setp_dict(dicttoii, 'smax' + strgelem, tdpy.retr_smaxkepl(dicttoii[strgperielem][0], dicttoii['masstotl'][0]))
        
        rsmaelem = (dicttoii[strgradielem][0] / dictfact['rsre'] + dicttoii['radistar'][0]) / (dictfact['aurs'] * dicttoii['smax'+strgelem][0])
        tdpy.setp_dict(dicttoii, 'rsma' + strgelem, rsmaelem)
        
        tdpy.setp_dict(dicttoii, 'irra', dicttoii['lumistar'][0] / dicttoii['smax'+strgelem][0]**2)
        
        
        # temp check if factor of 2 is right
        tmptcomp = dicttoii['tmptstar'][0] * np.sqrt(dicttoii['radistar'][0] / dicttoii['smax'+strgelem][0] / 2. / dictfact['aurs'])
        tdpy.setp_dict(dicttoii, 'tmpt' + strgelem, tmptcomp)
        
        # temp check if factor of 2 is right
        stdvtmptcomp = np.sqrt((dicttoii['stdvtmptstar'][0] / dicttoii['tmptstar'][0])**2 + \
                                                        0.5 * (dicttoii['stdvradistar'][0] / dicttoii['radistar'][0])**2) / np.sqrt(2.)
        tdpy.setp_dict(dicttoii, 'stdvtmpt' + strgelem, stdvtmptcomp)
        
        tdpy.setp_dict(dicttoii, 'dens' + strgelem, 5.51 * dicttoii[strgmasselem][0] / dicttoii[strgradielem][0]**3, 'g/cm^3')
        
        tdpy.setp_dict(dicttoii, 'booltran', np.ones(numbcomp, dtype=bool))
        
        tdpy.setp_dict(dicttoii, 'vesc', tdpy.retr_vesc(dicttoii[strgmasselem][0], dicttoii[strgradielem][0]))
        
        print('temp: vsiistar and projoblq are NaNs')
        tdpy.setp_dict(dicttoii, 'vsiistar', np.full(numbcomp, np.nan))
        tdpy.setp_dict(dicttoii, 'projoblq', np.full(numbcomp, np.nan))
        
        # replace confirmed planet features
        if boolreplexar:
            dictexar = retr_dictexar()
            listdisptess = objtexof['TESS Disposition'][indxcomp].values.astype(str)
            listdisptfop = objtexof['TFOPWG Disposition'][indxcomp].values.astype(str)
            indxexofcpla = np.where((listdisptfop == 'CP') & (listdisptess == 'PC'))[0]
            listticicpla = dicttoii['TICID'][0][indxexofcpla]
            numbticicpla = len(listticicpla)
            indxticicpla = np.arange(numbticicpla)
            for k in indxticicpla:
                indxexartici = np.where((dictexar['TICID'][0] == int(listticicpla[k])) & \
                                                    (dictexar['facidisc'] == 'Transiting Exoplanet Survey Satellite (TESS)'))[0]
                indxexoftici = np.where(dicttoii['TICID'][0] == int(listticicpla[k]))[0]
                for strg in dictexar.keys():
                    if indxexartici.size > 0:
                        dicttoii[strg][0] = np.delete(dicttoii[strg][0], indxexoftici)
                    dicttoii[strg][0] = np.concatenate((dicttoii[strg][0], dictexar[strg][indxexartici]))

        # derived quantities
        ## photometric noise in the TESS passband over one hour
        tdpy.setp_dict(dicttoii, 'noistess', tdpy.retr_noisphot(dicttoii['magtsystTESS'][0], 'TESS'))
        
        ## atmospheric characterization
        # calculate TSM and ESM
        calc_tsmmesmm(dicttoii, strgelem=strgelem)
    
        # turn zero TSM ACWG or ESM ACWG into NaN
        indx = np.where(dicttoii['tsmmacwg'][0] == 0)[0]
        dicttoii['tsmmacwg'][0][indx] = np.nan
        
        indx = np.where(dicttoii['esmmacwg'][0] == 0)[0]
        dicttoii['esmmacwg'][0][indx] = np.nan
        
        # surface gravity of the companion
        tdpy.setp_dict(dicttoii, 'logg' + strgelem, dicttoii[strgmasselem][0] / dicttoii[strgradielem][0]**2)

    return dicttoii


def retr_toiifstr():
    '''
    Return the TOI IDs that have been alerted by the FaintStar project
    '''
    
    dicttoii = retr_dicttoii(toiitarg=None, boolreplexar=False, \
                             # type of verbosity
                             ## -1: absolutely no text
                             ##  0: no text output except critical warnings
                             ##  1: minimal description of the execution
                             ##  2: detailed description of the execution
                             typeverb=1, \
                             )
    listtoiifstr = []
    for k in range(len(dicttoii['strgcomm'][0])):
        if isinstance(dicttoii['strgcomm'][0][k], str) and 'found in faint-star QLP search' in dicttoii['strgcomm'][0][k]:
            listtoiifstr.append(dicttoii['nametoii'][0][k][4:])
    listtoiifstr = np.array(listtoiifstr)

    return listtoiifstr


def retr_dictexar( \
                  strgexar=None, \
                  
                  # type of verbosity
                  ## -1: absolutely no text
                  ##  0: no text output except critical warnings
                  ##  1: minimal description of the execution
                  ##  2: detailed description of the execution
                  typeverb=1, \
                  
                  strgelem='comp', \
                 ):
    
    strgradielem = 'radi' + strgelem
    strgstdvradi = 'stdv' + strgradielem
    strgmasselem = 'mass' + strgelem
    strgstdvmass = 'stdv' + strgmasselem
    
    strgstrgrefrradielem = 'strgrefrradi' + strgelem
    strgstrgrefrmasselem = 'strgrefrmass' + strgelem
    
    # variables for which both the quantity and its uncertainty will be retrieved
    liststrgstdv = ['radistar', 'massstar', 'tmptstar', 'loggstar', strgradielem, strgmasselem, 'tmpt'+strgelem, 'tagestar', \
                    'magtsystVbnd', 'magtsystJbnd', 'magtsystHbnd', 'magtsystKbnd', 'magtsystTESS', 'metastar', 'distsyst', 'lumistar', \
                   ]
        
    strgnumbelemstar = 'numb%sstar' % strgelem
    strgnumbelemtranstar = 'numb%stranstar' % strgelem
    # get NASA Exoplanet Archive data
    path = _repository_data_path('EPHESOS', 'PSCompPars_2023.10.12_16.58.05.csv')
    if typeverb > 0:
        print('Reading from %s...' % path)
    objtexar = pd.read_csv(path, skiprows=318)
    if strgexar is None:
        indx = np.arange(objtexar['hostname'].size)
        #indx = np.where(objtexar['default_flag'].values == 1)[0]
    else:
        indx = np.where(objtexar['hostname'] == strgexar)[0]
        #indx = np.where((objtexar['hostname'] == strgexar) & (objtexar['default_flag'].values == 1))[0]
    
    dictfact = tdpy.retr_factconv()
    
    if indx.size == 0:
        
        print('The target name, %s, was **not** found in the NASA Exoplanet Archive composite table.' % strgexar)
        
        return None
    else:
        print('The target name, %s, has been found in the NASA Exoplanet Archive composite table.' % strgexar)
        dictexar = {}

        listnamefeat = ['namestar', 'nameplan', 'TICID', 'rascstar', 'declstar', 'TOIID', 'methdisc', 'eccecomp', 'facidisc', 'yeardisc', \
                        'irra', 'pericomp', 'smaxcomp', 'epocmtracomp', 'cosicomp', 'duratrantotl', 'booltran', 'strgprovmass', strgstrgrefrradielem, strgstrgrefrmasselem, \
                        'vesc', 'masstotl', 'densplan', 'vsiistar', 'projoblq', 'boolcibp', strgnumbelemstar, strgnumbelemtranstar, 'dcyc', 'lgalstar', 'bgalstar', \
                        'loecstar', 'laecstar', 'rratcomp', 'rsmacomp', 'esmm', 'tsmm', 'loggcomp', 'depttrancomp', \
                        
                        # derived
                        'lumibbodbolohost', \
                        'fracxrayboloaddihost', \
                        'fluxbbod0224host', \
                        'lumibbod0224host', \
                        'lumipred0224host', \
                        'fluxpred0224host', \

                        'esmm', 'tsmm', 'stdvesmm', 'stdvtsmm', \
                        ]
        for name in liststrgstdv:
            listnamefeat += [name]
            listnamefeat += ['stdv' + name]
        
        for name in listnamefeat:
            dictexar[name] = [[], '']
        
        dictexar['namestar'][0] = objtexar['hostname'][indx].values
        dictexar['nameplan'][0] = objtexar['pl_name'][indx].values
        
        numbplanexar = len(dictexar['nameplan'][0])
        indxplanexar = np.arange(numbplanexar)

        listticitemp = objtexar['tic_id'][indx].values
        dictexar['TICID'][0] = np.empty(numbplanexar, dtype=int)
        for k in indxplanexar:
            if isinstance(listticitemp[k], str):
                dictexar['TICID'][0][k] = listticitemp[k][4:]
            else:
                dictexar['TICID'][0][k] = 0
        
        dictexar['rascstar'][0] = objtexar['ra'][indx].values
        dictexar['declstar'][0] = objtexar['dec'][indx].values
        
        # err1 have positive values or zero
        # err2 have negative values or zero
        
        dictexar['TOIID'][0] = np.empty(numbplanexar, dtype=object)
        
        # discovery method
        dictexar['methdisc'][0] = objtexar['discoverymethod'][indx].values
        
        # discovery facility
        dictexar['facidisc'][0] = objtexar['disc_facility'][indx].values
        
        # discovery year
        dictexar['yeardisc'][0] = objtexar['disc_year'][indx].values
        
        dictexar['irra'][0] = objtexar['pl_insol'][indx].values
        dictexar['irra'][0][np.where(dictexar['irra'][0] <= 0.)] = np.nan
        
        dictexar['pericomp'][0] = objtexar['pl_orbper'][indx].values # [days]
        
        dictexar['smaxcomp'][0] = objtexar['pl_orbsmax'][indx].values # [AU]
        
        # eccentricity
        dictexar['eccecomp'][0] = objtexar['pl_orbeccen'][indx].values
        
        dictexar['epocmtracomp'][0] = objtexar['pl_tranmid'][indx].values # [BJD]
        dictexar['cosicomp'][0] = np.cos(objtexar['pl_orbincl'][indx].values / 180. * np.pi)
        dictexar['duratrantotl'][0] = objtexar['pl_trandur'][indx].values # [hour]
        dictexar['depttrancomp'][0] = 10. * objtexar['pl_trandep'][indx].values # ppt
        
        dictexar['booltran'][0] = objtexar['tran_flag'][indx].values.astype(bool)
        
        # mass provenance
        dictexar['strgprovmass'][0] = objtexar['pl_bmassprov'][indx].values

        # radius reference
        dictexar[strgstrgrefrradielem][0] = objtexar['pl_rade_reflink'][indx].values
        for a in range(dictexar[strgstrgrefrradielem][0].size):
            if isinstance(dictexar[strgstrgrefrradielem][0][a], float) and not np.isfinite(dictexar[strgstrgrefrradielem][0][a]):
                dictexar[strgstrgrefrradielem][0][a] = ''
        
        # mass reference
        dictexar[strgstrgrefrmasselem][0] = objtexar['pl_bmasse_reflink'][indx].values
        for a in range(dictexar[strgstrgrefrmasselem][0].size):
            if isinstance(dictexar[strgstrgrefrmasselem][0][a], float) and not np.isfinite(dictexar[strgstrgrefrmasselem][0][a]):
                dictexar[strgstrgrefrmasselem][0][a] = ''
        
        for strg in liststrgstdv:
            strgvarbexar = None
            if strg.startswith('magtsyst'):
                strgvarbexar = 'sy_'
                if strg.endswith('TESS'):
                    strgextn = 't'
                elif strg.endswith('bnd'):
                    strgextn = strg[-4].lower()
                strgvarbexar += '%smag' % strgextn
            if strg.endswith('syst'):
                strgvarbexar = 'sy_'
                if strg[:-4] == 'dist':
                    strgvarbexar += 'dist'
            if strg.endswith('star'):
                strgvarbexar = 'st_'
                if strg[:-4] == 'logg':
                    strgvarbexar += 'logg'
                if strg[:-4] == 'tage':
                    strgvarbexar += 'age'
                if strg[:-4] == 'meta':
                    strgvarbexar += 'met'
                if strg[:-4] == 'radi':
                    strgvarbexar += 'rad'
                if strg[:-4] == 'mass':
                    strgvarbexar += 'mass'
                if strg[:-4] == 'tmpt':
                    strgvarbexar += 'teff'
                if strg[:-4] == 'lumi':
                    strgvarbexar += 'lum'
            if strg.endswith('plan') or strg.endswith(strgelem):
                strgvarbexar = 'pl_'
                if strg[:-4].endswith('mag'):
                    strgvarbexar += '%smag' % strg[0]
                if strg[:-4] == 'tmpt':
                    strgvarbexar += 'eqt'
                if strg[:-4] == 'radi':
                    strgvarbexar += 'rade'
                if strg[:-4] == 'mass':
                    strgvarbexar += 'bmasse'
            if strgvarbexar is None:
                print('strg')
                print(strg)
                raise Exception('')
            
            stdv = (objtexar['%serr1' % strgvarbexar][indx].values - objtexar['%serr2' % strgvarbexar][indx].values) / 2.
            
            dictexar[strg][0] = objtexar[strgvarbexar][indx].values
            dictexar['stdv%s' % strg][0] = stdv
            if strg == strgradielem:
                dictexar[strg][1] = r'$R_\oplus$'
                dictexar['stdv%s' % strg][1] = r'$R_\oplus$'
            elif strg == strgmasselem:
                dictexar[strg][1] = r'$M_\oplus$'
                dictexar['stdv%s' % strg][1] = r'$M_\oplus$'
            
            if strg == 'tmptstar':
                dictexar[strg][1] = 'K'
            if strg == 'radistar':
                dictexar[strg][1] = r'R$_{\odot}$'
        
        
        dictexar['vesc'][0] = tdpy.retr_vesc(dictexar[strgmasselem][0], dictexar[strgradielem][0])
        dictexar['masstotl'][0] = dictexar['massstar'][0] + dictexar[strgmasselem][0] / dictfact['msme']
        
        dictexar['densplan'][0] = objtexar['pl_dens'][indx].values # [g/cm3]
        dictexar['vsiistar'][0] = objtexar['st_vsin'][indx].values # [km/s]
        dictexar['projoblq'][0] = objtexar['pl_projobliq'][indx].values # [deg]
        
        # Boolean flag indicating if the planet is part of a circumbinary planetary system
        dictexar['boolcibp'][0] = objtexar['cb_flag'][indx].values == 1
        
        dictexar[strgnumbelemstar][0] = np.empty(numbplanexar)
        dictexar[strgnumbelemtranstar][0] = np.empty(numbplanexar, dtype=int)
        boolfrst = np.zeros(numbplanexar, dtype=bool)
        #dictexar['booltrantotl'][0] = np.empty(numbplanexar, dtype=bool)
        for k, namestar in enumerate(dictexar['namestar'][0]):
            indxexarstar = np.where(namestar == dictexar['namestar'][0])[0]
            if k == indxexarstar[0]:
                boolfrst[k] = True
            dictexar['numb%sstar' % strgelem][0][k] = indxexarstar.size
            indxexarstartran = np.where((namestar == dictexar['namestar'][0]) & dictexar['booltran'][0])[0]
            dictexar['numb%stranstar' % strgelem][0][k] = indxexarstartran.size
            #dictexar['booltrantotl'][0][k] = dictexar['booltran'][0][indxexarstar].all()
        
        # transit duty cycle
        dictexar['dcyc'][0] = dictexar['duratrantotl'][0] / dictexar['pericomp'][0] / 24.
        
        # coordinates
        objticrs = astropy.coordinates.SkyCoord(ra=dictexar['rascstar'][0], \
                                               dec=dictexar['declstar'][0], frame='icrs', unit='deg')
        
        ## galactic longitude
        dictexar['lgalstar'][0] = np.array([objticrs.galactic.l])[0, :]
        
        ## galactic latitude
        dictexar['bgalstar'][0] = np.array([objticrs.galactic.b])[0, :]
        
        ## ecliptic longitude
        dictexar['loecstar'][0] = np.array([objticrs.barycentricmeanecliptic.lon.degree])[0, :]
        
        ## ecliptic latitude
        dictexar['laecstar'][0] = np.array([objticrs.barycentricmeanecliptic.lat.degree])[0, :]

        # radius ratio
        dictexar['rratcomp'][0] = dictexar[strgradielem][0] / dictexar['radistar'][0] / dictfact['rsre']
        
        # sum of the companion and stellar radii divided by the semi-major axis
        dictexar['rsmacomp'][0] = (dictexar[strgradielem][0] / dictfact['rsre'] + dictexar['radistar'][0]) / (dictexar['smaxcomp'][0] * dictfact['aurs'])
        
        # derived properties
        # bolometric flux of a star based on radius and temperature
        dictexar['lumibbodbolohost'] = tdpy.retr_fluxbolobbod(dictexar['tmptstar'], dictexar['radistar'])
        
        # X-ray flux of the host star
        dictxraystar = tdpy.retr_dictxraystar(dictexar['tmptstar'], dictexar['radistar'], dictexar['distsyst'])
        dictexar['fracxrayboloaddihost'] = dictxraystar['fracxrayboloaddi']
        dictexar['fluxbbod0224host'] = dictxraystar['fluxbbod0224']
        dictexar['lumibbod0224host'] = dictxraystar['lumibbod0224']
        dictexar['lumipred0224host'] = dictxraystar['lumipred0224']
        dictexar['fluxpred0224host'] = dictxraystar['fluxpred0224']
        
        # surface gravity of the companion
        dictexar['loggcomp'][0] = np.log10(9.8 * dictexar[strgmasselem][0] / dictexar[strgradielem][0]**2)

        # calculate TSM and ESM
        calc_tsmmesmm(dictexar, strgelem=strgelem)

        indxnonntran = np.where(~dictexar['booltran'][0])[0]
        dictexar['esmm'][0][indxnonntran] = np.nan
        dictexar['tsmm'][0][indxnonntran] = np.nan
        
    return dictexar
