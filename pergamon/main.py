import os
import sys

import numpy as np
import matplotlib.pyplot as plt

import ephesus
import tdpy
from tdpy import summgene 
import miletos

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


def main():
    
    pass

def cnfg_mock():

    pathimag = os.environ['MILETUS_DATA_PATH'] + '/'
    
    # conversion factors
    factrsrj, factmsmj, factrjre, factmjme, factaurj = ephesus.retr_factconv()

    numbsamp = 1
    indxsamp = np.arange(numbsamp)
    
    cade = 2. / 60. / 24. # [days]
    minmtime = 0.
    maxmtime = 30.
    
    numbtime = int((maxmtime - minmtime) / cade)
    time = np.linspace(minmtime, maxmtime, numbtime)
    
    listradistar = tdpy.icdf_powr(np.random.rand(numbsamp), 0.1, 10., -2.)
    listmassstar = tdpy.icdf_powr(np.random.rand(numbsamp), 0.1, 10., -2.)
    
    listlcur = np.empty((numbsamp, numbtime))
    listnumbplan = np.random.random_integers(1, 10, size=numbsamp)
    for n in indxsamp:
        print('listnumbplan')
        print(listnumbplan)
        listperi = tdpy.icdf_powr(np.random.rand(listnumbplan[n]), 1., 20., -2.)
        listepoc = tdpy.icdf_self(np.random.rand(listnumbplan[n]), minmtime, maxmtime)
        listincl = tdpy.icdf_self(np.random.rand(listnumbplan[n]), 89., 90.)
        listradiplan = tdpy.icdf_powr(np.random.rand(listnumbplan[n]), 0.5, 23., -2.) # R_E
        listmassplan = ephesus.retr_massfromradi(listradiplan / factrjre, boolinptsamp=False)
        listmasstotl = listmassplan + listmassstar[n]
        listsmax = ephesus.retr_smaxkepl(listperi, listmasstotl) * factaurj
        listrsma = (listradiplan + listradistar[n]) / listsmax
        listcosi = np.cos(np.pi * listincl / 180.)
    
        listlcur[n, :] = ephesus.retr_rflxtranmodl(time, listperi, listepoc, listradiplan, listradistar[n], listrsma, listcosi) - 1.
        strgextn = '%04d' % n
        ephesus.plot_lcur(pathimag, strgextn, timedata=time, lcurdata=listlcur[n, :])
        
        listlcur[n, :] += 1e-6 * np.random.randn(numbtime)
        
        # call miletos to analyze data
        listarrytser = dict()
        arry = np.empty((numbtime, 3))
        arry[:, 0] = time
        arry[:, 1] = listlcur[n, :]
        arry[:, 2] = 1e-6
        listarrytser['raww'] = [[[arry]], []]
        labltarg = 'Mock %d' % n
        pathtarg = pathimag + 'Mock%04d/' % n
        dictmile = miletos.main.init( \
                                  listarrytser=listarrytser, \
                                  labltarg=labltarg, \
                                  pathtarg=pathtarg, \
                                 )
        
        print('')

globals().get(sys.argv[1])()



