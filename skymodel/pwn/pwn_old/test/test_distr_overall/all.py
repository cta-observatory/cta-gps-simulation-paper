#!/usr/bin/env python3
"""
Created on Fri May 31 12:09:55 2020
@author: Andrea  ( andrea.giuliani@inaf.it )
"""

import sys, os
from numpy import *
from matplotlib.pyplot import *
from astropy.table import Table 
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
import aglib as ag

home = os.environ['HOME'] 
gpsmodels = home + '/GitHub/cta-gps-simulation-paper/skymodel/'

pwnefile  = gpsmodels+'pwn/PWNe_final_population.txt'
snrsfile  = gpsmodels+'snr/ALL_FILES_0/results_0.txt'
isnrsfile = gpsmodels+'int-snr/out/d.fits'

synt = Table()
 
# In[PWN]

pwne = Table.read(pwnefile, format='ascii') 

synt['Glon'] = pwne['GLON']
synt['Glat'] = pwne['GLAT']

synt['Distance'] = pwne['distance']
synt['Distance'].unit = 'kpc'

synt['Class'] = 'Pwn'
#synt['F01'] = pwne['ph_flux_above100GeV']
#synt['F01'].unit = 'cm-2 s-1'

ff0=[]
ff1=[]
ff10=[]

for pwn in pwne:
  
  sp=Table.read(gpsmodels+'/pwn/xml/'+pwn['filename'][10:],format='ascii')
  en = sp['col1']   #  MeV
  f  = sp['col2']   #  ph cm-2 s-1 MeV-1
  f0 =  ag.integ(en,f,emin=1e5,emax=1e9)[0]
  f1 =  ag.integ(en,f,emin=1e6,emax=1e9)[0]
  f10 =  ag.integ(en,f,emin=1e7,emax=1e9)[0]
  ff1.append( f1 )  
  ff10.append( f10 )  
  ff0.append( f0 )
  
synt['F01'] = ff0
synt['F01'].unit = 'cm-2 s-1'

synt['F1'] = ff1
synt['F1'].unit = 'cm-2 s-1'

synt['F10'] = ff10
synt['F10'].unit = 'cm-2 s-1'


# In[SNR]  
  
q = Table.read(snrsfile, format='ascii')
en=q['E[TeV]'][0:40].data               # TeV

for i in arange(int(len(q)/40)):

    sed=q['diff_spectrum'][i*40:(i+1)*40]   #  TeV cm-2 s-1  
    fl=sed / en**2.                         #  cm-2 s-1 TeV-1 
     
    #age = q['age[kyr]'][i*40]
    #typ = q['type'][i*40]
    #print(i, 'type:',typ, 'age', age.round(2))
    pos = SkyCoord( q['POS_Y'][i*40], -q['POS_X'][i*40], q['POS_Z'][i*40], unit='kpc',frame='galactocentric')

    synt['Glon'][-1] = pos.galactic.l.value
    synt['Glat'][-1] = pos.galactic.b.value

    synt.add_row()    
    synt['Class'][-1] = 'Snr'
     
    if (sum(fl) != 0):
     
      synt['F01'][-1] = ag.integ(en,fl, emin=0.1 ,emax=1000.)[0]      #   cm-2 s-1
      synt['F1'][-1]  = ag.integ(en,fl, emin=1.0 ,emax=1000.)[0]      #   cm-2 s-1
      synt['F10'][-1] = ag.integ(en,fl, emin=10. ,emax=1000.)[0]     #   cm-2 s-1

# In[iSNR]

isnrs = Table.read(isnrsfile)

for i in arange(len(isnrs)) :

  isnr=isnrs[i]   

  synt.add_row()    
  synt['Class'][-1] = 'iSR'
  synt['Glat'][-1] = isnr['Glat']
  synt['Glon'][-1] = isnr['Glon']
  synt['Distance'][-1] = isnr['Distance']
  
  sp=Table.read(gpsmodels+'/int-snr/out/isnr'+str(i)+'_spec.txt' ,format='ascii')
  en = sp['col1']   #  MeV
  f  = sp['col2']   #  ph cm-2 s-1 MeV-1 

  synt['F01'][-1]  =  ag.integ(en,f,emin=1e5,emax=1e9)[0]
  synt['F1'][-1]   =  ag.integ(en,f,emin=1e6,emax=1e9)[0]
  synt['F10'][-1]  =  ag.integ(en,f,emin=1e7,emax=1e9)[0]
  
  # synt['F01'][-1] = isnr['Flux_int_100GeV']
  # synt['F1'][-1]  = isnr['Flux_int_1TeV']
  # synt['F10'][-1] = isnr['Flux_int_10TeV']
  

# In[out]

synt.write('synt.fits')



