#import sys, os
from numpy import *
from matplotlib.pyplot import *
from astropy.table import Table
from astropy import units as u
#from astropy.io import fits
from crab import *


files=[]
files=files+ [ [ 'results_0.txt', 'current 0'    ] ]
files=files+ [ [ 'results_1.txt', 'current 1'    ] ]
files=files+ [ [ 'results_2.txt', 'current 2'    ] ]

#files=files+ [ [ '../200208/ALL_FILES_0/results_0.txt', 'Kep 1e-4, K 1000'    ] ]
 
for fil in files:

  q = Table.read(fil[0], format='ascii')
  en=q['E[TeV]'][0:40].data               # TeV

  emin=1. ;  emax=200.
  intf=[]
  intf100=[]
  aalp=[]


  for i in arange(int(len(q)/40)):

    sed=q['diff_spectrum'][i*40:(i+1)*40]   #  TeV cm-2 s-1  ??
    fl=sed / en**2.                         #  cm-2 s-1 TeV-1 ??
     
    age = q['age[kyr]'][i*40]
     
    if (sum(fl) != 0) :

      intf100=intf100+[ integ(en,fl, emin=.1,emax=100.,prec=1e4)[0] ]   #   cm-2 s-1
      
  intf100 = array(intf100)
  intf100 = intf100[~isnan(intf100)]
  
  logNlogS(intf100,label=fil[1])

  print('N. of objects brighter than .1 Crab : ', len(where(intf100 > 5.6e-11 )[0]) ) 



### Real SNRs

re=Table.read('real100.txt',format='ascii')
rflux100=re['Flux100GeV']

logNlogS(rflux100,label='Real SNR')
#logNlogS(intf100,label=fil[1])
xlabel('Flux above 100 GeV [ph cm-2 s-1]')  ; ylabel('Number of objects')
legend() ; loglog()
savefig('test_2.pdf')

