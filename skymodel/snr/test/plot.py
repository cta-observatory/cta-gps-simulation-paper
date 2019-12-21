#import sys, os
from numpy import *
from matplotlib.pyplot import *
from astropy.table import Table
from astropy import units as u
#from astropy.io import fits
from crab import *


files=[]
files=files+ [ [ '../ALL_FILES_0/results_0.txt', 'current test'    ] ]
#files=files+ [ [ '../ALL_FILES_0/results_1.txt', 'current test'    ] ]
#files=files+ [ [ 'results_2.txt', 'current test'    ] ]
#files=files+ [ [ 'results_3.txt', 'current test'    ] ]

#files=files+ [ [ '191211_b/results_0.txt', 'Efficiency= 0.005, alpha=2.25' ] ]
#files=files+ [ [ '191211_c/results_0.txt', 'Efficiency= 0.005, KNEE= 300.'  ] ]
#files=files+ [ [ '191211_d/results_0.txt', 'Efficiency= 0.005, alpha=2.3'  ] ]


for fil in files:

  q = Table.read(fil[0], format='ascii')
  en=q['E[TeV]'][0:40].data               # TeV

  emin=1. ;  emax=200.
  intf=[]

  for i in arange(int(len(q)/40)):

    sed=q['diff_spectrum'][i*40:(i+1)*40]   #  TeV cm-2 s-1  ??
    fl=sed / en**2.                         #  cm-2 s-1 TeV-1 ??
    #plot(en,fl,'y')
    
    intf=intf+[ integ(en,fl, emin=1.,emax=10.)[0] ]   #   cm-2 s-1
    
  intf=array(intf)
  logNlogS(intf,label=fil[1])

  print('N. of objects with flux larger than Crab : ',len(where(intf > 2e-11 )[0]))


### from Pierre (2018)
 
#p=Table.read('old/ctadc_skymodel_gps_sources_snr_3.ecsv',format='ascii.ecsv')
#logNlogS(p['flux_1_10'],label='P3')


### Real SNRs

re=Table.read('real.txt',format='ascii')
rflux=re['Flux_above_1TeV']
errorbar(rflux/1.6,re['n'],yerr=sqrt(re['n']),fmt='-o',label='Real SNR')

#icu = intf / crab( energy=[emin,emax])[0].to_value('cm-2 s-1')
#hist(intf)

#plot(en,crab(energy=en,giveme='sed').to('TeV cm-2 s-1')  )
#plot(en,crab(energy=en,giveme='flux') )

xlabel('Flux above 1 TeV [ph cm-2 s-1]')  ; ylabel('Number of objects')
legend()
loglog()
savefig('test_SNRs_vs_real.pdf')

