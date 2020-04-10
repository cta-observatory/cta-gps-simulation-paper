#import sys, os
from numpy import *
from matplotlib.pyplot import *
from astropy.table import Table
#from astropy import units as u
#from astropy.io import fits
#from crab import *
import aglib as ag

outfile= '../ALL_FILES_0/results_0.txt'   ;  label= 'Synt. population'
#outfile = sys.argv[1]  ;  label = sys.argv[2]

q = Table.read(outfile, format='ascii')
en=q['E[TeV]'][0:40].data               # TeV

emin=1. ;  emax=200.
intf=[]
intf100=[]
intf100tev=[]
aalp=[]

for i in arange(int(len(q)/40)):

    sed=q['diff_spectrum'][i*40:(i+1)*40]   #  TeV cm-2 s-1  ?
    fl=sed / en**2.                         #  cm-2 s-1 TeV-1 ?
     
    age = q['age[kyr]'][i*40]
    typ = q['type'][i*40]
    print(i, 'type:',typ, 'age', age.round(2))
     
    if (sum(fl) != 0)*(age < 100.):

      alp = log(fl[8]/fl[0])/log(en[8]/en[0])
      aalp=aalp+[alp]
  
      colo = ['pp','r','g','b','y']
  
      plot(en,sed,colo[typ])
    
      intf       = intf       +[ ag.integ(en,fl, emin=1.  ,emax=100.)[0] ]                  #   cm-2 s-1
      intf100    = intf100    +[ ag.integ(en,fl, emin=.1  ,emax=100.)[0] ]                  #   cm-2 s-1
      intf100tev = intf100tev +[ ag.integ(en,fl, emin=100.,emax=10000.)[0] ]                  #   cm-2 s-1
    
intf=array(intf)
intf=intf[~isnan(intf)]
  
intf100 = array(intf100)
intf100 = intf100[~isnan(intf100)]

intf100tev = array(intf100tev)
intf100tev = intf100tev[~isnan(intf100tev)]
  
#logNlogS(intf100,label=fil[1])

print('N. of objects brighter than 1 Crab (E >.1 and 1 TeV) : ', len(where(intf100 > 5.6e-10 )[0]), len(where(intf > 2e-11 )[0]) )
print('N. of objects brighter than 10^-16 (E>100 TeV) : ', len(where(intf100tev > 1e-16 )[0]) )

# Plot spectra

ylabel('Flux [TeV cm-2 s-1]')  ;  xlabel('Energy [TeV]')
xlim([0.1,2e2]) ; ylim([1e-24,1e-8])
loglog()  
savefig('Spectra.png')
show()

### Real SNRs

re=Table.read('real.txt',format='ascii.fixed_width')
rflux1=re['Flux1TeV']
rflux100=re['Flux100GeV']

ag.logNlogS(rflux1,label='Real SNR')
ag.logNlogS(intf,label=label)
xlabel('Flux above 1 TeV [ph cm-2 s-1]')  ; ylabel('Number of objects')
legend() ; loglog()
savefig('logNlogS_1TeV.png')
show()

ag.logNlogS(rflux100,label='Real SNR')
ag.logNlogS(intf100,label=label)
xlabel('Flux above 100 GeV [ph cm-2 s-1]')  ; ylabel('Number of objects')
legend() ; loglog()
savefig('logNlogS_100GeV.png')
show()

ag.logNlogS(intf100tev,label=label)
xlabel('Flux above 100 TeV [ph cm-2 s-1]')  ; ylabel('Number of objects')
legend() ; loglog()
savefig('logNlogS_100TeV.png')
show()

#sp Indices
 
snr = Table.read('snr_gc.ecsv')
hist(aalp,bins=arange(-4,-1,.1))
hist(-snr['spec_pl_index'][~isnan(snr['spec_pl_index'])],bins=arange(-4,-1,.2) )
xlabel('Index 100 GeV - 1 TeV')
savefig('slopes.png')
show()


