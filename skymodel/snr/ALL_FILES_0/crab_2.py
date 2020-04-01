#import sys, os
from numpy import *
from matplotlib.pyplot import *
#from astropy.table import Table
from astropy import units as u
#from astropy.io import fits





def logNlogS(fluxes,point='--.',label=''):

  f2=sort(fluxes)
  nn= arange(len(f2),0,-1)
  #errorbar(f2,nn,yerr=sqrt(nn),fmt=point,label=label)
  fill_between(f2,nn-sqrt(nn),nn+sqrt(nn),alpha=.3)
  plot(f2,nn,point,label=label)
  return f2,nn






def integ(enref,flref,emin=1.,emax=10., prec=1e4):

  #print('Integ. range : ',emin,emax)

  dd=log10(emax/emin)
  energy_b = 10**( arange(dd*prec+1)/prec )* emin
  de = energy_b[1:]-energy_b[:-1]
  en = sqrt(energy_b[1:]*energy_b[:-1])

  flux=exp(interp(log(en),log(enref),log(flref) ))
 
  #plot(en,flux,'.')
  #plot(enref,flref,'*')
  #loglog()
  #show()
 
  intFlux  = sum(flux*de )
  intFluxE = sum(flux*de *en)

  return [intFlux,intFluxE]



def pl(en,k=1.,index=2,e_piv=1.):
  
  flux = k * (en / e_piv)**(index)
  return flux



def crab(energy=[.1,10],model='hegra', giveme='intFlux'):

  try :
    print('Input energy units :',energy.unit)
  except:
    energy=energy*u.TeV
    print('Energy units set to Tev')


  def hegra(en):   # HEGRA

    k = 2.83e-11  *(u.cm**2 *u.s *u.TeV)**-1  # ph / cm2 s TeV
    index = -2.62
    e_piv = 1. * u.TeV
    flux = pl(en,k=k, index=index, e_piv=e_piv)
  
    return flux


  def amenomori(en):      # Amenomori et al. (PL)

    k = 1.49e-15  *(u.cm**2 *u.s *u.TeV)**-1  # ph / cm2 s TeV
    index = -2.91
    e_piv = 40. * u.TeV
    flux = pl(en,k=k, index=index, e_piv=e_piv)

    return flux


  def hawc(en):          # Hawc model (Log.Par)

    k = 2.35e-13  *(u.cm**2 *u.s*u.TeV)**-1  # typo nel paper ! -->  ph / TeV cm2 s
    alpha = 2.79
    beta = 0.06
    e_piv = 7. *u.TeV
    flux = k * (en / e_piv)**(-alpha-beta*log(en/e_piv))

    return flux
  

  def meyer(en):

    p = array([-10.2708, -0.53616, -0.179475, 0.0473174, 0, -0.00449161])

    ind=0.

    for i in arange(6):

      ind = ind+p[i]*log10(en.to_value('TeV'))**i
      sed = u.Quantity(10.**ind,'erg cm-2 s-1')
      flux = sed.to('TeV cm-2 s-1') / en.to('TeV')**2.

    return flux
    


  crabFlux={'hegra':hegra, 'amenomori':amenomori, 'hawc':hawc, 'meyer':meyer}


  try:
    flux=crabFlux[model](energy)
  except:
    print('Available models : ',crabFlux.keys())
    flux=0.
  
  sed=flux*energy*energy.to('erg')

  results={'flux':flux,'sed':sed}

  if giveme == 'intFlux' :

    emin=energy[0].to('TeV')
    emax=energy[1].to('TeV')

    dd=log10(emax/emin)

    energy_b = 10**( arange(dd*1000.+1.)/1000. )* emin      # TeV
    de = energy_b[1:]-energy_b[:-1]

    en = sqrt(energy_b[1:]*energy_b[:-1])

    flux=crabFlux[model](en)

    intFlux  = sum(flux*de )
    intFluxE = sum(flux*de *en.to('erg'))

    results['intFlux']= [intFlux,intFluxE]


  return results[giveme]

