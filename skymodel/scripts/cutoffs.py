import matplotlib.pyplot as plt
import numpy as np
from psrqpy import QueryATNF
from astropy.table import Table
from astropy.coordinates import SkyCoord
from scipy.optimize import fsolve
import astropy.units as u

# see document in known-sources/docs/Emax_SNR_PWN.pdf for references to equations

# load SNRCat
snrtable = Table.read('../known-sources/external-input/SNRcat20191008-SNR.csv',delimiter=';',comment='#')
snrobs = Table.read('../known-sources/external-input/SNRcat20191008-OBS.csv',delimiter=';',comment='#')

# Sedov-Taylor time in yr for SN type I and II, Eq. 13-14
t_ST = [235, 84.5]
# Sedov-Taylor radius in pc for SN type I and II, Eq. 8-9
R_ST = [2, 1]
# ejecta density profile power-law index for SN type I and II, Eq. 7
k = [7, 9]
# environmental parameter for SN type I and II, paragraph after Eq. 1
m = [0, 2]
# rest energy of proton in eV
E0 = 1.e9
# rest energy of electron in eV
E0e = 0.511e6
# expansion smoothing parameter
a = -5
# electron charge, esu
e = 4.8e-10
# speed of light, cm/s
c = 3e10
# CR acceleration efficiency
csi_CR = 0.1
# ISM density, g/cm3
rho_ISM = 1.6e-24
# mass loss in massive stars, Msun/yr
Mdot = 1.e-5
# massive star wind velocity, km/s
vW = 10
# solar mass, g
Msun = 1.99e33
# Thomson cross section, cm-2
sigmaT = 0.665e-24
# unit conversion
pc2cm = 3.086e18
yr2sec = 86400 * 365
eV2erg = 1.602e-12


#### methods to estimate cutoff

def lambda_exp(type):
    """
    Calculate remnant expansion indices, eq. 5
    :param type: type of SN, 1 or 2
    :return: lambda_ED, lambda_ST, expansion indices
    """

    # convert type 1 and 2 into list index
    s = type - 1

    lambda_ED = (k[s] - 3) / (k[s] - m[s])
    lambda_ST = 2 / (5 - m[s])

    return lambda_ED, lambda_ST


def radius(t, type):
    """
    Calculate radius of SNR as a function of time, Eq. 4
    :param t: time, yr
    :param type: type of SN, 1 or 2
    :return: R, radius of SNR in pc
    """

    # convert type 1 and 2 into list index
    s = type - 1

    # expansion indices
    lambda_ED, lambda_ST = lambda_exp(type)

    # radius
    R = (np.power(t / t_ST[s], a * lambda_ED) + np.power(t / t_ST[s], a * lambda_ST))
    R = np.power(R, 1. / a)
    R *= R_ST[s]

    return R


def vs(t, type):
    """
    Calculate SNR shock velocity as a function of time, Eq. 5
    :param t: time, yr
    :param type: type of SN, 1 or 2
    :return: v, shock velocity, in cm/s
    """
    # convert type 1 and 2 into list index
    s = type - 1

    # expansion indices
    lambda_ED, lambda_ST = lambda_exp(type)

    # radius
    R = radius(t, type)

    # shock velocity
    v = lambda_ED * np.power(t / t_ST[s], a * lambda_ED - 1) + \
        lambda_ST * np.power(t / t_ST[s], a * lambda_ST - 1)
    v *= np.power(R / R_ST[s], 1 - a) * R_ST[s] / t_ST[s]

    v *= pc2cm
    v /= yr2sec

    return v


def rho(t, type):
    """
    calculate medium density as a function of time
    :param t: time, yr
    :param type: type of SN, 1 or 2
    :return: dens, density in g/cm3
    """
    if type == 1:
        dens = rho_ISM  # ISM density (1 p/cm-3)
    elif type == 2:
        R = radius(t, type)
        # Eq 12
        dens = Mdot / (4 * np.pi * R ** 2 * vW)
        # convert to g/cm2
        # Mdot Msun/yr -> g/s
        dens *= Msun / yr2sec
        # R pc -> cm
        dens /= pc2cm ** 2
        # vW km/s -> cm/s
        dens /= 1.e5
    else:
        print('ERROR: SN of type {} not implemented'.format(type))

    return dens


def Emax_p(t, type, index):
    """
    calculate proton maximum energy
    :param t: time, yr
    :param type: type of SN, 1 or 2
    :param index: spectral index (positive)
    :return: emax, maximum energy in eV
    """
    # convert type 1 and 2 into list index
    s = type - 1

    # Equation 1
    Psi = 2 * e / ((4 - m[s]) * 5 * c * E0 * eV2erg) * csi_CR * np.power(vs(t, type), 2) * \
          np.sqrt(4 * np.pi * rho(t, type) * np.power(pc2cm * radius(t, type), 2))

    # Equation 2
    if index == 2:
        # case with index = 2
        # approximate ln(Emax/E0) = 13.8 for Emax = 10^15 eV
        emax = E0 * Psi / 13.8
    else:
        beta = index - 2
        # approximate 1 - (E0/Emax)^ beta as 1 - 1.e-6^beta
        emax = E0 * np.power(beta / (1 + beta) * (1. - np.power(1.e-6, beta)) * Psi,
                             1. / (1 + beta))

    return emax


def B(t, type):
    """
    Calculate SNR downstream magnetic field, Eq. 16
    :param t: time, yr
    :param type:  type of SN, 1 or 2
    :return: B, magnetic field intensity in G
    """

    b = 1.e-3 * np.power(csi_CR / 0.1, 0.5) * np.power(rho(t, type) / rho_ISM, 0.5) * np.power(
        vs(t, type) / 1.e9, 1.5)

    return b


def Emax_rad(t, type):
    """
    Calculate maximum energy from radiative losses, Eq. 21
    :param t: time, yr
    :param type:  type of SN, 1 or 2
    :return: emax, maximum energy in eV
    """
    emax = 27 * np.pi /16 * e / (sigmaT * B(t, type))
    emax = np.sqrt(emax)
    emax *= E0e * vs(t, type) / c

    return emax

def Emax_e(t,type,index):
    """
    Calculate maximum energy for electrons, Eq. 17
    :param t: time, yr
    :param type:  type of SN, 1 or 2
    :return: emax, maximum energy in eV
    """

    emax = np.minimum(Emax_p(t, type, index), Emax_rad(t, type))

    return emax

def Emax_PSR(Edot):
    """
    Calculate maximum energy of electrons accelerated by pulsar
    based on potential drop
    :param Edot: erg/s
    :return: maximum energy in eV
    """

    emax = 1.7e15 * np.power(Edot/1.e36,0.5)

    return emax

def BTS(t,Edot,tau=1.e3,q=2.5,eta=0.5):
    """
    Calculate PWN magnetic field at termination shock, Eq. 30
    :param t: yr, age of the system
    :param Edot: erg/s
    :param tau: yr, characteristic spin-down time
    :param q: function of braking index
    :param eta: fraction of total energy carried by PSr wind
    :return: magnetic field, G
    """

    Edot0 = Edot * np.power(1 + t/tau, q)
    bts = 3.4e-5 * np.power(eta/0.5,0.5) * np.power(t/1000,-13./10) * np.power(1 + 2 * t / (3 * tau), -7./10)
    return bts

def Emax_TS(t,Edot):
    """
    Radiation limited maximum energy at PWN termination shock, Eq. 32
    :param t: yr, age of the system
    :param Edot: current spindown power, erg/s
    :return: maximum energy, eV
    """

    emax = E0e * np.sqrt(9 * np.pi * e / (sigmaT * BTS(t,Edot)))
    return emax

def Emax_pwn(t,Edot):
    """
    Maximum particle energy in PWN
    Minimum between pulsar potential drop
    and radiation-limited maximum energy at TS
    :param t: yr, age of the system
    :param Edot: current spindown power, erg/s
    :return: maximum energy, eV
    """

    emax = np.minimum(Emax_PSR(Edot),Emax_TS(t,Edot))
    return emax


#### methods to set PWN cutoffs querying ATNF catalog

def get_atnf_version():

    return QueryATNF().get_version

#### cutoff methods

def get_random_cutoff(emin = 10, emax=100):
    # return random cutoff in range between emin and emax
    # flat log distribution
    rnd = np.random.random()
    logemin = np.log10(emin)
    logemax = np.log10(emax)
    ecut = np.power(10,logemin + rnd * (logemax - logemin))
    return ecut

def get_pwn_cutoff(ra,dec,rad_search=2., deathline = 1.e34):
    """
    Set cutoff for PWN based on association with pulsars in ATNF catalog
    Only pulsars with Edot > 1.e34 erg/s are considered, because cutoffs at ~1 TeV
    should have already been detected
    :param ra: R.A. of PWN center, deg
    :param dec: Dec. of PWN center, deg
    :param rad_search: search radius, deg
    :param deathline: minimum Edot, erg/s
    :return:
    """
    # convert coordinates to astropy SkyCoord
    c = SkyCoord(ra, dec, frame='icrs', unit='deg')
    # extract ra, dec in the format required by psrqpy
    ra_hms, dec_dms = c.to_string('hmsdms').split(' ')
    # define circular search region
    search_region = [ra_hms, dec_dms, rad_search]
    # query ATNF catalog
    psrs = QueryATNF(params=['JNAME', 'RAJ', 'DECJ', 'EDOT', 'AGE'], circular_boundary=search_region,
                     condition = 'EDOT > {}'.format(deathline))
    if len(psrs) == 0:
        # no PSR found, setting random cutoff
        ecut = get_random_cutoff()
    else:
        if len(psrs) == 1:
            # 1 pulsar found
            s = 0
        else:
            # multiple pulsars found
            # pulsars position in SkyCoord form
            cpsrs = SkyCoord(ra=psrs['RAJ'], dec=psrs['DECJ'], frame='icrs',
                             unit=(u.hourangle, u.deg))
            # calculate angular separation between pulsars and PWN
            sep = cpsrs.separation(c)
            # select closest pulsar
            s = np.where(sep == np.min(sep))[0][0]
        # assume E_gamma = 0.03 E_e in deep KN regime, convert to TeV
        ecut = 0.01 * 1.e-12 * Emax_pwn(psrs['AGE'][s], psrs['EDOT'][s])

    return ecut

def get_snr_cutoff(ra,dec,name=None, hess=False, index = 2.):

    # identify source in snrcat using coordinates if name not available
    if name == None:
        #print('select SNR by coordinates')
        c = SkyCoord(ra, dec, frame='icrs', unit='deg')
        csnrs = SkyCoord(ra=snrtable['J2000_ra (hh:mm:ss)'], dec = snrtable['J2000_dec (dd:mm:ss)'], frame='icrs',
                         unit=(u.hourangle, u.deg))
        sep = csnrs.separation(c)
        snr = snrtable[sep == np.min(sep)]
        name = snr['G']
    # otherwise, if name is available reformat it to query SNRCat directly
    else:
        #print('select SNR by name')
        # reformat name according to SNRcat conventions
        # drop initial string
        if name[:5] == 'SNR G':
            name = name[5:]
        elif name[0] == 'G':
            name = name[1:]
        # get lat sign, glon and glat
        if '+' in name:
            glatsign = '+'
        else:
            glatsign = '-'
        glon, glat = name.split(glatsign)
        glon = float(glon)
        glat = float(glat)
        name = 'G{:05.1f}{}{:04.1f}'.format(glon, glatsign, glat)
        snr = snrtable[snrtable['G'] == name]
    #print(snr)

    # determine if SN is type I or II
    # HESS objects are all indicated as interacting in the literature, set type 2 and radiation mechanism hadronic
    hadronic = False
    if hess:
        #print("literature says HESS object is interacting, thus type = II, emission hadronic")
        type = 2
        hadronic = True
    # otherwise check SNRcat table to see if interaction with MC or thermal composite emission is reported
    else:
        # check if the SNR type is thermal composite
        if 'thermal' in snr['type']:
            #print("it's a thermal composite SNR, thus type = II, emission hadronic")
            type = 2
            hadronic = True
        else:
            # otherwise check if the observations indicate interactions with MC
            obs = snrobs[snrobs['SNR_id'] == name]
            if 'cloud' in obs['source']:
                #print("interactions with MC, thus type = II, emission hadronic")
                type = 2
                hadronic = True
            else:
                #print("no evidence of ISM interaction, thus set randomly type I or II, emission leptonic")
                dice = np.random.random()
                if dice < 0.2:
                    type = 1
                else:
                    type = 2
                #print("it's type {}".format(type))

    # get age
    # if age not measured derive from size
    if np.isnan(float(snr['age_min (yr)'][0])) or np.isnan(float(snr['age_max (yr)'][0])):
        # get angular size in deg (diameter)
        angsize = snr['size_coarse (arcmin)'][0] / 60.
        # try to get measured distance
        dmin = float(snr['distance_min (kpc)'][0])
        dmax = float(snr['distance_max (kpc)'][0])
        if np.isnan(dmin) and np.isnan(dmax):
            #print('no distance information, set to 1 kpc')
            dist = 1
        elif not np.isnan(dmin) and np.isnan(dmax):
            #print('only dmin known, set to dmin + 1 kpc')
            dist = dmin + 1
        elif np.isnan(dmin) and not np.isnan(dmax):
            #print('only dmax known, set to dmax/2')
            dist = dmax / 2
        else:
            #print('measured distance, take average of dmin and dmax')
            dist = (dmin + dmax) / 2
        #print('angular size {} deg, distance {} kpc'.format(angsize,dist))
        # physical radius in pc
        size = 1.e3 * dist * np.tan(np.radians(angsize)) / 2
        #print('physical size {} pc'.format(size))
        # compare with radius from evolutionary model and derive age
        # define helper function that is 0 for size(age) matching observations
        f = lambda t: radius(t, type) - size
        # find zeros of f, initial guess 1000 years
        age = fsolve(f,1000.)[0]
        #print('inferred age {} yr'.format(age))
    else:
        #otherwise take measured age, select min age for max energy
        age = (float(snr['age_min (yr)'][0]) + float(snr['age_min (yr)'][0])) / 2
        #print('measured age {} yr'.format(age))

    # if estimated age is unrealistically large take a random number between 0.5 and 1 kyr
    if age >= 10000:
        age = 500. + 500. * np.random.random()
    #print('age used for calculation {} yr'.format(age))

    # set cutoff, get maximum particle energy
    if hadronic == True:
        #print('use max p energy')
        #print('particle spectral index {}'.format(index))
        emax = Emax_p(age,type,index)
    else:
        #print('use max e energy')
        # extracting index of electrons from index of gamma requires hypothesis such as Thomson regime etc.
        # use spectral index of 2
        emax = Emax_e(age,type,2.)
    #print('Emax: {} PeV'.format(emax*1.e-15))

    # assume E_gamma = 0.1 E_p (or E_e in KN regime), convert to TeV
    ecut = 0.1 * 1.e-12 * emax
    #print('Ecut: {} TeV'.format(ecut))

    return ecut

def get_cutoff_agn(z):
    """
    Calculate max energy for AGN based on EBL absorption
    :param z: redshift
    :return: cutoff energy, TeV
    """

    # if redshift not measured draw random number following distribution in 3FHL
    if np.isnan(z):
        z = np.random.exponential(scale=0.7)
        # low-z cut to avoid AGN next door
        z = np.maximum(0.1,z)

    # maximum energy allowed by EBL absorption
    # based on Figure 17 of 3FHL paper, analytical approx of max energy
    ecut = 0.04 * np.exp(z**-0.5)

    return ecut


def get_cutoff(ra,dec,classes,name=None, rad_search=4., hess=False, index = 2., z = None):
    if classes == 'PSR':
        ecut = get_pwn_cutoff(ra,dec,rad_search=rad_search)
    elif classes == 'SNR':
        ecut = get_snr_cutoff(ra,dec,name=name, hess = hess, index = index)
        ####### Formulas implemented so far are not valid for the type of old SNRs considered
        ####### set very high cutoff at 10 PeV
        ####### fix when new formulas available
        # ecut = 10000.
    elif classes == 'AGN':
        # based on Fig 17 of 3FHL catalog paper
        ecut = get_cutoff_agn(z=z)
    else:
        ecut = get_random_cutoff(emin=10.,emax=100.)

    return ecut


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Make some diagnostic plots to make sure implementation of formulas is correct
    age = np.logspace(0, 3, 100)  # age in yr

    fig = plt.figure()
    ax = plt.subplot()
    ax.set_xlabel('Age (yr)')
    ax.set_ylabel('Energy (eV)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title('Proton maximum energy')

    for type in [1, 2]:
        for index in [2., 2.3]:
            label = "type {}, spectral index {}".format(type, index)
            ax.plot(age, Emax_p(age, type, index), label=label)

    ax.legend()

    fig = plt.figure()
    ax = plt.subplot()
    ax.set_xlabel('Age (yr)')
    ax.set_ylabel('B (G)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title('Downstream magnetic field')

    for type in [1, 2]:
        label = "type {}".format(type)
        ax.plot(age, B(age, type), label=label)

    ax.legend()

    fig = plt.figure()
    ax = plt.subplot()
    ax.set_xlabel('Age (yr)')
    ax.set_ylabel('Energy (eV)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title('Maximum energy (type I)')

    ax.plot(age, Emax_p(age, 1, 2.), label='age')
    ax.plot(age, Emax_rad(age, 1), label='radiation')

    ax.legend()

    fig = plt.figure()
    ax = plt.subplot()
    ax.set_xlabel('Age (yr)')
    ax.set_ylabel('Energy (eV)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title('Maximum energy (type II)')

    ax.plot(age, Emax_p(age, 2, 2.), label='age')
    ax.plot(age, Emax_rad(age, 2), label='radiation')

    ax.legend()

    fig = plt.figure()
    ax = plt.subplot()
    ax.set_xlabel('Age (yr)')
    ax.set_ylabel('Energy (eV)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title('Electron maximum energy')

    for type in [1, 2]:
        for index in [2., 2.3]:
            label = "type {}, spectral index {}".format(type, index)
            ax.plot(age, Emax_e(age, type, index), label=label)

    ax.legend()

    age = np.logspace(0, 6, 200)  # age in yr
    edot = 1.e38/np.power(1+age/1.e3,2.5)

    fig = plt.figure()

    ax = plt.subplot()
    ax.set_xlabel('Age (yr)')
    ax.set_ylabel('Energy/me')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title('PWN maximum energy')

    ax.plot(age, Emax_PSR(edot)/E0e * np.ones(len(age)), label='PSR potential drop')
    ax.plot(age, Emax_TS(age,1.e38)/E0e, label='radiation')
    ax.plot(age, Emax_pwn(age, edot)/E0e, label='limit')

    ax.legend()

    plt.show()
