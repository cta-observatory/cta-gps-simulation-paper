import matplotlib.pyplot as plt
import numpy as np

# see document in known-sources/docs/FakePeVatron.pdf for references to equations

# Sedov-Taylor time in yr for SN type I and II, Eq. 12-14
t_ST = [235, 84.5]
# Sedov-Taylor radius in pc for SN type I and II, Eq. 7-8
R_ST = [2, 1]
# ejecta density profile power-law index for SN type I and II, Eq. 6
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
    Calculate radius of SNR as a function of time, Eq. 3
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
    Calculate SNR shock velocity as a function of time, Eq. 4
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
        # Eq 11
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
    Calculate downstream magnetic field, Eq. 15
    :param t: time, yr
    :param type:  type of SN, 1 or 2
    :return: B, magnetic field intensity in G
    """

    b = 1.e-3 * np.power(csi_CR / 0.1, 0.5) * np.power(rho(t, type) / rho_ISM, 0.5) * np.power(
        vs(t, type) / 1.e9, 1.5)

    return b


def Emax_rad(t, type):
    """
    Calculate maximum energy from radiative losses, Eq. 18
    :param t: time, yr
    :param type:  type of SN, 1 or 2
    :return: emax, maximum energy in eV
    """
    emax = 18 * np.pi * e / (sigmaT * B(t, type))
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
    :return: emax, maximum energy in eV
    """

    emax = 1.7e15 * np.power(Edot/1.e36,0.5)

    return emax


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Make some diagnostic plots to make sure implementation is correct
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

    plt.show()
