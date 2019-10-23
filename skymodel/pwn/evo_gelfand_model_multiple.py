import astropy.units as u
from astropy.constants import m_e, c, e, m_p
import matplotlib
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/fiori/GAMERA/lib') #To be modify with the location of GAMERA lib
import gappa as gp
import numpy as np
from scipy.integrate import cumtrapz
from tqdm import tqdm as tqdm


def launch_evo_m(pwn_n, age, l0, e0, mej2, nn, eta, t0, eps, dens, Tfir, Ufir, Tnir, Unir, binss, tmin, ebreak, alpha1, alpha2, dist): #function for launch all the computation (more info in every step of the computation)
    t = np.logspace(np.log10(1), np.log10(age+1.), 100)
    R_snr, R_rs, v_snr, v_rs = Evo_R_V_SNR(e0, mej2, dens, t)
    if age > 2000 and t0 < 1000 and l0 < 1e40: #This is needed since for some set of parameters the bohm diffusion breaks the GAMERA computation (need to be checked!). In such cases we neglect it for the moment. 
        no_escape = True
    else:
        no_escape = False
    R, V, V_ej, LT, B, M, EMAX, SED, P, F, E, DENS_EJ = Evo_pwn(t, l0, e0, mej2, nn, eta, t0, dens, Tfir, Ufir, Tnir, Unir, R_rs, v_rs, R_snr, v_snr, eps, ebreak, alpha1, alpha2, no_escape, pwn_n)
    if R_rs[len(R)-1] <= R[-1]:
        sed = 0
        tot = 0
    else:
        sed, tot = final_spectrum(t, age, LT, B, EMAX, R, V, dens, dist, Tfir, Ufir, Tnir, Unir, binss, tmin, ebreak, alpha1, alpha2, no_escape)
        #save all the results
        save_ph_sed(pwn_n, tot)
    return R_snr, R_rs, R

def save_ph_sed(pwn_n, spectra): #simple function for saving the results
    np.savetxt('ph_sed_pwn_{0}.txt'.format(pwn_n), np.transpose([spectra[:,0], spectra[:,1]]), header='en(TeV)\t\tN2(dN/dE)(erg)')

def Evo_R_V_SNR(e0, mej2, dens, t):
    """
    Self-similar solution for radius and expansion velocity of the forward
    and reverse shock of the SNR.
    Valid until the SNR goes radiative. (Appendix A - Gelfand et al. 2009)

    Parameters
    ----------
    e0 : float
        Energy of the SN explosion (tipically 1e51 erg)
    mej2 : float
        Ejected mass of the SN explosion (gr)
    dens : float
        Hidrogen particle density of the ISM (particle/cm**3)
    t : array-like
        Time steps (in year)
    Returns
    -------
    R_snr : array-like
        Array with the forward shock radius at each time step (Pc)
    R_rs : array-like
        Array with the reverse shock radius at each time step (Pc)
    v_snr : array-like
        Array with the forward shock expansion velocity at each time step (cm/s)
    v_rs : array-like
        Array with the reverse shock expansion velocity at each time step (cm/s)
    """
    rho = dens * (m_p.to('g').value) #density of ISM in gr/cm**3
    #Characteristic lengths
    R_ch = mej2**(1/3) * rho**(-1/3)
    t_ch = e0**(-1/2) * mej2**(5/6) * rho**(-1/3) / gp.yr_to_sec

    #Start calculation
    R_snr = []
    R_rs = []
    v_snr = []
    v_rs = []
    for i in range(len(t)):
        if t[i] <= (0.25*t_ch):
            R = 1.12 * R_ch * (t[i]/t_ch)**(2/3)
            rrs = (1/1.19) * R
            v = 0.75 * (R_ch/t_ch) * (t[i]/t_ch)**(-1/3) / gp.yr_to_sec
            vrs = (1/2.38) * v
            R_snr.append(np.copy(R))
            v_snr.append(np.copy(v))
            R_rs.append(np.copy(rrs))
            v_rs.append(np.copy(vrs))

        if t[i] > (0.25*t_ch) and t[i] <= (0.52*t_ch):
            R = 1.12 * R_ch * (t[i]/t_ch)**(2/3)
            rrs = (1.49-0.16*((t[i]-(0.25*t_ch))/t_ch)-0.46*np.log(t[i]/(0.25*t_ch))) * (R_ch/t_ch)*t[i]
            v = 0.75 * (R_ch/t_ch) * (t[i]/t_ch)**(-1/3) / gp.yr_to_sec
            vrs = (0.5+0.16*((t[i]-(0.25*t_ch))/t_ch))*((R_ch/t_ch)) / gp.yr_to_sec
            R_snr.append(np.copy(R))
            v_snr.append(np.copy(v))
            R_rs.append(np.copy(rrs))
            v_rs.append(np.copy(vrs))

        if t[i] > (0.52*t_ch):
            R = ((1.12*R_ch*(0.52**(2/3)))**(5/2) + (2.026*e0/rho)**(1/2) * (t[i]-(0.52*t_ch)) * gp.yr_to_sec)**(2/5)
            rrs = (1.49-0.16*((t[i]-(0.25*t_ch))/t_ch)-0.46*np.log(t[i]/(0.25*t_ch))) * (R_ch/t_ch)*t[i]
            v = (2/5) * (2.026*e0/rho)**(1/2) * R**(-3/2)
            vrs = (0.5+0.16*((t[i]-(0.25*t_ch))/t_ch))*((R_ch/t_ch)) / gp.yr_to_sec
            R_snr.append(np.copy(R))
            v_snr.append(np.copy(v))
            R_rs.append(np.copy(rrs))
            v_rs.append(np.copy(vrs))

    R_snr = np.array(R_snr) / gp.pc_to_cm
    R_rs = np.array(R_rs) / gp.pc_to_cm
    R_rs[R_rs < 0.0] = 0.0
    v_snr = np.array(v_snr)
    v_rs = np.array((v_rs))
    v_rs[R_rs <= 0.0] = 0.0
    return R_snr, R_rs, v_snr, v_rs

def vel_ej(r, t): #Balistic velocity of the material in the ejecta inside the SNR
    return r/(t*gp.yr_to_sec)

def dens_ej(r, r_snr, r_rs, t, e0, mej2, dens):
    """
    Evolution of the density inside the SNR (as in Blondin et al. 2001).
    Valid only up to the reverse shock of the SNR

    Parameters
    ----------
    r : array-like
        Array of radii where calculate the density (Pc)
    r_snr : array-like
        Array with the forward shock radius at each time step (Pc)
    r_rs : array-like
        Array with the reverse shock radius at each time step (Pc)
    t : array-like
        Time steps (in year)
    e0 : float
        Energy of the SN explosion (tipically 1e51 erg)
    mej2 : float
        Ejected mass of the SN explosion (gr)
    dens : float
        Hidrogen particle density of the ISM (particle/cm**3)
    Returns
    -------
    rho_ej : array-like
        Array of arrays with the density profile at each time step (gr/cm**3)
    """
    r *= gp.pc_to_cm
    r_rs *= gp.pc_to_cm
    r_snr *= gp.pc_to_cm
    v_t = np.sqrt((40./18.)*(e0/mej2)) #Transition velocity
    if (v_t*t*gp.yr_to_sec) <= r_snr:
        r_t = (v_t*t*gp.yr_to_sec)
    else:
        r_t = r_snr
    rho_ej = []
    if r_rs >= r_t:
        for j in range(len(r)):
            if r[j] < r_t:
                rho = (10./(9.*np.pi))*e0*(v_t**(-5))*((t*gp.yr_to_sec)**(-3))
            if r[j] >= r_t and r[j]< r_rs:
                rho = (10./(9.*np.pi))*e0*(v_t**(-5))*(r[j]/(v_t*(t*gp.yr_to_sec)))**(-9)*(((t*gp.yr_to_sec)**(-3)))
            if r[j] >= r_rs and r[j]< r_snr:
                rho = np.inf #just to break the computation if the radius of the PWN became bigger the the one of the RS
            rho_ej.append(rho)
    else:
        for j in range(len(r)):
            if r[j] < r_t:
                rho = (10/(9*np.pi))*e0*(v_t**(-5))*((t*gp.yr_to_sec)**(-3))
            if r[j] >= r_t and r[j]< r_rs:
                rho = (10/(9*np.pi))*e0*(v_t**(-5))*(r[j]/(v_t*(t*gp.yr_to_sec)))**(-9)*(((t*gp.yr_to_sec)**(-3)))
            if r[j] >= r_rs:
                rho = np.inf #just to break the computation if the radius of the PWN became bigger the the one of the RS
            rho_ej.append(rho)
    return np.array(rho_ej)

def rr(r): #simple function to calculate a grid of radii
    return np.logspace(np.log10(1e-5), np.log10(np.max(r)), 1000)

def broken_powerlaw(ebreak,index_low,index_high, emaxt, bins):
    #Particle spectrum injected from the Pulsar
    #(without normalization - GAMERA compute it from the spin down luminosity)
    e = np.logspace(np.log10(gp.m_e), np.log10(3*emaxt[-1]),bins)
    n = np.zeros(len(e))
    e_low = [e<ebreak]
    e_high = [e>=ebreak]
    n[tuple(e_low)] += (e[tuple(e_low)]/ebreak)**-index_low
    n[tuple(e_high)] += (e[tuple(e_high)]/ebreak)**-index_high
    return np.array(list(zip(e,n)))

def particle_spectrum(tmin, tmax, tt, p_spectrum, lt, b, emax, r, v, dens, Tfir, Ufir, Tnir, Unir, no_escape):
    """
    GAMERA computation of the particle spectrum
    http://libgamera.github.io/GAMERA/docs/time_dependent_modeling.html

    Procedure to do at each time step in the calculation of the PWN Radius

    Returns
    -------
    sed : array-like
        Array with the evolved particle spectrum (erg/cm**2/s vs TeV) at the
        last step
    energy : float
        Total particle energy content (in erg)
    """
    fp = gp.Particles()
    t = tt[tt<=tmax]
    fp.SetCustomInjectionSpectrum(p_spectrum)
    if no_escape == False:
        e = np.logspace(np.log10(gp.m_e),np.log10(3*np.max(emax)),100) #particle escape
        t_m, e_m = np.meshgrid(t, e) #particle escape
        fp.SetTimeAndEnergyDependentEscapeTime(t, e, t_esc(e_m, t_m, b, r)) #particle escape
    fp.SetLuminosity(list(zip(t,lt)))
    fp.SetBField(list(zip(t,b)))
    fp.SetEmax(list(zip(t,emax)))
    fp.SetRadius(list(zip(t,r)))
    fp.SetExpansionVelocity(list(zip(t,v)))
    fp.SetAmbientDensity(dens)
    fp.AddThermalTargetPhotons(2.7,0.25*gp.eV_to_erg) #CMB
    fp.AddThermalTargetPhotons(Tfir, Ufir) #FIR photon field
    fp.AddThermalTargetPhotons(Tnir, Unir) #NIR photon field
    fp.SetTmin(tmin)
    erad = np.logspace(-21,4.,250) * gp.TeV_to_erg # energies(in ergs) where radiation will be calculated
    fr = gp.Radiation()
    fr.AddThermalTargetPhotons(2.7,0.25*gp.eV_to_erg) #CMB
    fr.AddThermalTargetPhotons(Tfir, Ufir) #FIR photon field
    fr.AddThermalTargetPhotons(Tnir, Unir) #NIR photon field
    fr.SetAmbientDensity(dens)
    fp.SetAge(tmax)
    fp.ToggleQuietMode()
    fp.CalculateElectronSpectrum()
    sed = np.array(fp.GetParticleSED())
    energy = fp.GetParticleEnergyContent() * gp.TeV_to_erg
    return sed, energy

def particle_spectrum_start(tmin, tmax, p_spectrum, lt0, lt1, b0, b1, emax0, emax1, r0, r1, v0, v1, dens, Tfir, Ufir, Tnir, Unir):
    """
    GAMERA computation of the particle spectrum (for the first two time bins)
    http://libgamera.github.io/GAMERA/docs/time_dependent_modeling.html

    Procedure to do at the first two time steps in the calculation
    of the PWN Radius (seems that GAMERA needs at least three time steps to solve the advective equation)

    Returns
    -------
    sed : array-like
        Array with the evolved particle spectrum (erg/cm**2/s vs TeV) at the
        last step
    energy : float
        Total particle energy content (in erg)
    """
    t = np.linspace(tmin, tmax, 3)
    b = np.linspace(b0, b1, 3)
    lt = np.linspace(lt0, lt1, 3)
    emax = np.linspace(emax0, emax1, 3)
    r = np.linspace(r0, r1, 3)
    v = np.linspace(v0, v1, 3)
    fp = gp.Particles()
    fp.SetCustomInjectionSpectrum(p_spectrum)
    fp.SetLuminosity(list(zip(t,lt)))
    fp.SetBField(list(zip(t,b)))
    fp.SetEmax(list(zip(t,emax)))
    fp.SetRadius(list(zip(t,r)))
    fp.SetExpansionVelocity(list(zip(t,v)))
    fp.SetAmbientDensity(dens)
    fp.AddThermalTargetPhotons(2.7,0.25*gp.eV_to_erg) #CMB
    fp.AddThermalTargetPhotons(Tfir, Ufir) #FIR photon field
    fp.AddThermalTargetPhotons(Tnir, Unir) #NIR photon field
    fp.SetTmin(tmin)
    erad = np.logspace(-21,4.,250) * gp.TeV_to_erg # energies(in ergs) where radiation will be calculated
    fr = gp.Radiation()
    fr.AddThermalTargetPhotons(2.7,0.25*gp.eV_to_erg) #CMB
    fr.AddThermalTargetPhotons(Tfir, Ufir) #FIR photon field
    fr.AddThermalTargetPhotons(Tnir, Unir) #NIR photon field
    fr.SetAmbientDensity(dens)
    fp.SetAge(tmax)
    fp.ToggleQuietMode()
    fp.CalculateElectronSpectrum()
    sed = np.array(fp.GetParticleSED())
    energy = fp.GetParticleEnergyContent() * gp.TeV_to_erg
    return sed, energy

def final_spectrum(t, age, LT, B, EMAX, R, V, dens, dist, Tfir, Ufir, Tnir, Unir, binss, tmin, ebreak, alpha1, alpha2, no_escape):
    """
        GAMERA computation of the particle spectrum (for the extraction of the
        photon sed at the end of the evolution of the PWN)
        http://libgamera.github.io/GAMERA/docs/time_dependent_modeling.html

        Returns
        -------
        sed : array-like
            Array with the evolved particle spectrum (erg/cm**2/s vs TeV) at the
            last step
        tot : array-like
            Array with the total photon spectrum (erg/cm**2/s vs TeV)
        ic :  array-like
            Array with the inverse compton photon spectrum (erg/cm**2/s vs TeV)
        ic :  array-like
            Array with the inverse compton contribution to the total
            photon spectrum (erg/cm**2/s vs TeV)
        ic_cmb :  array-like
            Array with the cmb inverse compton contribution to the total
            photon spectrum (erg/cm**2/s vs TeV)
        ic_fir :  array-like
            Array with the fir inverse compton contribution to the total
            photon spectrum (erg/cm**2/s vs TeV)
        ic_nir :  array-like
            Array with the nir inverse compton contribution to the total
            photon spectrum (erg/cm**2/s vs TeV)
        ic_ssc :  array-like
            Array with the self-synchrotron compton contribution to the total
            photon spectrum (erg/cm**2/s vs TeV)
        ic_synch :  array-like
            Array with the synchrotron contribution to the total
            photon spectrum (erg/cm**2/s vs TeV)
        """
    fp = gp.Particles()
    p_spectrum = broken_powerlaw(ebreak,alpha1,alpha2,EMAX, 500)
    if no_escape == False:
        e = np.logspace(np.log10(gp.m_e),np.log10(3*np.max(EMAX)),100) #particle escape
        t_m, e_m = np.meshgrid(t, e) #particle escape
        fp.SetTimeAndEnergyDependentEscapeTime(t, e, t_esc(e_m, t_m, B, R)) #particle escape
    fp.SetCustomInjectionSpectrum(p_spectrum)
    fp.SetLuminosity(list(zip(t,LT)))
    fp.SetBField(list(zip(t,B)))
    fp.SetEmax(list(zip(t,EMAX)))
    fp.SetRadius(list(zip(t,R)))
    fp.SetExpansionVelocity(list(zip(t,V)))
    fp.SetAmbientDensity(dens)
    fp.AddThermalTargetPhotons(2.7,0.25*gp.eV_to_erg)
    fp.AddThermalTargetPhotons(Tfir, Ufir)
    fp.AddThermalTargetPhotons(Tnir, Unir)
    fp.SetTmin(tmin)
    erad = np.logspace(-2,4.,binss) * gp.TeV_to_erg # energies(in ergs) where radiation will be calculated
    fr = gp.Radiation()
    fr.SetDistance(dist)
    fr.AddThermalTargetPhotons(2.7,0.25*gp.eV_to_erg)
    fr.AddThermalTargetPhotons(Tfir, Ufir)
    fr.AddThermalTargetPhotons(Tnir, Unir)
    fr.SetAmbientDensity(dens)
    fp.SetAge(age)
    fp.ToggleQuietMode()
    fp.CalculateElectronSpectrum(binss)
    sed = np.array(fp.GetParticleSED())
    sp = np.array(fp.GetParticleSpectrum())
    fr.SetElectrons(sp[:])
    fr.SetBField(fp.GetBField())
    fr.AddSSCTargetPhotons(fp.GetRadius())
    fr.ToggleQuietMode()
    fr.CalculateDifferentialPhotonSpectrum(erad)
    tot = np.array(fr.GetTotalSED())
    return sed, tot

def spin_down_lum(l0, t0, t, nn): #simple function for the calculation of the spin down luminosity at each time step
    return  l0*((1+t/t0)**(-1.*(nn + 1.)/(nn - 1.)))

def p_emax(eps, eta, l0, t0, t, nn, b): #calculation of the maximum energy of the electrons at each time step.
    #compute it in two ways and check the smaller value.
    #1- from the condition that the Larmor radius of the electrons inside the PWN
    #   is smaller than the termination shock radius of the PWN
    #2- balancing synchrotron losses and acceleration.
    #   From de Jager et al. 1996 (alpha=1, <sin^2(theta)>=2/3)
    e1 = 3 * eps *  gp.el_charge * np.sqrt((eta*(spin_down_lum(l0, t0, t, nn))/((1-eta) * gp.c_speed)))
    e2 = 6.1e14 * np.sqrt((3*1e-3)/(2*b)) * gp.eV_to_erg
    if e1 < e2:
        return e1
    else:
        #print('Synchrotron losses very high!')
        return e2

def t_esc(e, t, b, r):#Particle escape timescale as in Bohm diffusion (Zhang et al. 2008)
    t_esc = 3.4e4*(b/1e-5)*((e/(10*gp.TeV_to_erg))**-1)*(r**2) * gp.yr_to_sec
    return t_esc

def Evo_pwn(t, l0, e0, mej2, nn, eta, t0, dens, Tfir, Ufir, Tnir, Unir, R_rs, v_rs, R_snr, v_snr, eps, ebreak, alpha1, alpha2, no_escape, pwn_n):
    """
    One-zone time-dependent leptonic model for the evolution of young PWNe.
    The code follow the work of Gelfand et al. 2008 (Sec. 2.2 - step 2 and 3 made by GAMERA)

    Parameters
    ----------
    t : array-like
        Time steps (in year)
    l0 : float
        Initial pulsar spin down luminosity (erg/s)
    e0 : float
        Energy of the SN explosion (tipically 1e51 erg)
    mej2 : float
        Ejected mass of the SN explosion (gr)
    nn : float
        Pulsar braking index
    eta : float
        Part of the spin down power that goes in magnetic field
    t0 :  float
        Initial spin-down timescale of the pulsar
    eps : float
        Containment factor (<1) for the calculation of the maximum electron energies
    ebreak : float
        Break energy of the broken power law of the injection spectrum (erg)
    alpha1 : float
        spectral index at E<E_break of the broken power law of the injection spectrum
    alpha2 : float
        spectral index at E>E_break of the broken power law of the injection spectrum
    dens : float
        Hidrogen particle density of the ISM (particle/cm**3)
    Tfir: float
        Far-infrared radiation temperature (for IC calculation - K)
    Ufir: float
        Far-infrared radiation energy density (for IC calculation - erg/cm**3)
    TNir: float
        Near-infrared radiation temperature (for IC calculation - k)
    Unir: float
        Near-infrared radiation energy density (for IC calculation - erg/cm**3)
    R_snr : array-like
        Array with the forward shock radius at each time step (Pc)
    R_rs : array-like
        Array with the reverse shock radius at each time step (Pc)
    v_snr : array-like
        Array with the forward shock expansion velocity at each time step (cm/s)
    v_rs : array-like
        Array with the reverse shock expansion velocity at each time step (cm/s)
    Returns
        -------
    R : array-like
        Array with the radius of the PWN at each time step (Pc)
    V : array-like
        Array with the expansion velocity of the PWN at each time step (cm/s)
    V_ej : array-like
        Array with the SNR ejecta velocity just outside the PWN at each time step (cm/s)
    LT : array-like
        Array with the spin down luminosity of the PSR at each time step (erg/s)
    B : array-like
        Array with the magnetic field strength inside the PWN at each time step (Gauss)
    M : array-like
        Array with the swept-up mass by the expanding PWN at each time step (gr)
    EMAX : array-like
        Array with the maximum electron energy at each time step (erg)
    SED : array-like
        Array of array of the electron spectrum at each time step (TeV vs erg/cm**2/s)
    F : array-like
        Array of the force applied on the mass shell outiside the PWN at each time step
    E : array-like
        Array of the total particle energy content in the PWN at each time step (erg)
    DENS_EJ : array-like
        Array of the SNR ejecta density just outside the PWN at each time step (gr/cm**3)
    """
    #Pre-Initial condition (needed for calculation of the particle spectrum with GAMERA)
    t00 = 1e-5 #year
    r00 = 1.44 * ((((l0**2)*(e0**3))/(mej2**5))**0.1) * ((t00 * gp.yr_to_sec)**(6/5)) #cm
    v00 = 1.2 * r00/(t00 * gp.yr_to_sec) #cm/s
    lt00 = (1-eta) * spin_down_lum(l0, t0, t00, nn)
    b00 = np.sqrt((gp.yr_to_sec*eta*6./(np.array(r00)**4)) * np.trapz(y=spin_down_lum(l0, t0, np.array([0,t00]), nn)*r00,x=np.array([0,t00])))
    emax00 = p_emax(eps, eta, l0, t0, t00, nn, b00)
    #Initial condition - First step - see APPENDIX B Gelfand et al. 2009
    r0 = 1.44 * ((((l0**2)*(e0**3))/(mej2**5))**0.1) * ((t[0] * gp.yr_to_sec)**(6/5))
    v0 = 1.2 * r0/(t[0] * gp.yr_to_sec)
    lt0 = (1-eta) * spin_down_lum(l0, t0, t[0], nn)
    b0 = np.sqrt((gp.yr_to_sec*eta*6./(np.array(r0)**4)) * np.trapz(y=spin_down_lum(l0, t0, np.array([0,t[0]]), nn)*r0,x=np.array([0,t[0]])))
    emax0 = p_emax(eps, eta, l0, t0, t[0], nn, b0)
    p_spectrum0 = broken_powerlaw(ebreak,alpha1,alpha2,[emax0],200)
    sed0, e0_p = particle_spectrum_start(t00, t[0], p_spectrum0, lt00, lt0, 0.0, 0.0, emax00, emax0, r00/gp.pc_to_cm, r0/gp.pc_to_cm, v0, v0, dens, Tfir, Ufir, Tnir, Unir) #Particle evolution
    p0_b = (b0**2)/(8*np.pi) #pressure from the magnetic field
    p0_p = (e0_p)/((4)*np.pi*(r0**3)) #pressure from the particles
    p0 = p0_b+p0_p #total pressure inside the PWN
    M_sw0 = (4*np.pi/3) * (r0**3) * (dens_ej(rr(r0/gp.pc_to_cm), R_snr[0], R_rs[0], t[0], e0, mej2, dens)[-1]) #Total mass in the thin shell outside the PWN
    f0 = 4*np.pi*(r0**2)*(p0) #force applied by the PWN on the mass shell (Pressure inside SNR = 0)

    #Storing variables
    EMAX = [emax0]
    R = [r0]
    V = [v00, v0]
    LT = [lt0]
    B = [np.float64(np.copy(b0))]
    M = [M_sw0]
    SED = [sed0]
    P = [p0]
    F = [f0]
    E = [e0_p]
    VEJ = [vel_ej(r0, t[0])]
    dens_ejecta = [dens_ej(rr(r0/gp.pc_to_cm), R_snr[0], R_rs[0], t[0], e0, mej2, dens)[-1]]
    
    with tqdm(total=len(t)-1, position=1, leave=True) as pbar:
        for i in tqdm(range(1, len(t)), position=1, leave=True, desc='PWN {0}'.format(pwn_n)):
            #print(' ')
            #print('######')
            #print(i+1, t[i])
            #print('######')
            #Start iteration
            r = r0 + (v0 * ((t[i]-t[i-1])*gp.yr_to_sec)) #Calculate new radius considering velocity at t-1
            R.append(np.copy(r))
            dens_ejecta.append(dens_ej(rr(r/gp.pc_to_cm), R_snr[i], R_rs[i], t[i], e0, mej2, dens)[-1]) #density of SNR ejecta at R_pwn
            LT.append((1-eta) * spin_down_lum(l0, t0, t[i], nn)) #spin down lum at time t[i]
            b = np.sqrt((gp.yr_to_sec*eta*6./(np.array(R)**4)) * np.concatenate(([0],  cumtrapz((LT*np.array(R)), x=t[:i+1])))) #mag field at t[i]
            B.append(np.float64(np.copy(b[-1])))
            EMAX.append(p_emax(eps, eta, l0, t0, t[i], nn, B[-1])) #maximum electron energy at t[i]
            p_spectrum0 = broken_powerlaw(ebreak,alpha1,alpha2,EMAX, 200)
            BB = np.copy(B)
            #Computation of particle spectrum with GAMERA (from t[0] to t[i])
            if i == 1: #This first step is needed like this since GAMERA seems to need at least three time steps to compute the evolution
                sed, e_p = particle_spectrum_start(t00, t[i], p_spectrum0, LT[-2], LT[-1], 0.0, 0.0, EMAX[-2], EMAX[-1], r0/gp.pc_to_cm, r/gp.pc_to_cm, V[-2], V[-1], dens, Tfir, Ufir, Tnir, Unir)
            else:
                #In order to obtain result in a reasonable amount of time the mag. field is "cutted" at a certain maximum value
                #The final SED compared to the one computed with the "uncutted" mag. have negligible differences
                if t[i] <= 5.:
                    #print('t<5yr')
                    mask = [BB>0.002] #max mag field 2000uG for t<5yr
                    BB[tuple(mask)] = 0.002
                    sed, e_p = particle_spectrum(t00, t[i], t, p_spectrum0, LT, BB, EMAX, np.copy(R)/gp.pc_to_cm, np.copy(V), dens, Tfir, Ufir, Tnir, Unir, no_escape)
                elif t[i] > 5. and t[i] <= 500.:
                    #print('t<500yr')
                    mask = [BB>0.0002] #max mag field 200uG for t<500yr
                    BB[tuple(mask)] = 0.0002
                    sed, e_p = particle_spectrum(t00, t[i], t, p_spectrum0, LT, BB, EMAX, np.copy(R)/gp.pc_to_cm, np.copy(V), dens, Tfir, Ufir, Tnir, Unir, no_escape)
                else:
                    #print('t>500yr')
                    sed, e_p = particle_spectrum(250., t[i], t, p_spectrum0, LT, BB, EMAX, np.copy(R)/gp.pc_to_cm, np.copy(V), dens, Tfir, Ufir, Tnir, Unir, no_escape) #tmin now is increased at 250 year to speed up the computation. To speed up more, increase this value
            #print('energy(erg), mag field(uG), mag field(uG, cut)',e_p, b[-1]*1e6, BB[-1]*1e6)
            p_b = (b[-1]**2)/(8*np.pi) #pressure from mag field
            p_p = (e_p)/(4*np.pi*(r**3)) #pressure from particles
            p = p_b+p_p #total pressure inside the PWN
            if r < R_rs[i]*gp.pc_to_cm: #computation valid only if R_PWN < R_rs for now
                VEJ.append(vel_ej(r, t[i]))
                if v0 >= vel_ej(r0, t[i-1]): #if expansion vel. of PWN is bigger then the velocity of ejecta inside the SNR the swept up mass increase
                    M_sw = M_sw0 + (((4/3)*np.pi) * ((r**3)-(r0**3)) * (dens_ej(rr(r/gp.pc_to_cm), R_snr[i], R_rs[i], t[i], e0, mej2, dens)[-1])) #swept up mass
                else:
                    M_sw = M_sw0
                f = 4*np.pi*(r**2)*(p) #force applied by the PWN on the mass shell (Pressure inside SNR = 0)
                v = (1/M_sw)*((M_sw0*(v0))+((M_sw-M_sw0)*vel_ej(r0, t[i-1]))+(f0*((t[i]-t[i-1])*gp.yr_to_sec))) #new velocity of the PWN
                #print('R_pwn, R_rs, R_snr (pc)',r/gp.pc_to_cm, R_rs[i], R_snr[i])
            if r >= R_rs[i]*gp.pc_to_cm: #End the computation in this case
                return np.array(R)/gp.pc_to_cm, np.array(V), np.array(VEJ), np.array(LT), np.array(B), np.array(M), np.array(EMAX), SED, np.array(P), np.array(F), np.array(E), np.array(dens_ejecta)
            #Storing variables
            V.append(np.copy(v))
            M.append(np.copy(M_sw))
            SED.append(np.copy(sed))
            P.append(np.copy(p))
            F.append(np.copy(f))
            E.append(np.copy(e_p))
            del r0, b0, v0, M_sw0, f0
            r0 = np.copy(r)
            b0 = np.copy(b)
            v0 = np.copy(v)
            M_sw0 = np.copy(M_sw)
            f0 = np.copy(f)
            #print('electron energy(TeV), Lum(erg/s)',(EMAX[-1]/gp.TeV_to_erg), LT[-1]/(1-eta))
            del r, b, v, M_sw, f
            pbar.update()
        pbar.close()
        return np.array(R)/gp.pc_to_cm, np.array(V), np.array(VEJ), np.array(LT), np.array(B), np.array(M), np.array(EMAX), SED, np.array(P), np.array(F), np.array(E), np.array(dens_ejecta)
