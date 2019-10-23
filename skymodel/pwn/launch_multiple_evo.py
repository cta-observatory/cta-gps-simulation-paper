import astropy.units as u
from astropy.constants import m_e, c, e, m_p
import sys
sys.path.append('/home/fiori/GAMERA/lib') #Da modificare in base a dove è installato GAMERA!
import gappa as gp
import numpy as np
from scipy.integrate import cumtrapz
from tqdm import tqdm as tqdm
import pandas as pd
from evo_gelfand_model_multiple import *
from astropy.coordinates import Galactocentric, ICRS

if len(sys.argv) > 1: #Specify the input file with all the SNRs
    snr_file = sys.argv[1]
else:
    snr_file = str(input('Input file: '))

def f_prop(dist_x): #probability density function to find a pulsar in the log(l0*)-log(tau0)* plane inside SNR (from Faucher-Gigue and Kaspi 2006)
    dist_x, dist_y = dist_x[:, 0], dist_x[:, 1]
    return np.exp((-5.9533*(dist_x+1.0324)**2) - (6.4263 * (dist_x + 1.0324) * (dist_y + 1.2977)) - (2.1164 * (dist_y + 1.2977)**2))

def metropolis_hastings(target_density, size=1): #Metropoli-hastings algorithm to compute the value of l0* and tau0*
    burnin_size = 1000
    x0 = np.array([[0, 0]])
    xt = x0
    b = 0
    a = 0
    samples = []
    while b < burnin_size:
        xt_candidate = np.array([np.random.multivariate_normal(xt[0], np.eye(2))])
        accept_prob = (target_density(xt_candidate))/(target_density(xt))
        if np.random.uniform(0, 1) < accept_prob:
                xt = xt_candidate
                b += 1
    while a < size:
    #for i in (range(size)):
        xt_candidate = np.array([np.random.multivariate_normal(xt[0], np.eye(2))])
        accept_prob = (target_density(xt_candidate))/(target_density(xt))
        if np.random.uniform(0, 1) < accept_prob:
                xt = xt_candidate
                a += 1
                #print(a)
                samples.append(xt)
    samples = np.array(samples[:])
    samples = np.reshape(samples, [samples.shape[0], 2])
    return samples

def skip(index): #simple function to extract one row for each SNR for the snr_file
    if (index-1) % 40 == 0:
        return False
    return True

def distance(r, theta, z): #distance of the system
    return np.sqrt((8.3-(r*np.cos(theta)))**2+((r**2)*(np.sin(theta)**2))+(z-0.027)**2)

def extract_SNR_and_pulsar_par(snr_file): #extract parameters of the SNRs from the relative file and compute all the values need for the computation of the PWNe  
    column_name = ['N', 'r', 'theta', 'z', 'nh', 'type', 'age']
    data_snr = pd.read_csv(snr_file, skiprows= lambda x: skip(x), delimiter='\t', usecols=[0, 1, 2, 3, 4, 5, 6], header=None, index_col=0, names=column_name)
    conv = u.year.to('s')
    data_snr.age = data_snr.age *1e3
    data_snr.loc[data_snr['type'] == 1, 'E_sn'] = 1e51
    data_snr.loc[data_snr['type'] == 4, 'E_sn'] = 3e51
    data_snr.loc[data_snr['type'] == 2, 'E_sn'] = 1e51
    data_snr.loc[data_snr['type'] == 3, 'E_sn'] = 1e51
    data_snr.loc[data_snr['type'] == 1, 'M_ej'] = ((1.4*u.M_sun).to('g')).value
    data_snr.loc[data_snr['type'] == 4, 'M_ej'] = ((1.*u.M_sun).to('g')).value
    data_snr.loc[data_snr['type'] == 2, 'M_ej'] = ((8.*u.M_sun).to('g')).value
    data_snr.loc[data_snr['type'] == 3, 'M_ej'] = ((2.*u.M_sun).to('g')).value
    data_snr['l_ch'] = data_snr['E_sn']**(3/2) * data_snr['M_ej']**(-5/6) * (data_snr['nh']*(m_p.to('g').value))**(1/3)
    data_snr['t_ch'] = data_snr['E_sn']**(-1/2) * data_snr['M_ej']**(5/6) * (data_snr['nh']*(m_p.to('g').value))**(-1/3) / conv
    samples = metropolis_hastings(f_prop, size=len(data_snr))
    data_snr['t*'] = samples[:,0] 
    data_snr['l*'] = samples[:,1]
    data_snr['prob'] = f_prop(samples)
    data_snr['l0'] = (10**(data_snr['l*'])) * data_snr['l_ch']
    data_snr['t0'] = (10**(data_snr['t*'])) * data_snr['t_ch']
    data_snr['distance'] = distance(data_snr['r'].values, data_snr['theta'].values, data_snr['z'].values)
    
    ##RANDOMLY GENERATION OF (UNKNOWN) PWN PARAMETERS
    data_snr['eta'] = np.random.uniform(0.01,0.4, data_snr.shape[0]) #Magnetization of the PWN
    data_snr['eps'] = np.random.uniform(0.066,0.66, data_snr.shape[0]) #Confinment factor
    data_snr['ebreak'] = np.random.uniform(0.01, 10, data_snr.shape[0])*gp.TeV_to_erg #energy break of the injection spectrum (Broken PWL)
    data_snr['alpha1'] = np.random.uniform(1.0, 1.7, data_snr.shape[0]) #first index of the injection spectrum (Broken PWL)
    data_snr['alpha2'] = np.random.uniform(2.0, 2.7, data_snr.shape[0]) #second index of the injection spectrum (Broken PWL)

    #Photon Background Info
    data_snr['Tfir'] = np.random.uniform(1, 50, data_snr.shape[0]) #Temperature of FIR field in K --> NEED TO INSERT A VALUE STARTING FROM GALPROP
    data_snr['Ufir'] = np.random.uniform(0.1, 4, data_snr.shape[0])*gp.eV_to_erg #density of FIR field in erg/cm^3 --> NEED TO INSERT A VALUE STARTING FROM GALPROP
    data_snr['Tnir'] = np.random.uniform(2000, 5000, data_snr.shape[0]) #Temperature of NIR field in K --> NEED TO INSERT A VALUE STARTING FROM GALPROP
    data_snr['Unir'] = np.random.uniform(0.5, 20, data_snr.shape[0])*gp.eV_to_erg #density of NIR field in erg/cm^3 --> NEED TO INSERT A VALUE STARTING FROM GALPROP
    
    #Initialize the column with the radii of the SNR (forward shock and reverse shock) and the PWN
    data_snr['R_fs'] = np.zeros(data_snr.shape[0]) 
    data_snr['R_rs'] = np.zeros(data_snr.shape[0])
    data_snr['R_pwn'] = np.zeros(data_snr.shape[0])
    
    del data_snr['type']
    del data_snr['t*']
    del data_snr['l*']
    del data_snr['l_ch']
    del data_snr['t_ch']
    pwn = data_snr['prob']>0.2
    return data_snr.loc[pwn]

def execute_multiple_evolutions(systems): #Compute the evolution of the PWNe, saving one file for every PWN and one file with the summary of the parameters of the PWNe 
    systems.to_csv('PWN_initial_parameters.txt', header=True, sep='\t', encoding='utf-8')
    for index, row in tqdm(systems.iterrows(), total=len(systems)):
        print('System n: {0}'.format(index))
        R_fs, R_rs, R_pwn = launch_evo_m(index, row.age, row.l0, row.E_sn, row.M_ej, 3, row.eta, row.t0, row.eps, row.nh, row.Tfir, row.Ufir, row.Tnir, row.Unir, 50, 0.1*row.age, row.ebreak, row.alpha1, row.alpha2, row.distance)
        if R_rs[len(R_pwn)-1] < R_pwn[-1]:
            print('PWN bigger than the SNR reverse shock!')
        row.R_fs = R_fs[-1]
        row.R_rs = R_rs[-1]
        row.R_pwn = R_pwn[-1]
        del R_fs, R_rs, R_pwn

    for index, row in systems.iterrows():
        if row.R_rs < row.R_pwn:
            systems = systems.drop(index=index)       
    systems.to_csv('PWN_results.txt', header=True, sep='\t', encoding='utf-8')
    
    
systems = extract_SNR_and_pulsar_par(snr_file)
execute_multiple_evolutions(systems)