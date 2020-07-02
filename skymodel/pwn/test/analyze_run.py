import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import astropy.units as u
from gammapy.catalog import SourceCatalogGammaCat
from gammapy.utils.coordinates import galactic
from gammapy.modeling.models import create_crab_spectral_model as CrabSpectrum
from scipy.stats import kde
import aglib as ag
import sys, os
from astropy.table import Table



def import_and_transform(filename, energy_cut1=1*u.TeV, energy_cut2=10*u.PeV):
    dato = np.loadtxt(filename).T
    en = dato[0]*u.MeV
    diff_flux = (dato[1]/u.MeV/u.cm**2/u.s)
    integral = np.trapz(diff_flux[(en>=energy_cut1) & (en<=energy_cut2)], en[(en>=energy_cut1) & (en<=energy_cut2)])
    integral2 = np.trapz(diff_flux[(en>=energy_cut1) & (en<=energy_cut2)] * en[(en>=energy_cut1) & (en<=energy_cut2)], en[(en>=energy_cut1) & (en<=energy_cut2)])
    return integral.value, integral2.to('erg/(s*cm**2)').value

def plotLogNlogS(fl, binning, density=False, label=None, color='C0'):
    fl = np.nan_to_num(fl)
    logs_min = -10
    logs_max = 1
    nbins = binning*(logs_max - logs_min)
    bins_lognlogs = np.logspace(logs_min, logs_max, nbins)
    if label==None:
        n, bins, patches = plt.hist(fl, bins=bins_lognlogs, density=density, histtype='step', cumulative=-1, lw=2)
    else:
        n, bins, patches = plt.hist(fl, bins=bins_lognlogs, density=density, histtype='step', cumulative=-1, lw=2)
    f = np.sqrt(bins[1:] * bins[:-1])
    plt.fill_between(f, n - np.sqrt(n), n + np.sqrt(n), color=color, label=label,
                     alpha=0.3)

    plt.loglog()
    return n, bins

def flux_from_gammacat(cat,emin=1,emax=1000, lattresh=2., lowlontresh=70., highlontresh=270.):
    # calculate integral flux in desired energy range from spectral model
    fluxes = np.array([])
    for source in cat:
        try:
            if np.abs(source.spatial_model().lat_0.value) <= lattresh:
                if (source.spatial_model().lon_0.value) <= lowlontresh or (source.spatial_model().lon_0.value) >= highlontresh:
                    try:
                        flux = source.spectral_model().integral(emin*u.TeV,emax*u.TeV)
                        fluxes = np.append(fluxes,flux.value)
                    except:
                        # sources without spectral model
                        fluxes = np.append(fluxes, np.nan)
                else:
                    fluxes = np.append(fluxes, np.nan)
            else:
                fluxes = np.append(fluxes, np.nan)
        except:
            fluxes = np.append(fluxes, np.nan)
        
    crab = CrabSpectrum('meyer')
    crab_flux = crab.integral(emin*u.TeV, emax*u.TeV).value
    crab_flux_1TeV = crab.integral(1*u.TeV, 1000.*u.TeV).value
    fluxes /= crab_flux
    return fluxes

def plot_figure_singlesrc(samples, fluxes_gammacat, fluxes_synt_pop, title, xlabel, figname):
    plt.figure()
    for s, sample in enumerate(samples[:]):
        mask = np.zeros(len(gammacat.table),dtype=bool)
        for c in sample['classes']:
            mask = np.logical_or(mask,gammacat.table['classes'] == c)
        # select sample
        flux_sample = fluxes_gammacat[mask==True]
        aa = plotLogNlogS(flux_sample, 10, label=sample['name'], color=color_cycle[s])
    aaa = plotLogNlogS(fluxes_synt_pop, 10, color='C3', label='Synthetic population')
    plt.legend()
    plt.xlim(1e-4, 0.3e1)
    plt.ylim(0.39806199042692636, 1000)
    plt.grid()
    plt.title(title)
    plt.xlabel(xlabel, fontsize=11)
    plt.ylabel('Number of sources (> Flux)', fontsize=11)
    plt.tight_layout()
    plt.savefig(figname, dpi=200, bbox_inches='tight')
    plt.close()
    
def plot_figure_allsrcs(fluxes_gammacat, fluxes_synt_pop1, fluxes_synt_pop2, fluxes_synt_pop3, title, xlabel, figname):
    plt.figure()
    aa = plotLogNlogS(fluxes_gammacat, 10, label='Gammacat sources', color='C0')
    aaa = plotLogNlogS(np.concatenate([fluxes_synt_pop1,fluxes_synt_pop2,fluxes_synt_pop3]), 10, color='C1', label='Synthetic populations')
    plt.legend()
    plt.xlim(1e-4, 0.3e1)
    plt.ylim(0.39806199042692636, 1000)
    plt.grid()
    plt.title(title)
    plt.xlabel(xlabel, fontsize=11)
    plt.ylabel('Number of sources (> Flux)', fontsize=11)
    plt.tight_layout()
    plt.savefig(figname, dpi=200, bbox_inches='tight')
    plt.close()
    
def plot_figure_pluscrab1(fluxes_gammacat, fluxes_pop, title, xlabel, filename):
    plt.figure()
    for s, sample in enumerate(samples[:]):
        mask = np.zeros(len(gammacat.table),dtype=bool)
        for c in sample['classes']:
            mask = np.logical_or(mask,gammacat.table['classes'] == c)
        # select sample
        flux_sample = fluxes_gammacat[mask==True]
        aa = plotLogNlogS(np.concatenate([flux_sample, np.array([0.6823, 0.6333])]), 10, label=sample['name'], color=color_cycle[s])
    aaa = plotLogNlogS(fluxes_pop, 10, color='C3', label='Synthetic population')
    plt.legend()
    plt.xlim(1e-4, 0.3e1)
    plt.ylim(0.39806199042692636, 1000)
    plt.grid()
    plt.title(title, fontsize=11)
    plt.xlabel(xlabel, fontsize=11)
    plt.ylabel('Number of sources (> Flux)', fontsize=11)
    plt.tight_layout()
    plt.savefig(filename, dpi=200, bbox_inches='tight')
    plt.close()
    
def plot_figure_pluscrab01(fluxes_gammacat, fluxes_pop, title, xlabel, filename):
    plt.figure()
    for s, sample in enumerate(samples[:]):
        mask = np.zeros(len(gammacat.table),dtype=bool)
        for c in sample['classes']:
            mask = np.logical_or(mask,gammacat.table['classes'] == c)
        # select sample
        flux_sample = fluxes_gammacat[mask==True]
        aa = plotLogNlogS(np.concatenate([flux_sample, np.array([0.7201, 0.077])]), 10, label=sample['name'], color=color_cycle[s])
    aaa = plotLogNlogS(fluxes_pop, 10, color='C3', label='Synthetic population')
    plt.legend()
    plt.xlim(1e-4, 0.3e1)
    plt.ylim(0.39806199042692636, 1000)
    plt.grid()
    plt.title(title, fontsize=11)
    plt.xlabel(xlabel, fontsize=11)
    plt.ylabel('Number of sources (> Flux)', fontsize=11)
    plt.tight_layout()
    plt.savefig(filename, dpi=200, bbox_inches='tight')
    plt.close()

col_names = ['N','filename', 'GLON', 'GLAT', 'R_pwn', 'R_pwn_deg', 'distance',
             'cr_fl_01', 'cr_fl_1', 'nh', 'age', 'X', 'Y', 'Z', 'v_X', 'v_Y', 'v_Z',
             'v_3d', 'E_sn', 'M_ej', 'l0', 't0', 'eta', 'eps', 'ebreak', 'alpha1',
             'alpha2', 'Tfir', 'Ufir', 'Tnir', 'Unir']

samples = [ {'name' : 'PWNe', 'classes' : ['pwn']},
            {'name' : 'PWNe + composites', 'classes' : ['pwn','pwn,snr']},
            {'name' : 'PWNe + composites + UNID', 'classes' : ['pwn','pwn,snr','unid']}
            ]

color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

final = pd.read_csv('../PWNe_final_population.txt', delimiter='\t', header=0, index_col=0, usecols=range(len(col_names)), names=col_names)
final['filename'] = '../xml/ph_sed_pwn_'+final.index.astype(str)+'.txt'

ph_fl_01 = []
for i in range(len(final)):
    a, b = import_and_transform(final.filename.iloc[i], energy_cut1=0.1*u.TeV, energy_cut2=1000*u.TeV)
    ph_fl_01.append(a)
    
ph_fl_1 = []
for i in range(len(final)):
    a, b = import_and_transform(final.filename.iloc[i], energy_cut1=1*u.TeV, energy_cut2=1000*u.TeV)
    ph_fl_1.append(a)
    
crab = CrabSpectrum('meyer')
emin01, emin1, emax = [0.1, 1, 1000] * u.TeV

crab_01 = crab.integral(emin01, emax).value
crab_1 = crab.integral(emin1, emax).value

flux_int_cu = (ph_fl_01 / crab_01)
flux_int_cu01 = (ph_fl_1 / crab_1)

final['cr_fl_1'] = flux_int_cu
final['cr_fl_01'] = flux_int_cu01
    
final2 = final[(final.GLAT<=2.) & (final.GLAT>=-2.) & (final.GLON<=130.) & (final.GLON>=-70)]
final3 = final[(final.GLAT<=2.) & (final.GLAT>=-2.)]

gammacat_file = '../../known-sources/external-input/gammacat.fits.gz'

gammacat = SourceCatalogGammaCat(gammacat_file)

gammacat_pwn_glon, gammacat_pwn_glat = [], []
for source in gammacat:
    if source.data.where.startswith('gal'):
        if 'pwn' or 'unid' in  source.data.classes:
            try:
                gammacat_pwn_glon.append(source.spatial_model().lon_0.value)
                gammacat_pwn_glat.append(source.spatial_model().lat_0.value)
            except:
                None

gammacat_pwn_glat  = np.array(gammacat_pwn_glat)
gammacat_pwn_glon  = np.array(gammacat_pwn_glon)
gammacat_pwn_glon = np.concatenate([gammacat_pwn_glon[gammacat_pwn_glon>180] - 360, gammacat_pwn_glon[gammacat_pwn_glon<180]])

k = kde.gaussian_kde(np.array([gammacat_pwn_glon, gammacat_pwn_glat]))
nbins=200
xi, yi = np.mgrid[gammacat_pwn_glon.min():gammacat_pwn_glon.max():nbins*1j, gammacat_pwn_glat.min():gammacat_pwn_glat.max():nbins*1j]
zi = k(np.vstack([xi.flatten(), yi.flatten()]))
zi /= zi.max()

glat = final.GLAT
glon = final.GLON
glon = np.concatenate([glon[glon>180] - 360, glon[glon<180]])

k1 = kde.gaussian_kde(np.array([glon, glat]))
nbins=200
xi1, yi1 = np.mgrid[glon.min():glon.max():nbins*1j, glat.min():glat.max():nbins*1j]
zi1 = k1(np.vstack([xi1.flatten(), yi1.flatten()]))
zi1 /= zi1.max()

fig1 = plt.figure()
plt.scatter(gammacat_pwn_glon, gammacat_pwn_glat, c='C3', alpha=0.33)
CS =plt.contour(xi, yi, zi.reshape(xi.shape), np.array([0.05, 0.2, 0.8]), cmap=plt.cm.jet_r)
CS.levels = [(1-val) for val in CS.levels]
plt.clabel(CS, CS.levels, inline=True, fmt=f'%.2f',  fontsize=10)
plt.gca().invert_xaxis()
plt.axhline(-4)
plt.axhline(3)
plt.xlabel('GLON')
plt.ylabel('GLAT')
plt.xlim(130, -120)
plt.ylim(-6,6)
plt.gca().set_xticks([100,50,0,-50, -100])
plt.gca().set_yticks([-6,-4, 0, 3, 6])
plt.gca().set_xticklabels(['100', '50' , '0', '350', '300'])
plt.title('Real Source Distribution (PWNe+Composite+UNID)')
plt.tight_layout()
fig1.savefig('real_source_distr.png', dpi=150, bbox_inches='tight')
plt.close()

fig2 = plt.figure()
CS =plt.contour(xi, yi, zi.reshape(xi.shape), np.array([0.05, 0.2, 0.8]), cmap=plt.cm.jet)
CS.levels = [(1-val) for val in CS.levels]
plt.clabel(CS, CS.levels, inline=True, fmt=f'%.2f',  fontsize=10)
plt.axvline(380, label='Real Sources', c='k', lw=1.33)
plt.axvline(380, label='Sim. Sources', c='k', lw=1.33, ls='--')
CS =plt.contour(xi1, yi1, zi1.reshape(xi1.shape), np.array([0.05, 0.2, 0.8]), cmap=plt.cm.jet, alpha=0.5, linestyles='--')
CS.levels = [(1-val) for val in CS.levels]
plt.clabel(CS, CS.levels, inline=True, fmt=f'%.2f',  fontsize=11)
plt.legend()
plt.gca().invert_xaxis()
plt.xlim(120, -120)
plt.ylim(-10, 10)
plt.xlabel('GLON')
plt.ylabel('GLAT')
plt.xlim(130, -120)
plt.gca().set_xticks([100,50,0,-50, -100])
plt.gca().set_xticklabels(['100', '50' , '0', '350', '300'])
plt.title('Real Source Distribution (PWNe+Composite+UNID) VS Simulated')
plt.tight_layout()
fig2.savefig('real_source_distr_vs_sim.png', dpi=150, bbox_inches='tight')
plt.close()

fluxes1_cut = flux_from_gammacat(gammacat, emin=1,emax=1000, lattresh=2., lowlontresh=70., highlontresh=270.)
fluxes01_cut = flux_from_gammacat(gammacat, emin=0.1,emax=1000, lattresh=2., lowlontresh=70., highlontresh=270.)

fluxes1_cut3 = flux_from_gammacat(gammacat, emin=1,emax=1000, lattresh=2., lowlontresh=27000., highlontresh=-270000.)
fluxes01_cut3 = flux_from_gammacat(gammacat, emin=0.1,emax=1000, lattresh=2., lowlontresh=27000., highlontresh=-27000.)

fluxes1 = flux_from_gammacat(gammacat, emin=1,emax=1000, lattresh=200., lowlontresh=1000., highlontresh=-1000.)
fluxes01 = flux_from_gammacat(gammacat, emin=1,emax=1000, lattresh=200., lowlontresh=1000., highlontresh=-1000.)


snrsfile  = '../../snr/ALL_FILES_0/results_0.txt'
isnrsfile = '../../int-snr/out/d.fits'

synt = Table()

# In[SNR]  
  
q = Table.read(snrsfile, format='ascii')
en=q['E[TeV]'][0:40].data               # TeV
glon, glat, f01, f1 = [], [], [], []

for i in np.arange(int(len(q)/40)):

    sed=q['diff_spectrum'][i*40:(i+1)*40]   #  TeV cm-2 s-1  
    fl=sed / en**2.                         #  cm-2 s-1 TeV-1 
         
    if (sum(fl) != 0):
        pos = galactic(x=q['POS_Y'][i*40]*u.kpc, y=-q['POS_X'][i*40]*u.kpc, z=q['POS_Z'][i*40]*u.kpc)
        glon.append(pos[1].value)
        glat.append(pos[2].value)
        f01.append(ag.integ(en,fl, emin=0.1 ,emax=1000.)[0]/crab_01)
        f1.append(ag.integ(en,fl, emin=1.0 ,emax=1000.)[0]/crab_1)

dic = {'GLON' : glon,
       'GLAT': glat,
       'F01': f01,
       'F1': f1}

snr = pd.DataFrame(dic)

# In[iSNR]

isnrs = Table.read(isnrsfile)
f01, f1 = [], []
for i in np.arange(len(isnrs)) :
    isnr=isnrs[i]   
  
    sp=Table.read('../../int-snr/out/isnr'+str(i)+'_spec.txt' ,format='ascii')
    en = sp['col1']   #  MeV
    f  = sp['col2']   #  ph cm-2 s-1 MeV-1 
    f01.append(ag.integ(en,f,emin=1e5,emax=1e9)[0]/crab_01)
    f1.append(ag.integ(en,f,emin=1e6,emax=1e9)[0]/crab_1)



dic1 = {'GLON' : np.array(isnrs['Glon'], dtype=float),
       'GLAT': np.array(isnrs['Glat'], dtype=float),
       'F01': f01,
       'F1': f1}

isr = pd.DataFrame(dic1)

snr2 = snr[(snr.GLAT<=2.) & (snr.GLAT>=-2.) & (snr.GLON<=130.) & (snr.GLON>=-70)]
isr2 = isr[(isr.GLAT<=2.) & (isr.GLAT>=-2.) & (isr.GLON<=130.) & (isr.GLON>=-70)]

snr3 = snr[(snr.GLAT<=2.) & (snr.GLAT>=-2.)]
isr3 = isr[(isr.GLAT<=2.) & (isr.GLAT>=-2.)]

samples2 = [ {'name' : 'SNR', 'classes' : ['snr']},
            {'name' : 'SNR + composites', 'classes' : ['snr','pwn,snr']},
            {'name' : 'SNR + composites + other', 'classes' : ['snr','pwn,snr','snr,mc', 'unid,snr,mc']}
            ]

plot_figure_singlesrc(samples, fluxes1_cut, final2.cr_fl_1,   'Sources in |GLAT|<2 and (GLON<70 | GLON>270)', "Flux > 1.0 TeV (Crab units)", 'PWNe_logNlogS_1TeV_HGPS_region.png') 
plot_figure_singlesrc(samples, fluxes01_cut, final2.cr_fl_01, 'Sources in |GLAT|<2 and (GLON<70 | GLON>270)', "Flux > 0.1 TeV (Crab units)", 'PWNe_logNlogS_0.1TeV_HGPS_region.png') 
plot_figure_singlesrc(samples, fluxes1_cut3, final3.cr_fl_1,   'Sources in |GLAT|<2', "Flux > 1.0 TeV (Crab units)", 'PWNe_logNlogS_1TeV_latcut.png') 
plot_figure_singlesrc(samples, fluxes01_cut3, final3.cr_fl_01, 'Sources in |GLAT|<2', "Flux > 0.1 TeV (Crab units)", 'PWNe_logNlogS_0.1TeV_latcut.png') 
plot_figure_singlesrc(samples, fluxes1, final.cr_fl_1,   'All PWNe', "Flux > 1.0 TeV (Crab units)", 'PWNe_logNlogS_1TeV_nocut.png') 
plot_figure_singlesrc(samples, fluxes01, final.cr_fl_01, 'All PWNe', "Flux > 0.1 TeV (Crab units)", 'PWNe_logNlogS_0.1TeV_nocut.png') 

plot_figure_singlesrc(samples2, fluxes1_cut, np.concatenate([snr2.F1[snr2.F1<3.], isr2.F1]),   'Sources in |GLAT|<2 and (GLON<70 | GLON>270)', "Flux > 1.0 TeV (Crab units)", 'SNR_logNlogS_1TeV_HGPS_region.png') 
plot_figure_singlesrc(samples2, fluxes01_cut, np.concatenate([snr2.F01[snr2.F01<3.], isr2.F01]), 'Sources in |GLAT|<2 and (GLON<70 | GLON>270)', "Flux > 0.1 TeV (Crab units)", 'SNR_logNlogS_0.1TeV_HGPS_region.png') 
plot_figure_singlesrc(samples2, fluxes1_cut3, np.concatenate([snr3.F1[snr3.F1<3.], isr3.F1]),   'Sources in |GLAT|<2', "Flux > 1.0 TeV (Crab units)", 'SNR_logNlogS_1TeV_latcut.png') 
plot_figure_singlesrc(samples2, fluxes01_cut3,np.concatenate([snr3.F01[snr3.F01<3.], isr3.F01]), 'Sources in |GLAT|<2', "Flux > 0.1 TeV (Crab units)", 'SNR_logNlogS_0.1TeV_latcut.png') 
plot_figure_singlesrc(samples2, fluxes1, np.concatenate([snr.F1[snr.F1<3.], isr.F1]),   'All SNRs', "Flux > 1.0 TeV (Crab units)", 'SNR_logNlogS_1TeV_nocut.png') 
plot_figure_singlesrc(samples2, fluxes01, np.concatenate([snr.F01[snr.F01<3.], isr.F01]), 'All SNRs', "Flux > 0.1 TeV (Crab units)", 'SNR_logNlogS_0.1TeV_nocut.png') 

plot_figure_allsrcs(fluxes1_cut, final2.cr_fl_1, snr2.F1[snr2.F1<3], isr2.F1, 'Sources in |GLAT|<2 and (GLON<70 | GLON>270)', "Flux > 1.0 TeV (Crab units)", 'All_logNlogS_1TeV_HGPS_region.png')
plot_figure_allsrcs(fluxes01_cut, final2.cr_fl_01, snr2.F01[snr2.F01<3], isr2.F01, 'Sources in |GLAT|<2 and (GLON<70 | GLON>270)', "Flux > 0.1 TeV (Crab units)", 'All_logNlogS_0.1TeV_HGPS_region.png')
plot_figure_allsrcs(fluxes1_cut3, final3.cr_fl_1, snr3.F1[snr3.F1<3], isr3.F1, 'Sources in |GLAT|<2', "Flux > 1.0 TeV (Crab units)", 'All_logNlogS_1TeV_latcut.png')
plot_figure_allsrcs(fluxes01_cut3, final3.cr_fl_01, snr3.F01[snr3.F01<3], isr3.F01, 'Sources in |GLAT|<2', "Flux > 0.1 TeV (Crab units)", 'All_logNlogS_0.1TeV_latcut.png')
plot_figure_allsrcs(fluxes1, final.cr_fl_1, snr.F1[snr.F1<3], isr.F1, 'All Sources', "Flux > 1.0 TeV (Crab units)", 'All_logNlogS_1TeV_nocut.png')
plot_figure_allsrcs(fluxes01, final.cr_fl_01, snr.F01[snr.F01<3], isr.F01, 'All Sources', "Flux > 0.1 TeV (Crab units)", 'All_logNlogS_0.1TeV_nocut.png')

plot_figure_pluscrab1(fluxes1_cut, final2.cr_fl_1, 'Sources in |GLAT|<2 and (GLON<70 | GLON>270) + Crab&others', "Flux > 1.0 TeV (Crab units)", 'PWNe_logNlogS_1TeV_HGPS_region_plus_outside_luminous_sources.png')
plot_figure_pluscrab01(fluxes01_cut, final2.cr_fl_01, 'Sources in |GLAT|<2 and (GLON<70 | GLON>270) + Crab&others', "Flux > 0.1 TeV (Crab units)", 'PWNe_logNlogS_0.1TeV_HGPS_region_plus_outside_luminous_sources.png')
plot_figure_pluscrab1(fluxes1_cut3, final3.cr_fl_1, 'Sources in |GLAT|<2 + Crab&others', "Flux > 1.0 TeV (Crab units)", 'PWNe_logNlogS_1TeV_latcut_plus_outside_luminous_sources.png')
plot_figure_pluscrab01(fluxes01_cut3, final3.cr_fl_01, 'Sources in |GLAT|<2 + Crab&others', "Flux > 0.1TeV (Crab units)", 'PWNe_logNlogS_0.1TeV_latcut_plus_outside_luminous_sources.png')

plt.figure()
aa = plotLogNlogS(np.concatenate([fluxes1_cut, np.array([0.6823, 0.7038, 0.6333])]), 10, label='Gammacat sources', color='C0')
aaa = plotLogNlogS(np.concatenate([final2.cr_fl_1, snr2.F1[snr2.F1<3], isr2.F1]), 10, color='C1', label='Synthetic populations')
plt.legend()
plt.xlim(1e-4, 0.3e1)
plt.ylim(0.39806199042692636, 1000)
plt.grid()
plt.title('Sources in |GLAT|<2 and (GLON<70 | GLON>270) + Crab&others')
plt.xlabel("Flux > 1.0 TeV (Crab units)", fontsize=11)
plt.ylabel('Number of sources (> Flux)', fontsize=11)
plt.tight_layout()
plt.savefig("All_logNlogS_1TeV_HGPS_region_plus_outside_luminous_sources.png", dpi=200, bbox_inches='tight')
plt.close()

plt.figure()
aa = plotLogNlogS(np.concatenate([fluxes01_cut, np.array([0.7201, 0.257, 0.077])]), 10, label='Gammacat sources', color='C0')
aaa = plotLogNlogS(np.concatenate([final2.cr_fl_01, snr2.F01[snr2.F01<3], isr2.F01]), 10, color='C1', label='Synthetic populations')
plt.legend()
plt.xlim(1e-4, 0.3e1)
plt.ylim(0.39806199042692636, 1000)
plt.grid()
plt.title('Sources in |GLAT|<2 and (GLON<70 | GLON>270) + Crab&others ')
plt.xlabel("Flux > 0.1 TeV (Crab units)", fontsize=11)
plt.ylabel('Number of sources (> Flux)', fontsize=11)
plt.tight_layout()
plt.savefig("All_logNlogS_0.1TeV_HGPS_region_plus_outside_luminous_sources.png", dpi=200, bbox_inches='tight')
plt.close()

plt.figure()
aa = plotLogNlogS(np.concatenate([fluxes1_cut3, np.array([0.6823, 0.6333])]), 10, label='Gammacat sources', color='C0')
aaa = plotLogNlogS(np.concatenate([final3.cr_fl_1, snr3.F1[snr3.F1<3], isr3.F1]), 10, color='C1', label='Synthetic populations')
plt.legend()
plt.xlim(1e-4, 0.3e1)
plt.ylim(0.39806199042692636, 1000)
plt.grid()
plt.title('Sources in |GLAT|<2 + Crab&others')
plt.xlabel("Flux > 1.0 TeV (Crab units)", fontsize=11)
plt.ylabel('Number of sources (> Flux)', fontsize=11)
plt.tight_layout()
plt.savefig("All_logNlogS_1TeV_latcut_plus_outside_luminous_sources.png", dpi=200, bbox_inches='tight')
plt.close()

plt.figure()
aa = plotLogNlogS(np.concatenate([fluxes01_cut3, np.array([0.7201, 0.077])]), 10, label='Gammacat sources', color='C0')
aaa = plotLogNlogS(np.concatenate([final3.cr_fl_01, snr3.F01[snr3.F01<3], isr3.F01]), 10, color='C1', label='Synthetic populations')
plt.legend()
plt.xlim(1e-4, 0.3e1)
plt.ylim(0.39806199042692636, 1000)
plt.grid()
plt.title('Sources in |GLAT|<2 + Crab&others ')
plt.xlabel("Flux > 0.1 TeV (Crab units)", fontsize=11)
plt.ylabel('Number of sources (> Flux)', fontsize=11)
plt.tight_layout()
plt.savefig("All_logNlogS_0.1TeV_latcut_plus_outside_luminous_sources.png", dpi=200, bbox_inches='tight')
plt.close()
