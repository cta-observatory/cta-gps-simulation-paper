import matplotlib.pyplot as plt
import gammalib
from gammapy.catalog import SourceCatalogGammaCat
import astropy.units as u
import numpy as np
import sys
sys.path.append('../../scripts/')
import utils
import pdb

############## input and configuration ###################

# data files
# PWN population
pwnpop_file = '../xml/pwn.xml'
# gamma-cat
gammacat_file = '../../known-sources/external-input/gammacat.fits.gz'

# configuration
# samples of real sources to compare with
samples = [ {'name' : 'PWNe', 'classes' : ['pwn']},
            {'name' : 'PWNe + composites', 'classes' : ['pwn','pwn,snr']},
            {'name' : 'PWNe + composites + UNID', 'classes' : ['pwn','pwn,snr','unid']}
            ]


# energy thresholds (TeV)
thresholds = [0.1,1.]

# binning for logN-logS
bins_lognlogs = np.logspace(-5, 1., 60)

# default maximum energy for spectral integration (TeV)
eup = 1000.

############################################################

def make_lognlogs(flux,bins = None, color = 'k',label=None):
    # histogram
    n, bins, patches = plt.hist(flux, bins=bins, color = color,
                                alpha=0.8, density=False, histtype='step', cumulative=-1)
    # geometric average of flux in bin
    f = np.sqrt(bins[1:] * bins[:-1])
    # shaded band with Poisson uncertainties
    plt.fill_between(f, n - np.sqrt(n), n + np.sqrt(n), color=color, label=label,
                     alpha=0.3)

def flux_from_gammacat(cat,emin,emax=eup):
    # calculate integral flux in desired energy range from spectral model
    fluxes = np.array([])
    for source in cat:
        try:
            flux = source.spectral_model().integral(emin*u.TeV,emax*u.TeV)
            fluxes = np.append(fluxes,flux.value)
        except:
            # sources without spectral model
            fluxes = np.append(fluxes, np.nan)

    # convert to Crab units
    # conversion consistent with utils for gammalib model containers
    # simple Crab TeV power law model
    crab = gammalib.GModelSpectralPlaw(5.7e-16, -2.48, gammalib.GEnergy(0.3, 'TeV'))
    # calculate crab flux over the desired energy range
    crab_flux = crab.flux(gammalib.GEnergy(emin,'TeV'), gammalib.GEnergy(emax,'TeV'))
    # remormalize crab flux so that it matches the Meyer model > 1 TeV (consistent with gamma-cat)
    crab_flux_1TeV = crab.flux(gammalib.GEnergy(1.,'TeV'), gammalib.GEnergy(eup,'TeV'))
    # (Meyer model, flux > 1 TeV in ph cm-2 s-1)
    crab_flux *= 2.0744340476909142e-11 / crab_flux_1TeV
    fluxes /= crab_flux

    return fluxes

# color cycle
color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

# open gammacat using gammapy
gammacat = SourceCatalogGammaCat(gammacat_file)

# read PWN population as gammalib model container
pwnpop = gammalib.GModels(pwnpop_file)

# loop over thresholds
for thresh in thresholds:
    # create figure
    fig1 = plt.figure()
    ax1 = plt.subplot()
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel("Flux > {} TeV (Crab units)".format(thresh), fontsize=14)
    ax1.set_ylabel('Number of sources (> Flux)', fontsize=14)
    utils.format_ax(ax1)

    # synthetic population
    # extract population fluxes
    lons, lats, radii, fluxes, names = utils.dist_from_gammalib(pwnpop,emin=thresh,emax=eup)
    # plots
    make_lognlogs(fluxes,bins=bins_lognlogs, color = 'k', label = 'synthetic population')

    # compare to gamma-cat samples
    # extract fluxes from gammacat
    fluxes = flux_from_gammacat(gammacat, emin=thresh, emax=eup)
    # samples
    for s, sample in enumerate(samples):
        # build mask to select desired sample
        mask = np.zeros(len(gammacat.table),dtype=bool)
        for c in sample['classes']:
            mask = np.logical_or(mask,gammacat.table['classes'] == c)
        # select sample
        flux_sample = fluxes[mask==True]
        # plots
        make_lognlogs(flux_sample, bins=bins_lognlogs, color=color_cycle[s],
                      label='gamma-cat ' + sample['name'])

    # add legend and save figure
    ax1.legend(loc='lower left')
    fig1.savefig('lognlogs-{}TeV.png'.format(thresh))


