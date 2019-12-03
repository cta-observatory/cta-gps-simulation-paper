import sys
import gammalib
from utils import *
import numpy as np
import matplotlib.pyplot as plt

# first input is XML file name
models = gammalib.GModels(sys.argv[1])
# second and third input are minimum and maximum energy in TeV
emin = float(sys.argv[2])
emax = float(sys.argv[3])

lons, lats, radii, fluxes, names = dist_from_gammalib(models, emin=emin,emax=emax)

# binning
logs_min = int(np.floor(np.log10(np.min(fluxes))))
logs_max = int(np.ceil(np.log10(np.max(fluxes))))
nbins = 10 * (logs_max - logs_min)
bins_lognlogs = np.logspace(logs_min, logs_max, nbins)

fig1 = plt.figure('LogNLogS')
ax1 = plt.subplot()
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel("Flux > {} TeV (Crab units)".format(emin), fontsize=14)
ax1.set_ylabel('Number of sources (> Flux)', fontsize=14)
format_ax(ax1)

ax1.hist(fluxes, bins=bins_lognlogs, density=False, histtype='step', cumulative=-1)

plt.show()