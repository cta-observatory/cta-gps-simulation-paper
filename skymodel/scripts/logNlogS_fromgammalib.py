import sys
import gammalib
from utils import *
import numpy as np
import matplotlib.pyplot as plt

bins_lognlogs = np.logspace(-4, 1., 40)
models = gammalib.GModels(sys.argv[1])
lons, lats, radii, fluxes, names = dist_from_gammalib(models)

fig1 = plt.figure('LogNLogS')
ax1 = plt.subplot()
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel("Flux > 1 TeV (Crab units)", fontsize=14)
ax1.set_ylabel('Number of sources (> Flux)', fontsize=14)
format_ax(ax1)

ax1.hist(fluxes, bins=bins_lognlogs, density=False, histtype='step', cumulative=-1)

plt.show()