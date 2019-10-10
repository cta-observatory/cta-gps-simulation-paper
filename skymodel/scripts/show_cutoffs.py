import matplotlib.pyplot as plt
import numpy as np

results = [{'filename': '../output/ecut_pwn.npy', 'name': 'PWN'},
           {'filename': '../output/ecut_snr.npy', 'name': 'SNR'},
           {'filename': '../output/ecut_agn.npy', 'name': 'AGN'},
           {'filename': '../output/ecut_unid.npy', 'name': 'UNID'}]

arrays = []
for result in results:
    arrays.append(np.load(result['filename']))

bins = np.logspace(-1, 3, 50)

fig = plt.figure()
ax = plt.subplot()
ax.set_xscale('log')
ax.set_xlabel("Cutoff energy (TeV)")
ax.set_ylabel("Number of sources")
ax.set_title("Artificial cutoffs for detected sources")

histos = []
ntot = np.zeros(len(bins) - 1)
for s, array in enumerate(arrays):
    n, bins_out, patches = plt.hist(array, bins=bins, density=False, label=results[s]['name'],
                                    bottom=ntot)
    ntot += n

ax.legend()

plt.savefig('../output/cutoffs.png')
