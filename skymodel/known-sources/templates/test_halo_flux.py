from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
import numpy as np
import matplotlib.pyplot as plt

# cross check, produces spectral plots for halos to compare with literature


filenames = ['gemingahalo_map.fits','psrb0656+14halo_map.fits']
psr = [SkyCoord([98.475637], [17.770253], frame="icrs", unit="deg"),
       SkyCoord([104.950762], [14.239317], frame="icrs", unit="deg")]
rad = [11,11]


def flux_plot(filename,c,rad):
    hdulist = fits.open(filename)

    cube = hdulist[0].data
    energies = hdulist[1].data['Energy']

    # for each energy calculate integrated flux in extraction from paper
    # identify pixel corresponding to center of extraction region
    w = wcs.WCS(hdulist[0].header)
    xp, yp = skycoord_to_pixel(c,w)
    # create mask to set flux to 0 outside extraction region
    nx = np.abs(hdulist[0].header['NAXIS1'])
    ny = np.abs(hdulist[0].header['NAXIS2'])
    mask = np.ones([ny,nx])
    binsize = np.abs(hdulist[0].header['CDELT1'])
    for y in range(ny):
        for x in range (nx):
            if (x - xp)** 2 + (y - yp)**2 < (rad/binsize)**2:
                pass
            else:
                mask[y,x] = 0.
    # multiply array by mask
    cube *= mask[np.newaxis,:,:]
    # multiply by pixels solid angle
    cube *= np.deg2rad(binsize)**2
    # sum over extraction region
    fluxes = np.sum(cube,axis=(1,2))

    # multiply flux by E2
    fluxes *= energies**2
    # energies from MeV to GeV
    energies *= 1.e-3
    # energy flux from MeV to GeV
    fluxes *= 1.e-3

    # plot
    fig = plt.figure()
    ax = plt.subplot()
    ax.loglog(energies,fluxes)
    ax.set_xlabel(r'Energy (GeV)')
    ax.set_ylabel(r'E$^{2}$ dN/dE (GeV cm$^{-2}$ s$^{-1}$)')
    ax.grid()

for s, filename in enumerate(filenames):
    flux_plot(filename,psr[s],rad[s])

plt.show()

