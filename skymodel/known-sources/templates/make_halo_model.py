import pandas as pd
import numpy as np
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u
import gammalib
import pdb

"""
make templates for Geminga and PSR B0656+14 halos
material provided by Ruben Lopez Coto, note that Ruben sent pasted below

Everything was calculated using EDGE (https://arxiv.org/abs/1709.07653).
Results are saved in  hdf5 format in two dataframes, one for energy [MeV] (20 GeV - 200 TeV)
the other with flux per angular bin [MeV^-1 cm^-2 s^-1 sr^-1].
"""

# positions of the two PSR at the halo centers as astropy SkyCoord object
pos_gem = SkyCoord([98.475637], [17.770253], frame="icrs", unit="deg")
pos_b06 = SkyCoord([104.950762], [14.239317], frame="icrs", unit="deg")

# profiles
profile_gem = 'energy_flux_Geminga.h5'
profile_b06 = 'energy_flux_PSR_B0656.h5'

# energy vectors
energy_gem = 'energy_Geminga.h5'
energy_b06 = 'energy_PSR_B0656.h5'


# generate map from radial profiles
def radprof2map(pos, energy, profile, nametag):
    # read energies and angular profiles
    energies = np.array(pd.read_hdf(energy))[:,0]
    fluxes = np.array(pd.read_hdf(profile))
    angles = np.array(pd.read_hdf(profile).keys())

    # determine number of pixels to cover the entire profile
    # reduce resolution by factor of 2
    # + 1 makes the final map centered on the pulsar
    npix = len(angles) + 1

    # create wcs for output map
    # output binning will cover the whole model with npix pixels
    out_res = 2 * angles[-1] / (npix - 1)
    out_wcs = wcs.WCS(naxis=3)
    out_wcs.wcs.crpix = [(npix - 1.) / 2 + 1., (npix - 1) / 2 + 1., 1]
    out_wcs.wcs.cdelt = [-out_res, out_res,
                         np.ediff1d(np.log10(energies))[0]]
    out_wcs.wcs.crval = [pos.ra.deg[0], pos.dec.deg[0], # centered on pulsar
                         np.log10(energies[0])]
    out_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN", "Log10(Energy/1 MeV)"]

    # create output map
    out_map = np.zeros([len(energies),npix,npix])
    # create pixel arrays
    ii = np.linspace(1,npix,npix)
    xv, yv = np.meshgrid(ii, ii)
    xx = xv.flatten()
    yy = yv.flatten()
    # fake array for energy
    zz = np.zeros(np.shape(xx))
    # world position array
    world = out_wcs.wcs_pix2world(np.array([xx,yy,zz]).T,1)
    # sky coordinate array
    sc = SkyCoord(world[:, 0] * u.deg, world[:, 1] * u.deg, frame='icrs')
    # angular distance array in deg
    sep = sc.separation(pos).base.base
    for k in range(len(energies)):
        # get flux interpolating over angular bins
        flux = np.interp(sep, angles, fluxes[k])
        # set flux to 0 beyond range covered by model
        flux[sep > angles[-1]] = 0.
        # reshape array and fill map
        out_map[k] = flux.reshape(npix,npix)

    # create the main hdu
    hdu = fits.PrimaryHDU(out_map)
    hdu.header = out_wcs.to_header()
    hdu.header.set('BUNIT', 'photon/cm2/s/MeV/sr', 'Photon flux', after='CRVAL3')
    hdu.verify('fix')

    # create the energy table
    ecol = fits.Column(name='Energy', format='D', unit='MeV', array=energies)
    tbhdu = fits.BinTableHDU.from_columns([ecol], name='ENERGIES')
    tbhdu.verify('fix')

    # write file
    hdulist = fits.HDUList([hdu, tbhdu])
    hdulist.writeto(nametag + '_map.fits', overwrite=True)

radprof2map(pos_gem,energy_gem,profile_gem,'gemingahalo')
radprof2map(pos_b06,energy_b06,profile_b06,'psrb0656+14halo')

# create xml models
def make_xml(nametag):
    spatial = gammalib.GModelSpatialDiffuseCube(gammalib.GFilename(nametag + '_map.fits'))
    spectral = gammalib.GModelSpectralConst(1.)
    model = gammalib.GModelSky(spatial, spectral)
    model.name(nametag)
    # fill to model container and write to disk
    models = gammalib.GModels()
    models.append(model)
    models.save(nametag + '.xml')

make_xml('gemingahalo')
make_xml('psrb0656+14halo')




