import gammalib
import numpy as np
from astropy import wcs
from astropy.io import fits
from reproject import reproject_interp
import os

# names
names = ['cygnuscocoon',
         'westerlund1']

# gas models
gasmaps = ['gas_templates/Cygnus_gas.fits',
           'gas_templates/Wd1_gas.fits']

# maximum radius in deg over which Aharonian et al determined the radial emissivity profile
rad_max = [2.5, 1.2]

# positions of the clusters from SIMBAD
glons = [80.22, 339.55]
glats = [0.80, -0.40]

# resolution of output maps in deg
out_res = 0.05

# set output in the templates directory
outdir = '../known-sources/templates/'

# spectral models
# the spectrum for the Cygnus cocoon comes from the joint ARGO+Fermi fit in Bartoli et al.  2014ApJ...790..152B
# with the cutoff at 40 TeV suggested by Milagro (Fig. 3)
spectrum_cyg = gammalib.GModelSpectralExpPlaw(3.5e-15, -2.16,
                                              gammalib.GEnergy(0.1, 'TeV'),
                                              gammalib.GEnergy(40, 'TeV'))
# the spectrum of Westerlund 1 comes from HESS collaboration 2011 A&A...537A.114A
spectrum_wd1 = gammalib.GModelSpectralPlaw(9.e-18, -2.19,
                                           gammalib.GEnergy(1, 'TeV'))
spectral_models = [spectrum_cyg, spectrum_wd1]

for s, name in enumerate(names):
    # first, create map

    # create wcs for output map
    npix = 1 + int(2 * rad_max[s] / out_res)
    out_wcs = wcs.WCS(naxis=2)
    out_wcs.wcs.crpix = [(npix - 1.) / 2 + 1., (npix - 1) / 2 + 1.]
    out_wcs.wcs.cdelt = [-out_res, out_res]
    out_wcs.wcs.crval = [glons[s], glats[s]]  # centered on cluster
    out_wcs.wcs.ctype = ["GLON-TAN", "GLAT-TAN"]

    # open gas map
    gasmap = fits.open(gasmaps[s])[0]

    # create empty output map and convert output wcs to header
    out_map = np.zeros([npix, npix])
    hdu = fits.PrimaryHDU(out_map)
    hdu.header = out_wcs.to_header()

    # reproject gas map on output grid
    array, footprint = reproject_interp(gasmap, hdu.header)
    # set output map to reporjected gas map
    hdu.data += array

    # calculate distance from center
    # since the area is small and we use a TAN projection
    # we do this in pixel coordinates
    pixnum = np.arange(npix) - ((npix - 1.) / 2)
    xv, yv = np.meshgrid(pixnum, pixnum)
    r = np.sqrt(xv ** 2 + yv ** 2)

    # central pixels will have r ~ 0
    # clip to avoid zero division error
    # clip at 0.25 deg, consistent with saturation of 1/r profiles in Aharonian et al 2019
    r[r * out_res < 0.25] = 0.25 / out_res

    # divide pixel value by distance (1/r profile as in Aharonian+ 2019)
    hdu.data /= r

    # clip tempate
    # calculate median flux at max rad
    low_flux = np.median(hdu.data[np.abs(r * out_res - rad_max[s]) < 4 * out_res])
    # clip at this flux level
    hdu.data[hdu.data < low_flux] = 0.

    # normalize map
    hdu.data /= np.sum(hdu.data)

    # create fits map file
    hdu.verify('fix')
    mapname = name + '_map.fits'
    hdu.writeto(mapname, overwrite=True)

    # create gammalib spatial model component based on map
    spatial = gammalib.GModelSpatialDiffuseMap(mapname)

    # create model and save as xml
    model = gammalib.GModelSky(spatial, spectral_models[s])
    model.name(name)
    models = gammalib.GModels()
    models.append(model)
    models.save(name + '.xml')

    # move output files to template directory
    os.system('mv {}* {}'.format(name,outdir))
