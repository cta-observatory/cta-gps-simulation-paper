import numpy as np
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.io import fits
import gammalib

"""
make templates for Geminga and PSR B0656+14 halos
material provided by Ruben Lopez Coto, note that Rube sent pasted below

The file containing the spectrum has:
Energy [TeV]		E^2 dN/dE [erg s^-1 cm^-2]
Everything was calculated using EDGE (https://arxiv.org/abs/1709.07653), and I extended the spectrum down to 1 GeV and up to ~500 TeV, but HAWC’s measurement was between 8 and 40 TeV (5 and 50 if you do not consider the upper/lower limits of the errors). This spectrum corresponds to the integral over the full area up to 2 * r_diff, where r_diff is the diffusion radius measured (5.5 deg for Geminga and 4.8 deg for PST B0656+14).
The file containing the surface brightness has:
Angular distance from the pulsar [deg]			Surface brightness [erg / (s cm^2 deg^2)]
The surface brightness, to comply with what we published in the HAWC paper, is calculated as the integral of \int{E dN/dE} between 5 and 50 TeV in that particular angular bin.
"""

# number of pixels in output map along each axis
# odd number to be have pulsar in the central pixel
npix = 201

# positions of the two PSR at the halo centers as astropy SkyCoord object
pos_gem = SkyCoord([98.475637], [17.770253], frame="icrs", unit="deg")
pos_b06 = SkyCoord([104.950762], [14.239317], frame="icrs", unit="deg")

# profiles
profile_gem = 'Gamma_Profile_Energy_flux_Geminga.txt'
profile_b06 = 'Gamma_Profile_Energy_flux_PSRB0656.txt'

# spectra
spectrum_gem = 'Gamma_Spectra_Geminga.txt'
spectrum_b06 = 'Gamma_Spectra_PSRB0656.txt'


# generate map from radial profiles
def radprof2map(pos, profile, nametag):
    # read profile
    prof = np.genfromtxt(profile)
    dist = prof[:, 0]  # distances
    val = prof[:, 1]  # surgace brightness values

    # create wcs for output map
    # output binning will cover the whole model with npix pixels
    out_res = 2 * dist[-1] / (npix - 1)
    out_wcs = wcs.WCS(naxis=2)
    out_wcs.wcs.crpix = [(npix - 1.) / 2 + 1., (npix - 1) / 2 + 1.]
    out_wcs.wcs.cdelt = [-out_res, out_res]
    out_wcs.wcs.crval = [pos.ra.deg[0], pos.dec.deg[0]]  # centered on ~GC
    out_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    # create output map and fill
    out_map = np.zeros([npix,npix])
    for j in range(npix):
        for i in range(npix):
            world = out_wcs.wcs_pix2world([[i,j]], 1)[0]
            pos_pix = SkyCoord([world[0]], [world[1]], frame="icrs", unit="deg")
            sep = pos_pix.separation(pos).deg[0]
            out_map[j,i] = np.interp(sep,dist,val)

    # create the main hdu
    hdu = fits.PrimaryHDU(out_map)
    hdu.header = out_wcs.to_header()
    hdu.verify('fix')

    # create fits map file
    hdu.writeto(nametag + '_map.fits', overwrite=True)

radprof2map(pos_gem,profile_gem,'gemingahalo')
radprof2map(pos_b06,profile_b06,'psrb0656+14halo')

# reformat spectrum
def reformat_spectrum(spectrum,nametag):
    spec = np.genfromtxt(spectrum)
    en = spec[:, 0]
    flux = spec[:, 1]
    # convert energy from TeV to MeV
    en *= 1.e6
    # convert SED to flux per MeV
    # erg to MeV
    flux *= 6.24151e5
    # SED to flux
    flux /= en**2
    out_spec = np.zeros([len(en),2])
    out_spec[:,0] = en
    out_spec[:,1] = flux
    np.savetxt(nametag + '_spectrum.txt', out_spec)

reformat_spectrum(spectrum_gem,'gemingahalo')
reformat_spectrum(spectrum_b06,'psrb0656+14halo')

# create xml models
def make_xml(nametag):
    spatial = gammalib.GModelSpatialDiffuseMap(nametag + '_map.fits')
    spectral = gammalib.GModelSpectralFunc(gammalib.GFilename(nametag + '_spectrum.txt'),1.)
    model = gammalib.GModelSky(spatial,spectral)
    model.name(nametag)
    models = gammalib.GModels()
    models.append(model)
    models.save(nametag + '.xml')

make_xml('gemingahalo')
make_xml('psrb0656+14halo')




