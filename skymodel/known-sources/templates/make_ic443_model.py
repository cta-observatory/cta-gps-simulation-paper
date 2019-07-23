import gammalib

from normalize_template import normalize_template

"""
make XML model for IC 443
based on data from VERITAS and Fermi-LAT
Humensky, B. et al., 34th ICRC (2015); arXiv:1512.01911
"""

# input map
in_fits = 'IC443cutout.fits'
in_hdu = 1

# input spectrum
# smoothly broken power law
index1 = 2.212
index2 = 2.941
ebreak = gammalib.GEnergy(5.82e1, 'GeV')
eref = gammalib.GEnergy(10, 'GeV')
beta = 0.1  # smoothness parameter
flux = 11.05e-12  # ph cm-2 s-1
flux_ethresh = gammalib.GEnergy(200., 'GeV')

# output files
out_fits = 'map_ic443.fits'
out_xml = 'model_ic443.xml'

###############################################################

# normalize template
normalize_template(in_fits, in_hdu, out_fits)

# create gammalib model
# spatial model
spatial = gammalib.GModelSpatialDiffuseMap(out_fits)
# spectral model
spectral = gammalib.GModelSpectralSmoothBrokenPlaw(1., index1, eref, index2, ebreak, beta)
# change prefactor to get correct flux
norm = spectral.flux(flux_ethresh,gammalib.GEnergy(100,'TeV'))
spectral['Prefactor'].value(flux/norm)

# create source model
source = gammalib.GModelSky(spatial,spectral)

# create model container and append source
models = gammalib.GModels()
models.append(source)

# write to disk
models.save(out_xml)