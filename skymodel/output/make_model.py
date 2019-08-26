import gammalib
import numpy as np
from gammapy.catalog import SourceCatalogGammaCat

# inputs from external sources
gammacat_file = '/Users/ltibaldo/Software/GitHub/gamma-cat/output/gammacat.fits.gz'

# create model container
models = gammalib.GModels()

# add sources from gamma-cat, keep tracks of their IDs
gammacat_ids = []
gammacat = SourceCatalogGammaCat(gammacat_file).table
for source in gammacat:
    # retain only sources with known spectral model
    # and centered within 10 degrees from the galactic pla
    if not source['spec_type'] == 'none' and (np.abs(source['glat']) < 10.):
        skip = False
        # retrieve source spatial model
        if source['morph_type'] == 'point' or source['morph_type'] == 'none':
            # source with morphology none are treated as point sources
            ra = np.double(source['ra'])
            dec = np.double(source['dec'])
            spatial = gammalib.GModelSpatialPointSource(ra, dec)
        else:
            print(
                'WARNING: source {} from gamma-cat has spatial model of type {} which is not implemented'.format(
                    source['common_name'], source['morph_type']))
            skip = True
        # retrieve source spectral model

        # put together and append to model container
        if not skip:
            gammacat_ids.append(source['source_id'])


# add CTA background
# power law spectral correction with pivot energy at 1 TeV
spectral = gammalib.GModelSpectralPlaw(1, 0, gammalib.GEnergy(1, 'TeV'))
# Irf Background
bkgmodel = gammalib.GCTAModelIrfBackground(spectral)
bkgmodel.name('Background')
bkgmodel.instruments('CTA')
# append to models
models.append(bkgmodel)

# save models
models.save('models_gps.xml')
