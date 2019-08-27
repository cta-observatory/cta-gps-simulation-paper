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
        # source direction
        ra = np.double(source['ra'])
        dec = np.double(source['dec'])
        src_dir = gammalib.GSkyDir()
        src_dir.radec_deg(ra, dec)
        if source['morph_type'] == 'point' or source['morph_type'] == 'none':
            # source with morphology none are treated as point sources
            spatial = gammalib.GModelSpatialPointSource(src_dir)
        elif source['morph_type'] == 'gauss':
            # Gaussian morphology
            sigma = np.double(source['morph_sigma'])
            if np.isnan(source['morph_sigma2']):
                # symmetric Gaussian
                spatial = gammalib.GModelSpatialRadialGauss(src_dir, sigma)
            else:
                # elliptical Gaussian
                # WARNING: not sure position angle definition in gamma-cat matches gammalib, emailed Christoph on 27/08/2019 to ask (LT)
                if source['morph_pa_frame'] == 'radec':
                    # morphology defined in celestial coordinates
                    sigma2 = np.double(source['morph_sigma2'])
                    pa = np.double(source['morph_pa'])
                    spatial = gammalib.GModelSpatialEllipticalGauss(src_dir, sigma, sigma2, pa)
                elif source['morph_pa_frame'] == 'galactic':
                    # need to implement!!!
                    pass
                else:
                    print('WARNING: source {} from gamma-cat has spatial model frame of type {} which is not implemented'.format(
                            source['common_name'], source['morph_pa_frame']))
                    skip = True
        elif source['morph_type'] == 'shell':
            pass
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
models.save('../output/models_gps.xml')
