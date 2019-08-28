import gammalib
import matplotlib.pyplot as plt
import numpy as np
from gammapy.catalog import SourceCatalogGammaCat
from utils import *
import os

# inputs from external sources
gammacat_file = '/Users/ltibaldo/Software/GitHub/gamma-cat/output/gammacat.fits.gz'

# initialize diagnostic plots
fig1 = plt.figure('LogNLogS')
ax1 = plt.subplot()
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel("Flux (Crab units)", fontsize=14)
ax1.set_ylabel('Number of sources (> Flux)', fontsize=14)
format_ax(ax1)

fig2 = plt.figure('GLON')
ax2 = plt.subplot()
ax2.set_xlabel('Galactic longitude (deg)', fontsize=14)
ax2.set_ylabel('Number of sources', fontsize=14)
format_ax(ax2)

fig3 = plt.figure('GLAT')
ax3 = plt.subplot()
ax3.set_xlabel('Galactic latitude (deg)', fontsize=14)
ax3.set_ylabel('Number of sources', fontsize=14)
format_ax(ax3)

# define binning to make distributions
bins_lognlogs = np.logspace(-4, 1., 40)
bins_lon = np.linspace(-180, 180, 90)
bins_lat = np.linspace(-10, 10, 100)

# create model container
models = gammalib.GModels()

# add sources from gamma-cat, keep tracks of their IDs, longitudes, latitudes, and fluxes
gammacat_ids = []
gammacat_lons = []
gammacat_lats = []
gammacat_flux = []
gammacat = SourceCatalogGammaCat(gammacat_file).table
for source in gammacat:
    # retain only sources with known spectral model
    # and centered within 10 degrees from the galactic pla
    if not source['spec_type'] == 'none' and (np.abs(source['glat']) <= 10.):
        skip = False
        # retrieve source spatial model
        # source direction
        if np.isnan(source['pos_ra']):
            # if morphology not measured in gamma rays use value from SIMBAD
            ra = np.double(source['ra'])
            dec = np.double(source['dec'])
            lon = source['glon']
            lat = source['glat']
        else:
            # otherwise use gamma-ray measurement
            ra = np.double(source['pos_ra'])
            dec = np.double(source['pos_dec'])
            lon = source['pos_glon']
            lat = source['pos_glat']
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
                # WARNING: not sure position angle definition in gamma-cat matches gammalib
                # this looks consistent with values for Vela X, emailed Christoph on 27/08/2019 to ask (LT)
                if source['morph_pa_frame'] == 'radec':
                    # morphology defined in celestial coordinates
                    sigma2 = np.double(source['morph_sigma2'])
                    pa = np.double(source['morph_pa'])
                    spatial = gammalib.GModelSpatialEllipticalGauss(src_dir, sigma, sigma2, pa)
                elif source['morph_pa_frame'] == 'galactic':
                    # need to implement!!!
                    sigma2 = np.double(source['morph_sigma2'])
                    pa = np.double(source['morph_pa'])
                    spatial = gammalib.GModelSpatialEllipticalGauss(src_dir, sigma, sigma2, pa)
                else:
                    print(
                        'WARNING: elliptical source {} from gamma-cat has spatial model frame of type {} which is not implemented'.format(
                            source['common_name'], source['morph_pa_frame']))
                    skip = True
        elif source['morph_type'] == 'shell':
            if np.isnan(source['morph_sigma2']):
                # only average radius known, assume with of shell is 30% of radius
                rad = np.double(source['morph_sigma'])
                r_inner = rad * 0.85
                width = 0.3 * rad
            else:
                r_inner = np.double(source['morph_sigma'])
                r_outer = np.double(source['morph_sigma2'])
                width = r_outer - r_inner
            spatial = gammalib.GModelSpatialRadialShell(src_dir, r_inner, width)
        else:
            print(
                'WARNING: source {} from gamma-cat has spatial model of type {} which is not implemented'.format(
                    source['common_name'], source['morph_type']))
            skip = True
        # retrieve source spectral model
        if source['spec_type'] == 'pl' or source['spec_type'] == 'pl2':
            # treat together to handle fake PeVatron correction
            # index
            if source['spec_type'] == 'pl':
                index = source['spec_pl_index']
                norm = source['spec_pl_norm']
                eref = gammalib.GEnergy(np.double(source['spec_pl_e_ref']), 'TeV')
            else:
                index = source['spec_pl2_index']
                flux = source['spec_pl2_flux']
                emin = source['spec_pl2_e_min']
                emax = source['spec_pl2_e_max']
                if np.isnan(emax):
                    emax = 1.e6
                # for sources modeled by PL2 take pivot energy at 1 TeV
                eref = gammalib.GEnergy(1., 'TeV')
                norm = flux * (-index + 1) / (
                        np.power(emax, -index + 1) - np.power(emin, -index + 1))
            # convert norm from ph cm-2 s-1 TeV-1 (gamma-cat) to ph cm-2 s-1 MeV-1 (gammalib)
            norm *= 1.e-6
            if False:
                # implement fake PeVatron correction
                pass
            else:
                index = np.double(index)
                norm = np.double(norm)
                spectral = gammalib.GModelSpectralPlaw(norm, -index, eref)
        elif source['spec_type'] == 'ecpl':
            index = source['spec_ecpl_index']
            norm = source['spec_ecpl_norm']
            eref = gammalib.GEnergy(np.double(source['spec_ecpl_e_ref']), 'TeV')
            ecut = gammalib.GEnergy(np.double(source['spec_ecpl_e_cut']), 'TeV')
            # convert norm from ph cm-2 s-1 TeV-1 (gamma-cat) to ph cm-2 s-1 MeV-1 (gammalib)
            norm *= 1.e-6
            index = np.double(index)
            norm = np.double(norm)
            spectral = gammalib.GModelSpectralExpPlaw(norm, -index, eref, ecut)
        else:
            print(
                'WARNING: source {} from gamma-cat has spectral model of type {} which is not implemented'.format(
                    source['common_name'], source['spec_type']))
            skip = True
        # put together and append to model container
        if not skip:
            model = gammalib.GModelSky(spatial, spectral)
            model.name(source['common_name'])
            models.append(model)
            gammacat_ids.append(source['source_id'])
            gammacat_lons.append(lon)
            gammacat_lats.append(lat)
            gammacat_flux.append(source['spec_flux_1TeV_crab'])

print('Added {} gamma-cat sources'.format(len(gammacat_ids)))

# renormalize gamma-cat flux so that they are in Crab units
gammacat_flux = np.array(gammacat_flux)
gammacat_flux *= 1.e-2

# change lon range from 0...360 to -180...180
gammacat_lons = np.array(gammacat_lons)
gammacat_lons[gammacat_lons > 180] = gammacat_lons[gammacat_lons > 180] - 360.

# make first distributions
ax1.hist(gammacat_flux, bins=bins_lognlogs, density=False, histtype='step', cumulative=-1,
         label='gamma-cat', alpha=0.5, linewidth=2)
ax2.hist(gammacat_lons, bins=bins_lon, density=False, histtype='step',
         label='gamma-cat', alpha=0.5, linewidth=2)
ax3.hist(gammacat_lats, bins=bins_lat, density=False, histtype='step',
         label='gamma-cat', alpha=0.5, linewidth=2)

# make distributions from gammalib model container
lons, lats, fluxes = dist_from_gammalib(models)
# change lon range from 0...360 to -180...180
lons = np.array(lons)
lons[lons > 180] = lons[lons > 180] - 360.
ax1.hist(fluxes, bins=bins_lognlogs, density=False, histtype='step', cumulative=-1,
         label='converted gamma-cat', alpha=0.5, linewidth=2, linestyle=':')
ax2.hist(lons, bins=bins_lon, density=False, histtype='step',
         label='converted gamma-cat', alpha=0.5, linewidth=2, linestyle=':')
ax3.hist(lats, bins=bins_lat, density=False, histtype='step',
         label='converted gamma-cat', alpha=0.5, linewidth=2, linestyle=':')

# add templates

# read template list
template_list = open('../known-sources/templates/templates.dat').readlines()

replaced = 0
added = 0
for template in template_list:
    if template[0] == '#':
        # skip header and commented lines
        pass
    else:
        name, id = template.split(',')
        # remove gamma-cat model if present
        if 'NONE' in id:
            # source not included in gammacat, pass
            added += 1
        else:
            # remove gamma-cat model
            id = int(id)
            # if source included in gammacat with spectral information, remove gammacat model
            source = gammacat[gammacat['source_id']==id][0]
            if source['spec_type'] == 'none' or np.abs(source['glat']) > 10:
                added +=1
            else:
                models.remove(source['common_name'])
                replaced += 1
        # add templates
        models_template = gammalib.GModels('../known-sources/templates/{}.xml'.format(name))
        for model in models_template:
            models.append(model)
        # copy FITS maps to output repository
        os.system('cp ../known-sources/templates/{}*_map.fits ../output/'.format(name))

print('Replaced {} gamma-cat sources with templates. Added {} sources as templates'.format(replaced,added))

# re-make distributions from gammalib model container
lons, lats, fluxes = dist_from_gammalib(models)
# change lon range from 0...360 to -180...180
lons = np.array(lons)
lons[lons > 180] = lons[lons > 180] - 360.
ax1.hist(fluxes, bins=bins_lognlogs, density=False, histtype='step', cumulative=-1,
         label='gamma-cat + templates', alpha=0.5, linewidth=2, linestyle=':')
ax2.hist(lons, bins=bins_lon, density=False, histtype='step',
         label='gamma-cat + templates', alpha=0.5, linewidth=2, linestyle=':')
ax3.hist(lats, bins=bins_lat, density=False, histtype='step',
         label='gamma-cat + templates', alpha=0.5, linewidth=2, linestyle=':')

# add binaries

# add pulsars

# add SFR

# add HAWC

# add 3FHL

# add CHECKS for duplicated sources !!!!!!!!!!!!!!!!!!

# add synthetic PWNe and SNRs

# add IEM

# add CTA background
# power law spectral correction with pivot energy at 1 TeV
spectral = gammalib.GModelSpectralPlaw(1, 0, gammalib.GEnergy(1, 'TeV'))
# Irf Background
bkgmodel = gammalib.GCTAModelIrfBackground(spectral)
bkgmodel.name('Background')
bkgmodel.instruments('CTA')
# append to models
models.append(bkgmodel)
print('Added background model')

# save models
models.save('../output/models_gps.xml')

# save diagnostic plots
ax1.legend()
fig1.savefig('../output/logNlogS.png', dpi=300)
ax2.legend()
ax2.set_xlim(180, -180)
fig2.savefig('../output/glon.png', dpi=300)
ax3.legend()
fig3.savefig('../output/glat.png', dpi=300)
