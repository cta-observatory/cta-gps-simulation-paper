import gammalib
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil
from astropy.table import Table
from add_fhl import *
from add_hawc import *
from utils import *
from cutoffs import get_cutoff, get_atnf_version

# inputs from external sources
gammacat_file = '../known-sources/external-input/gammacat.fits.gz'

# maximum latitude to include in the model
bmax = 10.

# go to output directory as working directory
# this simplifies file path handling
os.chdir('../output')

# create report file
outfile = open('report.txt', 'w')

# record version of ATNF catalog used
msg = 'Using ATNF catalog version {}\n'.format(get_atnf_version())
print(msg)
outfile.write(msg)

# initialize diagnostic plots
fig1 = plt.figure('LogNLogS')
ax1 = plt.subplot()
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel("Flux > 1 TeV (Crab units)", fontsize=14)
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
bins_lat = np.linspace(-bmax, bmax, 10 * bmax)

# create model container
models = gammalib.GModels()

print('\n')

# add sources from gamma-cat, keep tracks of their IDs, longitudes, latitudes, and fluxes
gammacat_ids = []
gammacat_lons = []
gammacat_lats = []
gammacat_flux = []
# keep track also of artificial cutoffs
ecut_pwn = []
ecut_snr = []
ecut_unid = []
gammacat = Table.read(gammacat_file)
for source in gammacat:
    # retain only sources with known spectral model
    # and centered within 10 degrees from the galactic pla
    if not source['spec_type'] == 'none' and (np.abs(source['glat']) <= bmax):
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
                ##############################################
                # temporary fix for missing sigma of Westerlund 1
                # to be removed once Westerlund 1 dedicated model is added!!!
                if source['common_name'] == 'Westerlund 1':
                    sigma = 1.1
                # symmetric Gaussian
                spatial = gammalib.GModelSpatialRadialGauss(src_dir, sigma)
            else:
                # elliptical Gaussian
                # WARNING: not sure position angle definition in gamma-cat matches gammalib
                # this looks consistent with values for Vela X, emailed Christoph on 27/08/2019 to ask (LT)
                sigma2 = np.double(source['morph_sigma2'])
                pa = np.double(source['morph_pa'])
                if source['morph_pa_frame'] == 'radec':
                    # morphology defined in celestial coordinates
                    spatial = gammalib.GModelSpatialEllipticalGauss(src_dir, sigma, sigma2, pa)
                elif source['morph_pa_frame'] == 'galactic':
                    # rotate position angle because gammalib only accepts elliptical models in celestial coordinates
                    gpole = gammalib.GSkyDir()
                    # North Galactic pole
                    gpole.radec_deg(192.858333, 27.128333)
                    pa = pa + src_dir.posang_deg(gpole)
                    if pa < 0:
                        pa = 360. + pa
                    spatial = gammalib.GModelSpatialEllipticalGauss(src_dir, sigma, sigma2, pa)
                else:
                    msg = 'WARNING: elliptical source {} from gamma-cat has spatial model frame of type {} which is not implemented\n'.format(
                        source['common_name'], source['morph_pa_frame'])
                    print(msg)
                    outfile.write(msg)
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
            msg = 'WARNING: source {} from gamma-cat has spatial model of type {} which is not implemented\n'.format(
                source['common_name'], source['morph_type'])
            print(msg)
            outfile.write(msg)
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
            index = np.double(index)
            norm = np.double(norm)
            # Set artifical cutoffs for sources measured with hard PL spectra without measured cutoff
            # Binaries and pulsars are skipped because they are addressed later
            # We also leave Westerlund 1 and HESS J1641-463 as PeVatron candidates
            if index < 2.4 and not source['classes'] == 'bin'\
                    and not source['classes'] == 'psr'\
                    and not source['common_name'] == 'Westerlund 1'\
                    and not source['common_name'] == 'HESS J1641-463':
                #print('set artificial cutoff for source {} of type {}'.format(source['common_name'], source['classes']))
                if 'pwn' in source['classes']:
                    # dummy model to obtain search radius
                    mod = gammalib.GModelSky(spatial,gammalib.GModelSpectralPlaw())
                    # search radius
                    rad = get_model_radius(mod) + 0.2
                    # set cutoff
                    ecut = get_cutoff(ra,dec,'PSR',rad_search=rad)
                    ecut_pwn.append(ecut)
                # elif 'snr' in source['classes']:
                #     # try to get Green name
                #     onames = source['other_names'].split(',')
                #     gname = [name for name in test if name[:5] == 'SNR G' or name[0] == 'G']
                #     if len(gname)>0:
                #         gname = gname[0]
                #     else:
                #         gname = None
                #     # verify if interaction with molecular clouds is listed
                #     if 'mc' in source['classes']:
                #         interacting = True
                #     else:
                #         interacting = False
                else:
                    ecut = get_cutoff(ra, dec, 'X')
                    ecut_unid.append(ecut)
                ecut = gammalib.GEnergy(np.double(ecut), 'TeV')
                spectral = gammalib.GModelSpectralExpPlaw(norm, -index, eref, ecut)
            else:
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
            msg = 'WARNING: source {} from gamma-cat has spectral model of type {} which is not implemented\n'.format(
                source['common_name'], source['spec_type'])
            print(msg)
            outfile.write(msg)
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

msg = 'Added {} gamma-cat sources\n'.format(len(gammacat_ids))
print(msg)
outfile.write(msg)


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
            # remove gamma-cat model if present
            id = int(id)
            if id in gammacat_ids:
                source = gammacat[gammacat['source_id'] == id][0]
                models.remove(source['common_name'])
                replaced += 1
            else:
                added += 1
        # add templates
        models_template = gammalib.GModels('../known-sources/templates/{}.xml'.format(name))
        for model in models_template:
            # find model map name and path
            filename = model.spatial().filename().file()
            filepath = model.spatial().filename().path()
            # copy file to output directory
            shutil.copy(filepath + filename, './')
            # replace file with the one in output directory
            model.spatial(gammalib.GModelSpatialDiffuseMap(filename))
            models.append(model)

msg = 'Replaced {} gamma-cat sources with templates. Added {} sources as templates\n'.format(
    replaced, added)
print(msg)
outfile.write(msg)

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

# read binary list
binary_list = open('../bin/binaries.dat').readlines()

# read models
binary_models = gammalib.GModels('../bin/models_binaries.xml')

replaced = 0
added = 0
for binary in binary_list:
    if binary[0] == '#':
        # skip header and commented lines
        pass
    else:
        id, bmodels = binary.split(',')
        if 'None' in id:
            # source not included in gammacat, pass
            added += 1
        else:
            # remove gamma-cat model if present
            id = int(id)
            if id in gammacat_ids:
                source = gammacat[gammacat['source_id'] == id][0]
                models.remove(source['common_name'])
                replaced += 1
            else:
                added += 1
        # add binary model
        # find names
        try:
            # check if multiple models are present
            model_names = bmodels.split(':')
        except:
            # otherwise set single name as list
            model_names = [models]
        for model_name in model_names:
            model = binary_models[model_name.strip()]
            # find phasecurve name and path
            filename = model.temporal().filename().file()
            filepath = model.temporal().filename().path()
            # copy file to output directory
            shutil.copy(filepath + filename, './')
            # replace file with the one in output directory
            model.temporal().filename(filename)
            models.append(model)

msg = 'Replaced {} gamma-cat sources with binaries. Added {} sources as binaries\n'.format(
    replaced, added)
print(msg)
outfile.write(msg)

# add pulsars

# re-make distributions from gammalib model container
lons, lats, fluxes = dist_from_gammalib(models)
# change lon range from 0...360 to -180...180
lons = np.array(lons)
lons[lons > 180] = lons[lons > 180] - 360.
ax1.hist(fluxes, bins=bins_lognlogs, density=False, histtype='step', cumulative=-1,
         label='gamma-cat + templates + bin', alpha=0.5, linewidth=2, linestyle=':')
ax2.hist(lons, bins=bins_lon, density=False, histtype='step',
         label='gamma-cat + templates + bin', alpha=0.5, linewidth=2, linestyle=':')
ax3.hist(lats, bins=bins_lat, density=False, histtype='step',
         label='gamma-cat + templates + bin', alpha=0.5, linewidth=2, linestyle=':')

# add SFR

# add 3FHL

models, newpt, newext = append_fhl(models,bmax,dist_sigma=3.)

msg = 'Added {} FHL sources, of which {} as pointlike and {} as extended.\n'.format(
    newpt+newext, newpt,newext)
print(msg)
outfile.write(msg)

# add HAWC

models, newpt, newext = append_hawc(models,bmax,dist_sigma=3.)

msg = 'Added {} HAWC sources, of which {} as pointlike and {} as extended.\n'.format(
    newpt+newext, newpt,newext)
print(msg)
outfile.write(msg)

# save cutoff values
np.save('ecut_pwn.npy',ecut_pwn)
np.save('ecut_snr.npy',ecut_snr)
np.save('ecut_unid.npy',ecut_unid)

# re-make distributions from gammalib model container
lons, lats, fluxes = dist_from_gammalib(models)
# change lon range from 0...360 to -180...180
lons = np.array(lons)
lons[lons > 180] = lons[lons > 180] - 360.
ax1.hist(fluxes, bins=bins_lognlogs, density=False, histtype='step', cumulative=-1,
         label='gamma-cat + templates + bin + FHL + HAWC', alpha=0.5, linewidth=2, linestyle=':')
ax2.hist(lons, bins=bins_lon, density=False, histtype='step',
         label='gamma-cat + templates + bin + FHL + HAWC', alpha=0.5, linewidth=2, linestyle=':')
ax3.hist(lats, bins=bins_lat, density=False, histtype='step',
         label='gamma-cat + templates + bin + FHL + HAWC', alpha=0.5, linewidth=2, linestyle=':')

# CHECKS for duplicated sources
for model1 in models:
    dir1 = get_model_dir(model1)
    for model2 in models:
        if model2.name() == model1.name():
            pass
        else:
            dir2 = get_model_dir(model2)
            dist_deg = dir1.dist_deg(dir2)
            if dist_deg < 0.2:
                msg = 'WARNING: sources {} and {} are at a distance of {} deg lower than 0.2 deg\n'.format(
                    model1.name(), model2.name(), dist_deg)
                print(msg)
                outfile.write(msg)

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
outfile.write('Added background model\n')

# save models
models.save('models_gps.xml')

# close report file
outfile.close()

# save diagnostic plots
ax1.legend(fontsize=5)
fig1.savefig('logNlogS.png', dpi=300)
ax2.legend(fontsize=5)
ax2.set_xlim(180, -180)
fig2.savefig('glon.png', dpi=300)
ax3.legend(fontsize=5)
fig3.savefig('glat.png', dpi=300)
