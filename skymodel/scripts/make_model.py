import gammalib
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil
import datetime
from astropy.table import Table
from add_fhl import *
from add_hawc import *
from utils import *
from cutoffs import get_cutoff, get_atnf_version
from binpop import get_binpop_models

# inputs from external sources
gammacat_file = '../known-sources/external-input/gammacat.fits.gz'

# maximum latitude to include in the model
bmax = 10.

# minimum flux for synthetic source (mCrab)
fmin = 0.001

# minimum radius (deg) below which sources are treated as pointlike
# for association and plotting
radmin = 0.05

# go to output directory as working directory
# this simplifies file path handling
os.chdir('../output')

# create report file
outfile = open('report.txt', 'w')

# record time and machine on which the run is done
timenow = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
user = os.environ['HOME'].split('/')[-1]
machine = os.uname()[1]
msg = 'Model built by {} at {} on {}\n'.format(user,machine,timenow)
print(msg)
outfile.write(msg)

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

fig5 = plt.figure('LogNLogS-solidang')
ax5 = plt.subplot()
ax5.set_xscale('log')
ax5.set_yscale('log')
ax5.set_xlabel("Flux per solid angle > 1 TeV (Crab units)", fontsize=14)
ax5.set_ylabel('Number of sources (> Flux)', fontsize=14)
format_ax(ax5)

fig2 = plt.figure('GLON')
ax2 = plt.subplot()
ax2.set_xlabel('Galactic longitude (deg)', fontsize=14)
ax2.set_ylabel('Number of sources', fontsize=14)
ax2.set_ylim(0,20)
format_ax(ax2)

fig3 = plt.figure('GLAT')
ax3 = plt.subplot()
ax3.set_xlabel('Galactic latitude (deg)', fontsize=14)
ax3.set_ylabel('Number of sources', fontsize=14)
format_ax(ax3)

fig0 = plt.figure('Radius')
ax0 = plt.subplot()
ax0.set_yscale('log')
ax0.set_xlabel('Radius (deg)', fontsize=14)
ax0.set_ylabel('Number of sources', fontsize=14)
format_ax(ax0)

# define binning to make distributions
bins_lognlogs = np.logspace(np.log10(1.e-3 * fmin), 1., 50)
bins_lon = np.linspace(-180, 180, 90)
bins_lat = np.linspace(-bmax, bmax, 10 * bmax)
bins_rad = np.linspace(0, 2, 30)

# load synthetic populations, so that high flux members can be dropped as we add real sources

# binaries
bin_models, bin_dict = get_binpop_models('../binpop',fmin,'./',bmax)
msg = 'Loaded {} synthetic binaries\n'.format(bin_models.size())
print(msg)
outfile.write(msg)
# keep track of properties of binaries deleted
bin_distx = []
bin_disty = []
bin_radr = []
bin_frlog = []

# pwn
pwn_models, pwn_dict = get_syn_model('../pwn/xml/pwn.xml',
                                       1.e-3*fmin,bmax,0.1,1000.)
msg = 'Loaded {} synthetic PWNe\n'.format(pwn_models.size())
print(msg)
outfile.write(msg)
# keep track of properties of binaries deleted
pwn_distx = []
pwn_disty = []
pwn_radr = []
pwn_frlog = []

# snr
snr_models, snr_dict = get_syn_model('../snr/OUTPUT_FILES_1/ctadc_skymodel_gps_sources_pevatron_0.xml',
                                     1.e-3*fmin,bmax,emin=0.1,emax=1000.)
msg = 'Loaded {} synthetic young SNRs\n'.format(snr_models.size())
print(msg)
outfile.write(msg)
# keep track of properties of snrs deleted
snr_distx = []
snr_disty = []
snr_radr = []
snr_frlog = []

# isnr
isnr_models, isnr_dict = get_syn_model('../int-snr/out/isnr.xml',
                                       1.e-3*fmin,bmax,0.1,1000.)
msg = 'Loaded {} synthetic interacting SNRs\n'.format(isnr_models.size())
print(msg)
outfile.write(msg)
# keep track of properties of binaries deleted
isnr_distx = []
isnr_disty = []
isnr_radr = []
isnr_frlog = []

# set composite SNR/PWN systems
comp_dict, pwn_dict, snr_dict = set_composites(pwn_dict,snr_dict)
msg = '{} pairs PWN/SNR constitute composite systems\n'.format(len(comp_dict['name']))
print(msg)
outfile.write(msg)
# keep track of properties of composites deleted
comp_distx = []
comp_disty = []
comp_radr = []
comp_frlog = []

# create final model container
models = gammalib.GModels()
# create container for synthetic sources only
models_syn = gammalib.GModels()

print('\n')

# add sources from gamma-cat, keep tracks of their IDs, longitudes, latitudes, and fluxes
gammacat_ids = []
gammacat_lons = []
gammacat_lats = []
gammacat_rads = []
gammacat_flux = []
# keep track also of artificial cutoffs
ecut_pwn = []
ecut_snr = []
ecut_unid = []
ecut_agn = []
n_ecut_pwn = 0
n_ecut_snr = 0
n_ecut_unid = 0
# keep track of synthetic sources that are be deleted
n_bin_del = 0
n_snr_del = 0
n_isnr_del = 0
n_comp_del = 0
n_pwn_del = 0

# read file
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
            radius = 0.
        elif source['morph_type'] == 'gauss':
            # Gaussian morphology
            sigma = np.double(source['morph_sigma'])
            radius = 2 * sigma
            if np.isnan(source['morph_sigma2']):
                ##############################################
                # temporary fix for missing sigma of Westerlund 1
                # to be removed once Westerlund 1 dedicated model is added!!!
                if source['common_name'] == 'Westerlund 1':
                    sigma = 1.1
                # symmetric Gaussian
                # rad = 2 * sigma
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
            radius = r_inner + width
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
                # if source may be PWN treat as such
                if 'pwn' in source['classes']:
                    # dummy model to obtain search radius
                    mod = gammalib.GModelSky(spatial,gammalib.GModelSpectralPlaw())
                    # search radius
                    rad = get_model_radius(mod) + 0.2
                    # set cutoff
                    ecut = get_cutoff(ra,dec,'PSR',rad_search=rad)
                    ecut_pwn.append(ecut)
                    n_ecut_pwn+=1
                # otherwise consider SNR
                elif 'snr' in source['classes']:
                    # try to get Green name
                    onames = source['other_names'].split(',')
                    gname = [name for name in onames if name[:5] == 'SNR G' or name[0] == 'G']
                    if len(gname)>0:
                        gname = gname[0]
                    else:
                        gname = None
                    hess = True
                    # compute cutoff
                    # we use default value for particle spectral index because measurements are highly uncertain
                    ecut = get_cutoff(ra, dec, 'SNR', name = gname, hess=hess)
                    ecut_snr.append(ecut)
                    n_ecut_snr+=1
                # otherwise if unidentified
                else:
                    if 'unid' in source['classes']:
                        pass
                    else:
                        # set warning if we have hard source of unexpected type
                        msg = 'Gamma-cat source {} of type {} has an unxepctedly hard spectrum ' \
                              'with index {}. We are setting a random artificial cutoff\n'.format(source['common_name'],source['classes'],index)
                        print(msg)
                        outfile.write(msg)
                    ecut = get_cutoff(ra, dec, 'UNID')
                    ecut_unid.append(ecut)
                    n_ecut_unid +=1
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
            gammacat_rads.append(radius)
            gammacat_flux.append(source['spec_flux_1TeV_crab'])
            # find which synthetic source need to be deleted to account for the source added
            if source['classes'] == 'bin':
                rname, bin_dict, distx, disty, radr, frlog = find_source_to_delete(bin_dict,src_dir.l_deg(),src_dir.b_deg(), get_model_radius(model), 1.e-2*source['spec_flux_1TeV_crab'],radmin=radmin)
                bin_models.remove(rname)
                n_bin_del +=1
                bin_distx.append(distx)
                bin_disty.append(disty)
                bin_radr.append(radr)
                bin_frlog.append(frlog)
            elif source['classes'] == 'pwn,snr':
                rname, comp_dict, distx, disty, radr, frlog = find_source_to_delete(comp_dict,src_dir.l_deg(),src_dir.b_deg(),get_model_radius(model),1.e-2 *source['spec_flux_1TeV_crab'],radmin=radmin)
                pwn_models.remove(rname[0])
                snr_models.remove(rname[1])
                n_comp_del += 1
                comp_distx.append(distx)
                comp_disty.append(disty)
                comp_radr.append(radr)
                comp_frlog.append(frlog)
            elif source['classes'] == 'snr,mc':
                rname, isnr_dict, distx, disty, radr, frlog = find_source_to_delete(isnr_dict,src_dir.l_deg(),src_dir.b_deg(), get_model_radius(model),1.e-2 * source['spec_flux_1TeV_crab'],radmin=radmin)
                isnr_models.remove(rname)
                n_isnr_del += 1
                isnr_distx.append(distx)
                isnr_disty.append(disty)
                isnr_radr.append(radr)
                isnr_frlog.append(frlog)
            elif source['classes'] == 'snr':
                rname, snr_dict, distx, disty, radr, frlog = find_source_to_delete(snr_dict,src_dir.l_deg(),src_dir.b_deg(), get_model_radius(model),1.e-2 * source['spec_flux_1TeV_crab'],radmin=radmin)
                snr_models.remove(rname)
                n_snr_del += 1
                snr_distx.append(distx)
                snr_disty.append(disty)
                snr_radr.append(radr)
                snr_frlog.append(frlog)
            elif 'pwn' in source['classes'] or source['classes'] == 'unid':
                # all sources possibly associated with PWNe or unidentified are supposed to be PWNe
                rname, pwn_dict, distx, disty, radr, frlog = find_source_to_delete(pwn_dict,src_dir.l_deg(),src_dir.b_deg(),get_model_radius(model),1.e-2 *source['spec_flux_1TeV_crab'],radmin=radmin)
                pwn_models.remove(rname)
                n_pwn_del += 1
                pwn_distx.append(distx)
                pwn_disty.append(disty)
                pwn_radr.append(radr)
                pwn_frlog.append(frlog)
            else:
                pass

msg = 'Added {} gamma-cat sources\n'.format(len(gammacat_ids))
print(msg)
outfile.write(msg)
msg = 'Set estimated cutoffs for {} PWN, {} SNR, {} UNID\n'.format(n_ecut_pwn,n_ecut_snr,n_ecut_unid)
print(msg)
outfile.write(msg)
msg = 'Deleted {} synthetic binaries\n'.format(n_bin_del)
print(msg)
outfile.write(msg)
msg = 'Deleted {} synthetic young SNRs\n'.format(n_snr_del)
print(msg)
outfile.write(msg)
msg = 'Deleted {} synthetic interacting SNRs\n'.format(n_isnr_del)
print(msg)
outfile.write(msg)
msg = 'Deleted {} synthetic PWNe\n'.format(n_pwn_del)
print(msg)
outfile.write(msg)
msg = 'Deleted {} synthetic composite PWNe/SNRs\n'.format(n_comp_del)
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
ax5.hist((1 - np.cos(np.deg2rad(radmin))) * gammacat_flux / (1 - np.cos(np.deg2rad(np.maximum(gammacat_rads,radmin)))),
         bins=bins_lognlogs, density=False, histtype='step', cumulative=-1,
         label='gamma-cat', alpha=0.5, linewidth=2)
ax2.hist(gammacat_lons, bins=bins_lon, density=False, histtype='step',
         label='gamma-cat', alpha=0.5, linewidth=2)
ax3.hist(gammacat_lats, bins=bins_lat, density=False, histtype='step',
         label='gamma-cat', alpha=0.5, linewidth=2)
ax0.hist(gammacat_rads, bins=bins_rad, density=False, histtype='step',
         label='gamma-cat', alpha=0.5, linewidth=2)

# make distributions from gammalib model container
lons, lats, radii, fluxes, names = dist_from_gammalib(models)
# change lon range from 0...360 to -180...180
lons = np.array(lons)
lons[lons > 180] = lons[lons > 180] - 360.
ax1.hist(fluxes, bins=bins_lognlogs, density=False, histtype='step', cumulative=-1,
         label='converted gamma-cat', alpha=0.5, linewidth=2, linestyle=':')
ax5.hist((1 - np.cos(np.deg2rad(radmin))) * np.array(fluxes) / (1 - np.cos(np.deg2rad(np.maximum(radii,radmin)))),
         bins=bins_lognlogs, density=False, histtype='step', cumulative=-1,
         label='converted gamma-cat', alpha=0.5, linewidth=2, linestyle=':')
ax2.hist(lons, bins=bins_lon, density=False, histtype='step',
         label='converted gamma-cat', alpha=0.5, linewidth=2, linestyle=':')
ax3.hist(lats, bins=bins_lat, density=False, histtype='step',
         label='converted gamma-cat', alpha=0.5, linewidth=2, linestyle=':')
ax0.hist(radii, bins=bins_rad, density=False, histtype='step',
         label='converted gamma-cat', alpha=0.5, linewidth=2, linestyle=':')

# add templates

# read template list
template_list = open('../known-sources/templates/templates.dat').readlines()

replaced = 0
added = 0
# keep track of synthetic sources that are be deleted
n_bin_del = 0
n_snr_del = 0
n_isnr_del = 0
n_comp_del = 0
n_pwn_del = 0
for template in template_list:
    new = False
    if template[0] == '#':
        # skip header and commented lines
        pass
    else:
        name, id, cl = template.split(',')
        # strip leading/trailing spaces from class name
        cl = cl.strip()
        # remove gamma-cat model if present, mark source as new
        if 'NONE' in id:
            # source not included in gammacat, mark as new
            added += 1
            new = True
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
            # if model contains spatial map take care of it
            if model.spatial().type() == 'DiffuseMap' or model.spatial().type() == 'DiffuseMapCube':
                # find model map name and path
                filename = model.spatial().filename().file()
                filepath = model.spatial().filename().path()
                # copy file to output directory
                shutil.copy(filepath + filename, './')
                # replace file with the one in output directory
                if model.spatial().type() == 'DiffuseMap':
                    model.spatial(gammalib.GModelSpatialDiffuseMap(filename))
                elif model.spatial().type() == 'DiffuseMapCube':
                    model.spatial(gammalib.GModelSpatialDiffuseCube(gammalib.GFilename(filename),
                                                                    model.spatial()['Normalization'].value()))
            # if model contains spatial map take care of it
            if model.spectral().type() == 'FileFunction':
                # find spectrum file name and path
                filename = model.spectral().filename().file()
                filepath = model.spectral().filename().path()
                # copy file to output directory
                shutil.copy(filepath + filename, './')
                # replace file with the one in output directory
                model.spectral(gammalib.GModelSpectralFunc(gammalib.GFilename(filename),
                               model.spectral()['Normalization'].value()))
            # append model to container
            models.append(model)
            # if new remove synthetic source
            if new:
                src_dir = get_model_dir(model)
                src_flux = flux_Crab(model,1.,1000.)
                if cl == 'bin':
                    rname, bin_dict, distx, disty, radr, frlog = find_source_to_delete(bin_dict,src_dir.l_deg(),src_dir.b_deg(), get_model_radius(model), src_flux,radmin=radmin)
                    bin_models.remove(rname)
                    n_bin_del +=1
                    bin_distx.append(distx)
                    bin_disty.append(disty)
                    bin_radr.append(radr)
                    bin_frlog.append(frlog)
                elif cl == 'comp':
                    rname, comp_dict, distx, disty, radr, frlog = find_source_to_delete(comp_dict,src_dir.l_deg(),src_dir.b_deg(),get_model_radius(model),src_flux,radmin=radmin)
                    pwn_models.remove(rname[0])
                    snr_models.remove(rname[1])
                    n_comp_del += 1
                    comp_distx.append(distx)
                    comp_disty.append(disty)
                    comp_radr.append(radr)
                    comp_frlog.append(frlog)
                elif cl == 'isnr':
                    rname, isnr_dict, distx, disty, radr, frlog = find_source_to_delete(isnr_dict,src_dir.l_deg(),src_dir.b_deg(), get_model_radius(model),src_flux,radmin=radmin)
                    isnr_models.remove(rname)
                    n_isnr_del += 1
                    isnr_distx.append(distx)
                    isnr_disty.append(disty)
                    isnr_radr.append(radr)
                    isnr_frlog.append(frlog)
                elif cl == 'snr':
                    rname, snr_dict, distx, disty, radr, frlog = find_source_to_delete(snr_dict,src_dir.l_deg(),src_dir.b_deg(), get_model_radius(model),src_flux,radmin=radmin)
                    snr_models.remove(rname)
                    n_snr_del += 1
                    snr_distx.append(distx)
                    snr_disty.append(disty)
                    snr_radr.append(radr)
                    snr_frlog.append(frlog)
                elif cl == 'pwn' or cl == 'unid':
                    # all sources possibly associated with PWNe or unidentified are supposed to be PWNe
                    rname, pwn_dict, distx, disty, radr, frlog = find_source_to_delete(pwn_dict,src_dir.l_deg(),src_dir.b_deg(),get_model_radius(model),src_flux,radmin=radmin)
                    pwn_models.remove(rname)
                    n_pwn_del += 1
                    pwn_distx.append(distx)
                    pwn_disty.append(disty)
                    pwn_radr.append(radr)
                    pwn_frlog.append(frlog)

msg = 'Replaced {} gamma-cat sources with templates. Added {} sources as templates\n'.format(
    replaced, added)
print(msg)
outfile.write(msg)
msg = 'Deleted {} synthetic binaries\n'.format(n_bin_del)
print(msg)
outfile.write(msg)
msg = 'Deleted {} synthetic young SNRs\n'.format(n_snr_del)
print(msg)
outfile.write(msg)
msg = 'Deleted {} synthetic interacting SNRs\n'.format(n_isnr_del)
print(msg)
outfile.write(msg)
msg = 'Deleted {} synthetic PWNe\n'.format(n_pwn_del)
print(msg)
outfile.write(msg)
msg = 'Deleted {} synthetic composite PWNe/SNRs\n'.format(n_comp_del)
print(msg)
outfile.write(msg)

# re-make distributions from gammalib model container
lons, lats, radii, fluxes, names = dist_from_gammalib(models)
# change lon range from 0...360 to -180...180
lons = np.array(lons)
lons[lons > 180] = lons[lons > 180] - 360.
ax1.hist(fluxes, bins=bins_lognlogs, density=False, histtype='step', cumulative=-1,
         label='gamma-cat + templates', alpha=0.5, linewidth=2, linestyle=':')
ax5.hist((1 - np.cos(np.deg2rad(radmin))) * np.array(fluxes) / (1 - np.cos(np.deg2rad(np.maximum(radii,radmin)))),
         bins=bins_lognlogs, density=False, histtype='step', cumulative=-1,
         label='gamma-cat + templates', alpha=0.5, linewidth=2, linestyle=':')
ax2.hist(lons, bins=bins_lon, density=False, histtype='step',
         label='gamma-cat + templates', alpha=0.5, linewidth=2, linestyle=':')
ax3.hist(lats, bins=bins_lat, density=False, histtype='step',
         label='gamma-cat + templates', alpha=0.5, linewidth=2, linestyle=':')
ax0.hist(radii, bins=bins_rad, density=False, histtype='step',
         label='gamma-cat + templates', alpha=0.5, linewidth=2, linestyle=':')

# add binaries

# read binary list
binary_list = open('../bin/binaries.dat').readlines()

# read models
binary_models = gammalib.GModels('../bin/models_binaries.xml')

replaced = 0
added = 0
# keep track of synthetic sources that are being deleted
n_bin_del = 0
for binary in binary_list:
    if binary[0] == '#':
        # skip header and commented lines
        pass
    else:
        newsrc = False
        id, bmodels = binary.split(',')
        if 'None' in id:
            # source not included in gammacat, pass
            added += 1
            newsrc = True
        else:
            # remove gamma-cat model if present
            id = int(id)
            if id in gammacat_ids:
                source = gammacat[gammacat['source_id'] == id][0]
                models.remove(source['common_name'])
                replaced += 1
            else:
                added += 1
                newsrc = True
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
        # get rid of a synthetic binary for each newly added one
        if newsrc:
            src_dir = get_model_dir(model)
            flux = 0.
            for s in range(len(model_names)):
                flux += flux_Crab(models[-1-s],1.,1000.)
            rname, bin_dict, distx, disty, radr, frlog = find_source_to_delete(bin_dict, src_dir.l_deg(), src_dir.b_deg(), get_model_radius(model),
                                                    flux,radmin=radmin)
            bin_models.remove(rname)
            n_bin_del += 1
            bin_distx.append(distx)
            bin_disty.append(disty)
            bin_radr.append(radr)
            bin_frlog.append(frlog)


msg = 'Replaced {} gamma-cat sources with binaries. Added {} sources as binaries\n'.format(
    replaced, added)
print(msg)
outfile.write(msg)
msg = 'Deleted {} synthetic binaries\n'.format(n_bin_del)
print(msg)
outfile.write(msg)

# add pulsars

# read models
psr_models = gammalib.GModels('../psr/psrs_gps_models.xml')
# list of psr to delete from gamma-cat
psr_del = np.genfromtxt('../psr/del_psr.txt',delimiter=',',dtype=None,names=True)

replaced = 0
added = 0
for psr in psr_models:
        # check if pulsar needs to replace an existing gamma-cat model
        if psr.name() in psr_del['name'].astype('str'):
            # remove gamma-cat model
            id = psr_del['gammacatid'][psr_del['name'].astype('str') == psr.name()][0]
            if id in gammacat_ids:
                source = gammacat[gammacat['source_id'] == id][0]
                models.remove(source['common_name'])
                replaced += 1
            else:
                added += 1
        else:
            added += 1
        # find phasecurve name and path
        filename = psr.temporal().filename().file()
        filepath = psr.temporal().filename().path()
        # copy file to output directory
        shutil.copy(filepath + filename, './')
        # replace file with the one in output directory
        psr.temporal().filename(filename)
        models.append(psr)


msg = 'Replaced {} gamma-cat sources with pulsars. Added {} sources as pulsars\n'.format(
    replaced, added)
print(msg)
outfile.write(msg)

# re-make distributions from gammalib model container
lons, lats, radii, fluxes, names = dist_from_gammalib(models)
# change lon range from 0...360 to -180...180
lons = np.array(lons)
lons[lons > 180] = lons[lons > 180] - 360.
ax1.hist(fluxes, bins=bins_lognlogs, density=False, histtype='step', cumulative=-1,
         label='gamma-cat + templates + bin + psr', alpha=0.5, linewidth=2, linestyle=':')
ax5.hist((1 - np.cos(np.deg2rad(radmin))) * np.array(fluxes) / (1 - np.cos(np.deg2rad(np.maximum(radii,radmin)))),
         bins=bins_lognlogs, density=False, histtype='step', cumulative=-1,
         label='gamma-cat + templates + bin + psr', alpha=0.5, linewidth=2, linestyle=':')
ax2.hist(lons, bins=bins_lon, density=False, histtype='step',
         label='gamma-cat + templates + bin + psr', alpha=0.5, linewidth=2, linestyle=':')
ax3.hist(lats, bins=bins_lat, density=False, histtype='step',
         label='gamma-cat + templates + bin + psr', alpha=0.5, linewidth=2, linestyle=':')
ax0.hist(radii, bins=bins_rad, density=False, histtype='step',
         label='gamma-cat + templates + bin + psr', alpha=0.5, linewidth=2, linestyle=':')

# add 3FHL

result_fhl = append_fhl(models,bmax,
                        bin_models, bin_dict, bin_distx, bin_disty, bin_radr, bin_frlog,
                        snr_models, snr_dict, snr_distx, snr_disty, snr_radr, snr_frlog,
                        isnr_models, isnr_dict, isnr_distx, isnr_disty, isnr_radr, isnr_frlog,
                        pwn_models, pwn_dict, pwn_distx, pwn_disty, pwn_radr, pwn_frlog,
                        comp_dict, comp_distx, comp_disty, comp_radr, comp_frlog,
                        dist_sigma=3.)
models = result_fhl['models']
msg = 'Added {} FHL sources, of which {} as pointlike and {} as extended.\n'.format(
    result_fhl['newpt']+result_fhl['newext'], result_fhl['newpt'],result_fhl['newext'])
print(msg)
outfile.write(msg)

ecut_pwn.extend(result_fhl['ecut_pwn'])
ecut_snr.extend(result_fhl['ecut_snr'])
ecut_agn.extend(result_fhl['ecut_agn'])
ecut_unid.extend(result_fhl['ecut_unid'])
msg = 'Set estimated cutoffs for {} PWN, {} SNR, {} AGN, {} UNID\n'.format(result_fhl['n_ecut_pwn'],result_fhl['n_ecut_snr'],result_fhl['n_ecut_agn'],result_fhl['n_ecut_unid'])
print(msg)
outfile.write(msg)

bin_models = result_fhl['bin_models']
bin_dict = result_fhl['bin_dict']
bin_distx = result_fhl['bin_distx']
bin_disty = result_fhl['bin_disty']
bin_radr = result_fhl['bin_radr']
bin_frlog = result_fhl['bin_frlog']
msg = 'Deleted {} synthetic binaries\n'.format(result_fhl['n_bin_del'])
print(msg)
outfile.write(msg)

snr_models = result_fhl['snr_models']
snr_dict = result_fhl['snr_dict']
snr_distx = result_fhl['snr_distx']
snr_disty = result_fhl['snr_disty']
snr_radr = result_fhl['snr_radr']
snr_frlog = result_fhl['snr_frlog']
msg = 'Deleted {} synthetic young SNRs\n'.format(result_fhl['n_snr_del'])
print(msg)
outfile.write(msg)

isnr_models = result_fhl['isnr_models']
isnr_dict = result_fhl['isnr_dict']
isnr_distx = result_fhl['isnr_distx']
isnr_disty = result_fhl['isnr_disty']
isnr_radr = result_fhl['isnr_radr']
isnr_frlog = result_fhl['isnr_frlog']
msg = 'Deleted {} synthetic interacting SNRs\n'.format(result_fhl['n_isnr_del'])
print(msg)
outfile.write(msg)

pwn_models = result_fhl['pwn_models']
pwn_dict = result_fhl['pwn_dict']
pwn_distx = result_fhl['pwn_distx']
pwn_disty = result_fhl['pwn_disty']
pwn_radr = result_fhl['pwn_radr']
pwn_frlog = result_fhl['pwn_frlog']
msg = 'Deleted {} synthetic PWNe\n'.format(result_fhl['n_pwn_del'])
print(msg)
outfile.write(msg)

comp_dict = result_fhl['comp_dict']
comp_distx = result_fhl['comp_distx']
comp_disty = result_fhl['comp_disty']
comp_radr = result_fhl['comp_radr']
comp_frlog = result_fhl['comp_frlog']
msg = 'Deleted {} synthetic composite PWNe/SNRs\n'.format(result_fhl['n_comp_del'])
print(msg)
outfile.write(msg)

if len(result_fhl['msg']) > 0:
    print(result_fhl['msg'])
    outfile.write(result_fhl['msg'])

# add HAWC

models, newpt, newext, warning, pwn_models, n_pwn_del, pwn_distx, pwn_disty, pwn_radr, pwn_frlog = append_hawc(models,bmax,
                                                                                                               pwn_models, pwn_dict, pwn_distx, pwn_disty, pwn_radr, pwn_frlog,
                                                                                                               dist_sigma=3,radmin=radmin)

msg = 'Added {} HAWC sources, of which {} as pointlike and {} as extended.\n'.format(
    newpt+newext, newpt,newext)
print(msg)
outfile.write(msg)

msg = 'Deleted {} synthetic PWNe\n'.format(n_pwn_del)
print(msg)
outfile.write(msg)

if len(warning) > 0:
    print(warning)
    outfile.write(warning)

# save cutoff values
np.save('ecut_pwn.npy',ecut_pwn)
np.save('ecut_snr.npy',ecut_snr)
np.save('ecut_agn.npy',ecut_agn)
np.save('ecut_unid.npy',ecut_unid)

# re-make distributions from gammalib model container
lons, lats, radii, fluxes, names = dist_from_gammalib(models)
# change lon range from 0...360 to -180...180
lons = np.array(lons)
lons[lons > 180] = lons[lons > 180] - 360.
ax1.hist(fluxes, bins=bins_lognlogs, density=False, histtype='step', cumulative=-1,
         label='gamma-cat + templates + bin + psr + FHL + HAWC', alpha=0.5, linewidth=2, linestyle=':')
ax5.hist((1 - np.cos(np.deg2rad(radmin))) * np.array(fluxes) / (1 - np.cos(np.deg2rad(np.maximum(radii,radmin)))),
         bins=bins_lognlogs, density=False, histtype='step', cumulative=-1,
         label='gamma-cat + templates + bin + psr + FHL + HAWC', alpha=0.5, linewidth=2, linestyle=':')
ax2.hist(lons, bins=bins_lon, density=False, histtype='step',
         label='gamma-cat + templates + bin + psr + FHL + HAWC', alpha=0.5, linewidth=2, linestyle=':')
ax3.hist(lats, bins=bins_lat, density=False, histtype='step',
         label='gamma-cat + templates + bin + psr + FHL + HAWC', alpha=0.5, linewidth=2, linestyle=':')
ax0.hist(radii, bins=bins_rad, density=False, histtype='step',
         label='gamma-cat + templates + bin + psr + FHL + HAWC', alpha=0.5, linewidth=2, linestyle=':')

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

# add synthetic binaries, PWNe and SNRs

# binaries
for model in bin_models:
    models.append(model)
    models_syn.append(model)

msg = 'Added {} synthetic binaries\n'.format(bin_models.size())
print(msg)
outfile.write(msg)

# re-make distributions from gammalib model container
lons, lats, radii, fluxes, names = dist_from_gammalib(models)
# change lon range from 0...360 to -180...180
lons = np.array(lons)
lons[lons > 180] = lons[lons > 180] - 360.
ax1.hist(fluxes, bins=bins_lognlogs, density=False, histtype='step', cumulative=-1,
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin', alpha=0.5, linewidth=2, linestyle=':')
ax5.hist((1 - np.cos(np.deg2rad(radmin))) * np.array(fluxes) / (1 - np.cos(np.deg2rad(np.maximum(radii,radmin)))),
         bins=bins_lognlogs, density=False, histtype='step', cumulative=-1, linestyle=':',
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin', alpha=0.5, linewidth=2)
ax2.hist(lons, bins=bins_lon, density=False, histtype='step',
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin', alpha=0.5, linewidth=2, linestyle=':')
ax3.hist(lats, bins=bins_lat, density=False, histtype='step',
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin', alpha=0.5, linewidth=2, linestyle=':')
ax0.hist(radii, bins=bins_rad, density=False, histtype='step',
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin', alpha=0.5, linewidth=2, linestyle=':')

# young SNRs
for model in snr_models:
    models.append(model)
    models_syn.append(model)

msg = 'Added {} synthetic young SNRs, of which {} in composite systems\n'.format(snr_models.size(),len(comp_dict['name']))
print(msg)
outfile.write(msg)

# re-make distributions from gammalib model container
lons, lats, radii, fluxes, names = dist_from_gammalib(models)
# change lon range from 0...360 to -180...180
lons = np.array(lons)
lons[lons > 180] = lons[lons > 180] - 360.
ax1.hist(fluxes, bins=bins_lognlogs, density=False, histtype='step', cumulative=-1,
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin + synth SNR', alpha=0.5, linewidth=2, linestyle=':')
ax5.hist((1 - np.cos(np.deg2rad(radmin))) * np.array(fluxes) / (1 - np.cos(np.deg2rad(np.maximum(radii,radmin)))),
         bins=bins_lognlogs, density=False, histtype='step', cumulative=-1, linestyle=':',
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin + synth SNR', alpha=0.5, linewidth=2)
ax2.hist(lons, bins=bins_lon, density=False, histtype='step',
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin + synth SNR', alpha=0.5, linewidth=2, linestyle=':')
ax3.hist(lats, bins=bins_lat, density=False, histtype='step',
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin + synth SNR', alpha=0.5, linewidth=2, linestyle=':')
ax0.hist(radii, bins=bins_rad, density=False, histtype='step',
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin + synth SNR', alpha=0.5, linewidth=2, linestyle=':')

# PWNe
for model in pwn_models:
    models.append(model)
    models_syn.append(model)

msg = 'Added {} synthetic PWNe, of which {} in composite systems\n'.format(pwn_models.size(),len(comp_dict['name']))
print(msg)
outfile.write(msg)

# re-make distributions from gammalib model container
lons, lats, radii, fluxes, names = dist_from_gammalib(models)
# change lon range from 0...360 to -180...180
lons = np.array(lons)
lons[lons > 180] = lons[lons > 180] - 360.
ax1.hist(fluxes, bins=bins_lognlogs, density=False, histtype='step', cumulative=-1,
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin + synth SNR + synth PWNe', alpha=0.5, linewidth=2, linestyle=':')
ax5.hist((1 - np.cos(np.deg2rad(radmin))) * np.array(fluxes) / (1 - np.cos(np.deg2rad(np.maximum(radii,radmin)))),
         bins=bins_lognlogs, density=False, histtype='step', cumulative=-1, linestyle=':',
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin + synth SNR + synth PWNe', alpha=0.5, linewidth=2)
ax2.hist(lons, bins=bins_lon, density=False, histtype='step',
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin + synth SNR + synth PWNe', alpha=0.5, linewidth=2, linestyle=':')
ax3.hist(lats, bins=bins_lat, density=False, histtype='step',
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin + synth SNR + synth PWNe', alpha=0.5, linewidth=2, linestyle=':')
ax0.hist(radii, bins=bins_rad, density=False, histtype='step',
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin + synth SNR + synth PWNe', alpha=0.5, linewidth=2, linestyle=':')


# interacting SNRs
for model in isnr_models:
    if model.spectral().type() == 'FileFunction':
        # find spectrum file name and path
        filename = model.spectral().filename().file()
        filepath = model.spectral().filename().path()
        # copy file to output directory
        shutil.copy(filepath + filename, './')
        # replace file with the one in output directory
        model.spectral(gammalib.GModelSpectralFunc(gammalib.GFilename(filename), model.spectral()['Normalization'].value()))
    models.append(model)
    models_syn.append(model)

msg = 'Added {} synthetic interacting SNRs\n'.format(isnr_models.size())
print(msg)
outfile.write(msg)

# re-make distributions from gammalib model container
lons, lats, radii, fluxes, names = dist_from_gammalib(models)
# change lon range from 0...360 to -180...180
lons = np.array(lons)
lons[lons > 180] = lons[lons > 180] - 360.
ax1.hist(fluxes, bins=bins_lognlogs, density=False, histtype='step', cumulative=-1,
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin + synth SNR + synth PWNe + synth iSNR', alpha=0.5, linewidth=2, linestyle=':')
ax5.hist((1 - np.cos(np.deg2rad(radmin))) * np.array(fluxes) / (1 - np.cos(np.deg2rad(np.maximum(radii,radmin)))),
         bins=bins_lognlogs, density=False, histtype='step', cumulative=-1,
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin + synth SNR + synth PWNe + synth iSNR', alpha=0.5, linewidth=2, linestyle=':')
ax2.hist(lons, bins=bins_lon, density=False, histtype='step',
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin + synth SNR + synth PWNe + synth iSNR', alpha=0.5, linewidth=2, linestyle=':')
ax3.hist(lats, bins=bins_lat, density=False, histtype='step',
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin + synth SNR + synth PWNe + synth iSNR', alpha=0.5, linewidth=2, linestyle=':')
ax0.hist(radii, bins=bins_rad, density=False, histtype='step',
         label='gamma-cat + templates + bin + psr + FHL + HAWC + synth bin + synth SNR + synth PWNe + synth iSNR', alpha=0.5, linewidth=2, linestyle=':')

# add IEM

# read template list
component_list = open('../iem/components.dat').readlines()
ncomp = 0

for name in component_list:
    if name[0] == '#':
        # skip header and commented lines
        pass
    else:
        # add component
        name = name.split('\n')[0]
        models_template = gammalib.GModels('../iem/{}.xml'.format(name))
        for model in models_template:
            # if model contains spatial map take care of it
            if model.spatial().type() == 'DiffuseMap':
                # find model map name and path
                filename = model.spatial().filename().file()
                filepath = model.spatial().filename().path()
                # copy file to output directory
                shutil.copy(filepath + filename, './')
                # replace file with the one in output directory
                model.spatial(gammalib.GModelSpatialDiffuseMap(filename))
            # append model to container
            models.append(model)
            ncomp +=1

msg = 'Added {} interstellar emission components\n'.format(ncomp)
print(msg)
outfile.write(msg)

# CHECKS for problems with spectral models
for model in models:
    flux = flux_Crab(model,1.,1000.)
    if flux > 1.:
        msg = 'WARNING: source {} has a flux of {} Crab\n'.format(model.name(), flux)
        print(msg)
        outfile.write(msg)

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
models_syn.save('models_gps_synthetic.xml')

# close report file
outfile.close()

# save diagnostic plots
ax1.legend(fontsize=5,loc='lower left')
fig1.savefig('logNlogS.png', dpi=300)
ax2.legend(fontsize=5)
ax2.set_xlim(180, -180)
fig2.savefig('glon.png', dpi=300)
ax3.legend(fontsize=5)
fig3.savefig('glat.png', dpi=300)
ax0.legend(fontsize=5)
fig0.savefig('radius.png', dpi=300)
ax5.legend(fontsize=5)
fig5.savefig('logNlogS-solidang.png', dpi=300)

# make more diagnostic plots with properties of synthetic sources deleted
plot_del_sources(bin_distx,bin_disty,bin_radr,bin_frlog,'bin','binaries')
plot_del_sources(snr_distx,snr_disty,snr_radr,snr_frlog,'snr','young SNRs')
plot_del_sources(pwn_distx,pwn_disty,pwn_radr,pwn_frlog,'pwn','PWNe')
plot_del_sources(comp_distx,comp_disty,comp_radr,comp_frlog,'comp','composite PWN/SNR')
plot_del_sources(isnr_distx,isnr_disty,isnr_radr,isnr_frlog,'isnr','interacting SNRs')