import gammalib
import numpy as np
from astropy.io import fits
from utils import *
from cutoffs import get_cutoff

# retrieve Fermi high-energy catalog data
fhl = fits.getdata('../known-sources/external-input/gll_psch_v13.fit', 1)
ext_fhl = fits.getdata('../known-sources/external-input/gll_psch_v13.fit', 2)

# table of interacting SNRs
snr_class = fits.getdata('../known-sources/external-input/iSNRin3FHL.fits',1)

# filter sources wth known TeV association, they are included elsewhere
# this avoids issues due to model building choices (i.e., numerical precision, source for position etc.)
fhl = fhl[fhl['ASSOC_TEV'] == ' ']


def append_fhl(models, bmax,
                bin_models, bin_dict, bin_distx, bin_disty, bin_radr, bin_frlog,
                snr_models, snr_dict, snr_distx, snr_disty, snr_radr, snr_frlog,
                isnr_models, isnr_dict, isnr_distx, isnr_disty, isnr_radr, isnr_frlog,
                dist_sigma=3., sig50_thresh=3., eph_thresh=100.):
    """
    Append missing models from Fermi high-energy catalog to gammalib model container
    :param models: ~gammalib.GModels, gammalib model container
    :param bmax: float, maximum latitude of models to include
    :param dist_sigma: float, minimum distance in sigma's from pre-existing source to add as new
    :param sig50_thresh: float, threshold sigma on significance above 50 GeV to apply
    :param eph_thresh: float, minimum energy of detected photons required
    :return: result: dictionary
    """
    # filter FHL table based on latitude and hardness
    # hard sources are seleted based on significance above 50 GeV and highest photon energy
    # calculate significance above 50 GeV
    sig50 = np.sqrt(np.sum(np.power(fhl['Sqrt_TS_Band'][:,2:], 2),axis=1))
    # filter
    fsources = fhl[(np.abs(fhl['GLAT']) <= bmax) & (sig50 >= sig50_thresh) & (
                fhl['HEP_Energy'] >= eph_thresh)]

    # prepare new gammalib model container
    newpt = 0
    newext = 0
    n_bin_del = 0
    n_snr_del = 0
    n_isnr_del = 0
    newmodels = gammalib.GModels()
    # keep track also of artificial cutoffs
    ecut_pwn = []
    ecut_snr = []
    ecut_unid = []
    ecut_agn = []
    n_ecut_pwn = 0
    n_ecut_snr = 0
    n_ecut_agn = 0
    n_ecut_unid = 0

    msg = ''

    # loop over fermi sources
    for fsource in fsources:
        # assume source is new
        new = 1
        # fermi source position as gammalib.GSkyDir
        fdir = gammalib.GSkyDir()
        ra = np.double(fsource['RAJ2000'])
        dec = np.double(fsource['DEJ2000'])
        fdir.radec_deg(ra, dec)
        ######################################### treat pointlike sources for Fermi ###########
        if fsource['Extended_Source_Name'] == '':
            # case of sources pointlike for Fermi
            # loop over gammalib container and determine closest neighbor
            # initialize at 1000
            dist_min = 1000.
            for source in models:
                dir = get_model_dir(source)
                dist = fdir.dist_deg(dir)
                dist_min = np.minimum(dist, dist_min)
            if dist_min < dist_sigma * fsource['Conf_95_SemiMajor'] / 2:
                # closeby source foud in container, source will not be added
                new = 0
            else:
                # source will be added to container, set spatial model
                spatial = gammalib.GModelSpatialPointSource(fdir)
                newpt += 1
        ######################################### treat extended sources for Fermi ############
        else:
            # retrieve Fermi extended source radius
            ext_fsource = ext_fhl[ext_fhl['Source_Name'] == fsource['Extended_Source_Name']][0]
            fradius = np.double(ext_fsource['Model_SemiMajor'])
            # "corrected" distance, i.e, centre distance - radius of source
            # initialize at 1000
            dist_corr = 1000.
            for source in models:
                dir = get_model_dir(source)
                dist = fdir.dist_deg(dir)
                # get source radius according to model used
                radius = get_model_radius(source)
                # use as radius the max between Fermi and model
                radius = np.maximum(fradius,radius)
                # subtract from distance the max radius
                dist -= radius
                dist_corr = np.minimum(dist,dist_corr)
            if dist_corr < 0.:
                # overlapping source foud in container, source will not be added
                new = 0
            else:
                # source will be added to container, set spatial model
                fradius2 = np.double(ext_fsource['Model_SemiMinor'])
                fpangle = np.double(ext_fsource['Model_PosAng'])
                # retrieve Fermi spatial model to set it in container
                if ext_fsource['Model_Form'] == 'Disk':
                    if fradius2 == fradius:
                        spatial = gammalib.GModelSpatialRadialDisk(fdir,fradius)
                    else:
                        spatial = gammalib.GModelSpatialEllipticalDisk(fdir, fradius, fradius2, fpangle)
                elif ext_fsource['Model_Form'] == '2D Gaussian':
                    if fradius2 == fradius:
                        spatial = gammalib.GModelSpatialRadialGauss(fdir,fradius)
                    else:
                        spatial = gammalib.GModelSpatialEllipticalGauss(fdir, fradius, fradius2, fpangle)
                elif ext_fsource['Model_Form'] == 'Ring':
                    spatial = gammalib.GModelSpatialRadialShell(fdir,fradius, fradius2)
                elif ext_fsource['Model_Form'] == 'Map':
                    print('{} modeled by spatial template, which is not implemented, skip'.format(fsource['Source_Name']))
                    new = 0
                else:
                    print('{} modeled by model type {}, which is not implemented, skip'.format(fsource['Source_Name'],ext_fsource['Model_Form']))
                    new = 0
                if new == 1:
                    newext +=1
        ######################################### spectra #####################################
        if new == 1:
            # spectral model
            eref = gammalib.GEnergy(np.double(fsource['Pivot_Energy']), 'GeV')
            norm = np.double(fsource['Flux_Density'])
            norm *= 1.e-3  # ph GeV-1 -> ph MeV-1
            # use power law model only if signifcance of curvature < 1
            # this avoids extrapolating hard power laws not justified by the data
            # some source with low significant curvature have beta < 0, they will be modelled as power laws
            if fsource['Signif_Curve'] < 1 or fsource['beta'] < 0.:
                index = np.double(fsource['PowerLaw_Index'])
                if index < 2.4:
                    # correction for fake pevatrons
                    if fsource['CLASS'] == 'PWN' or fsource['CLASS'] == 'pwn':
                        # dummy model to obtain search radius
                        mod = gammalib.GModelSky(spatial, gammalib.GModelSpectralPlaw())
                        # search radius
                        rad = get_model_radius(mod) + 0.2
                        # set cutoff
                        ecut = get_cutoff(ra, dec, 'PSR', rad_search=rad)
                        ecut_pwn.append(ecut)
                        n_ecut_pwn += 1
                    elif fsource['CLASS'] == 'SNR' or fsource['CLASS'] == 'snr':
                        gname = fsource['ASSOC1']
                        if 'SNR' in gname:
                            pass
                        else:
                            gname = None
                        # compute cutoff
                        # we use default value for particle spectral index because measurements are highly uncertain
                        ecut = get_cutoff(ra, dec, 'SNR', name=gname)
                        ecut_snr.append(ecut)
                        n_ecut_snr += 1
                    elif fsource['CLASS'] == 'bll' or fsource['CLASS'] == 'bcu' or fsource['CLASS'] == 'fsrq':
                        ecut = get_cutoff(ra, dec, 'AGN', z = fsource['Redshift'])
                        ecut_agn.append(ecut)
                        n_ecut_agn += 1
                    else:
                        if fsource['CLASS'] == '' or fsource['CLASS'] == 'unknown':
                            pass
                        else:
                            # set warning if we have hard source of unexpected type
                            msg += 'FHL source {} of type {} has an unxepctedly hard spectrum ' \
                                  'with index {}. We are setting a random artificial cutoff\n'.format(
                                fsource['Source_Name'], fsource['CLASS'], index)
                            print(msg)
                        ecut = get_cutoff(ra, dec, 'UNID')
                        ecut_unid.append(ecut)
                        n_ecut_unid += 1
                    ecut = gammalib.GEnergy(np.double(ecut), 'TeV')
                    spectral = gammalib.GModelSpectralExpPlaw(norm, -index, eref, ecut)
                else:
                    spectral = gammalib.GModelSpectralPlaw(norm, -index, eref)
            else:
                index = np.double(fsource['Spectral_Index'])
                curvature = np.double(fsource['beta'])
                spectral = gammalib.GModelSpectralLogParabola(norm, -index, eref, -curvature)
            # assemble model and append to container
            model = gammalib.GModelSky(spatial, spectral)
            model.name(fsource['Source_Name'])
            newmodels.append(model)
            # delete synthetic sources as needed
            if fsource['CLASS'] == 'HMB' or fsource['CLASS'] == 'hmb' or fsource['CLASS'] == 'BIN' or fsource['CLASS'] == 'bin':
                rname, bin_dict, distx, disty, radr, frlog = find_source_to_delete(bin_dict, fdir.l_deg(),
                                                        fdir.b_deg(), get_model_radius(model),
                                                        flux_Crab(model,1.,1000.))
                bin_models.remove(rname)
                n_bin_del += 1
                bin_distx.append(distx)
                bin_disty.append(disty)
                bin_radr.append(radr)
                bin_frlog.append(frlog)
            elif fsource['CLASS'] == 'PWN' or fsource['CLASS'] == 'pwn' or fsource['CLASS'] == 'spp':
                # implement synthetic PWNe here
                pass
            elif fsource['CLASS'] == 'SNR' or fsource['CLASS'] == 'snr':
                # determine if snr is interacting
                assoc = snr_class[snr_class['ASSOC1'] == fsource['ASSOC1']]
                if assoc[0]['is int'] == 'yes':
                    # interacting
                    rname, isnr_dict, distx, disty, radr, frlog = find_source_to_delete(isnr_dict,
                                                                                 fdir.l_deg(),
                                                                                 fdir.b_deg(), get_model_radius(model),
                                                                                 flux_Crab(model, 1.,1000.))
                    isnr_models.remove(rname)
                    n_isnr_del += 1
                    isnr_distx.append(distx)
                    isnr_disty.append(disty)
                    isnr_radr.append(radr)
                    isnr_frlog.append(frlog)
                else:
                    # young
                    rname, snr_dict, distx, disty, radr, frlog = find_source_to_delete(snr_dict,
                                                                                  fdir.l_deg(),
                                                                                  fdir.b_deg(), get_model_radius(model),
                                                                                  flux_Crab(model,1.,1000.))
                    snr_models.remove(rname)
                    n_snr_del += 1
                    snr_distx.append(distx)
                    snr_disty.append(disty)
                    snr_radr.append(radr)
                    snr_frlog.append(frlog)
            else:
                pass
        else:
            # source already present, skip
            pass

    # add new models to original container
    for model in newmodels:
        models.append(model)

    # assemble output in dictionary
    result = { 'models' : models, 'newpt' : newpt, 'newext' : newext,
               'ecut_pwn' : ecut_pwn, 'ecut_snr' : ecut_snr, 'ecut_agn' : ecut_agn,
               'ecut_unid' : ecut_unid, 'n_ecut_pwn' : n_ecut_pwn, 'n_ecut_snr' : n_ecut_snr,
               'n_ecut_agn' : n_ecut_agn, 'n_ecut_unid' : n_ecut_unid, 'msg' : msg,
               'bin_models' : bin_models, 'bin_dict' : bin_dict, 'n_bin_del' : n_bin_del,
               'bin_distx' : bin_distx, 'bin_disty' : bin_disty, 'bin_radr' : bin_radr, 'bin_frlog' : bin_frlog,
               'snr_models': snr_models, 'snr_dict': snr_dict, 'n_snr_del': n_snr_del,
               'snr_distx': snr_distx, 'snr_disty': snr_disty, 'snr_radr': snr_radr, 'snr_frlog': snr_frlog,
               'isnr_models': isnr_models, 'isnr_dict': isnr_dict, 'n_isnr_del': n_isnr_del,
               'isnr_distx': isnr_distx, 'isnr_disty': isnr_disty, 'isnr_radr': isnr_radr, 'isnr_frlog': isnr_frlog
    }

    return result
