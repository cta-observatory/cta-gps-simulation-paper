import gammalib
import numpy as np
from astropy.io import fits
from utils import *

# retrieve Fermi high-energy catalog data
fhl = fits.getdata('../known-sources/external-input/gll_psch_v13.fit', 1)
ext_fhl = fits.getdata('../known-sources/external-input/gll_psch_v13.fit', 2)

# filter sources wth known TeV association, they are included elsewhere
# this avoids issues due to model building choices (i.e., numerical precision, source for position etc.)
fhl = fhl[fhl['ASSOC_TEV'] == ' ']


def append_fhl(models, bmax, dist_sigma=3., sig50_thresh=3., eph_thresh=100.):
    """
    Append missing models from Fermi high-energy catalog to gammalib model container
    :param models: ~gammalib.GModels, gammalib model container
    :param bmax: float, maximum latitude of models to include
    :param dist_sigma: float, minimum distance in sigma's from pre-existing source to add as new
    :param sig50_thresh: float, threshold sigma on significance above 50 GeV to apply
    :param eph_thresh: float, minimum energy of detected photons required
    :return: models: ~gammalib.GModels, gammalib model container with new sources added
    :return: newpt: int, number of new pointlike sources
    :return: neext: int, number of new extended sources
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
    newmodels = gammalib.GModels()

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
            if fsource['Signif_Curve'] < 1:
                index = np.double(fsource['PowerLaw_Index'])
                spectral = gammalib.GModelSpectralPlaw(norm, -index, eref)
            else:
                index = np.double(fsource['Spectral_Index'])
                curvature = np.double(fsource['beta'])
                spectral = gammalib.GModelSpectralLogParabola(norm, -index, eref, curvature)
            # assemble model and append to container
            model = gammalib.GModelSky(spatial, spectral)
            model.name(fsource['Source_Name'])
            newmodels.append(model)
        else:
            # source already present, skip
            pass

    # add new models to original container
    for model in newmodels:
        models.append(model)

    return models, newpt, newext
