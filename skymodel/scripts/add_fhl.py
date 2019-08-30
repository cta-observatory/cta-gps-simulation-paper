import gammalib
import numpy as np
from astropy.io import fits
from utils import *

# retrieve Fermi high-energy catalog data
fhl = fits.getdata('../known-sources/external-input/gll_psch_v13.fit', 1)

# filter sources wth known TeV association, they are included elsewhere
# this avoids issues due to model building choices
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
    :return: nsrc: int, number of new sources
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
        if fsource['Extended_Source_Name'] == '':
            # case of sources pointlike for Fermi
            # fermi source position as gammalib.GSkyDir
            fdir = gammalib.GSkyDir()
            ra = np.double(fsource['RAJ2000'])
            dec = np.double(fsource['DEJ2000'])
            fdir.radec_deg(ra, dec)
            # loop over gammalib container and determine closest neighbor
            dist_min = 1000.
            for source in models:
                dir = get_model_dir(source)
                dist = fdir.dist_deg(dir)
                dist_min = np.minimum(dist, dist_min)
            if dist_min < dist_sigma * fsource['Conf_95_SemiMajor'] / 2:
                # closeby source foud in container, source will not be added
                new = 0
            else:
                # sources will be added to container, set spatial model
                src_dir = gammalib.GSkyDir()
                ra = np.double(fsource['RAJ2000'])
                dec = np.double(fsource['DEJ2000'])
                src_dir.radec_deg(ra, dec)
                spatial = gammalib.GModelSpatialPointSource(src_dir)
                newpt += 1
        else:
            # for the moment ignore extended sources
            new = 0
        if new == 1:
            # spectral model
            eref = gammalib.GEnergy(np.double(fsource['Pivot_Energy']), 'GeV')
            norm = np.double(fsource['Flux_Density'])
            norm *= 1.e-3  # ph GeV-1 -> ph MeV-1
            if fsource['SpectrumType'] == 'PowerLaw':
                index = np.double(fsource['PowerLaw_Index'])
                spectral = gammalib.GModelSpectralPlaw(norm, -index, eref)
            elif fsource['SpectrumType'] == 'LogParabola':
                index = np.double(fsource['Spectral_Index'])
                curvature = np.double(fsource['beta'])
                spectral = gammalib.GModelSpectralLogParabola(norm, -index, eref, curvature)
            else:
                print(
                    'WARNING: source {} has spectral model of type {} which is not implemented'.format(
                        fsource['Source_Name'], fsource['SpectrumType']))
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
