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
        # fermi source position as gammalib.GSkyDir
        fdir = gammalib.GSkyDir()
        ra = np.double(fsource['RAJ2000'])
        dec = np.double(fsource['DEJ2000'])
        fdir.radec_deg(ra, dec)
        if fsource['Extended_Source_Name'] == '':
            # case of sources pointlike for Fermi
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
                # source will be added to container, set spatial model
                src_dir = gammalib.GSkyDir()
                ra = np.double(fsource['RAJ2000'])
                dec = np.double(fsource['DEJ2000'])
                src_dir.radec_deg(ra, dec)
                spatial = gammalib.GModelSpatialPointSource(src_dir)
                newpt += 1
        else:
            # retrieve Fermi extended source radius
            ext_fsource = ext_fhl[ext_fhl['Source_Name'] == fsource['Extended_Source_Name']][0]
            fradius = np.double(ext_fsource['Model_SemiMajor'])
            dist_corr = 10.
            for source in models:
                dir = get_model_dir(source)
                dist = fdir.dist_deg(dir)
                if source.spatial().type() == 'PointSource':
                    radius = 0.
                elif source.spatial().type() == 'RadialGaussian':
                    radius = 2 * source['Sigma'].value()
                elif source.spatial().type() == 'EllipticalGaussian':
                    radius = 2 * source['MajorRadius'].value()
                elif source.spatial().type() == 'RadialShell':
                    radius = source['Radius'].value() + source['Width'].value()
                elif source.spatial().type() == 'RadialDisk':
                    radius = source['Radius'].value()
                elif source.spatial().type() == 'DiffuseMap':
                    radius = source.spatial().region().radius()
                else:
                    print('source {} has spatial model type {} which is not implemented'.format(source.name(),source.spatial().type()))
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
                if ext_fsource['Model_Form'] == 'Disk':
                    if fpangle == 0.:
                        spatial = gammalib.GModelSpatialRadialDisk(fdir,fradius)
                    else:
                        spatial = gammalib.GModelSpatialEllipticalDisk(fdir, fradius, fradius2, fpangle)
                elif ext_fsource['Model_Form'] == '2D Gaussian':
                    if fpangle == 0.:
                        spatial = gammalib.GModelSpatialRadialGauss(fdir,fradius)
                    else:
                        spatial = gammalib.GModelSpatialEllipticalGauss(fdir, fradius, fradius2, fpangle)
                elif ext_fsource['Model_Form'] == '2D Gaussian':
                    if fpangle == 0.:
                        spatial = gammalib.GModelSpatialRadialShell(fdir,fradius, fradius2)
                    else:
                        print('{} modeled by elliptical ring, which is not implemented, skip'.format(fsource['Source_Name']))
                        new = 0.
                elif ext_fsource['Model_Form'] == 'Map':
                    print('{} modeled by spatial template, which is not implemented, skip'.format(fsource['Source_Name']))
                    new = 0
                else:
                    print('{} modeled by model type {}, which is not implemented, skip'.format(fsource['Source_Name'],ext_fsource['Model_Form']))
                    new = 0
                if new == 1:
                    newext +=1
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
