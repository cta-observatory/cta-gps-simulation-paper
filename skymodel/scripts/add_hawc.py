import numpy as np
from utils import *
import gammalib

hawc_table2 = '../known-sources/external-input/HAWC_cat_table2_data.txt'
hawc_table3 = '../known-sources/external-input/HAWC_table3data.txt'

# load HAWC catalog table 2 and 3
hawc_sources = np.genfromtxt(hawc_table2,
                             dtype = [('name','<U15'),
                                      ('radius',np.double),
                                      ('ts',np.float),
                                      ('ra',np.double),
                                      ('dec',np.double),
                                      ('glon',np.double),
                                      ('glat',np.double),
                                      ('pos_err',np.float),
                                      ('dist_tevcat',np.float),
                                      ('assoc','<U15')]
                             )

hawc_spectra = np.genfromtxt(hawc_table3,
                             dtype=[('name', '<U15'),
                                    ('radius', np.double),
                                    ('index', np.double),
                                    ('index_error', np.double),
                                    ('F7', np.double),
                                    ('F7_error', np.double),
                                    ('assoc','<U15')]
                             )

# The HAWC Catalog gives spectra as power laws with reference energy at 7 TeV
eref = gammalib.GEnergy(7., 'TeV')

def append_hawc(models, bmax, pwn_models, pwn_dict, pwn_distx, pwn_disty, pwn_radr, pwn_frlog, dist_sigma=3., radmin = 0.05):
    """
    Append missing models from HAWC catalog to gammalib model container
    :param models: input gammalib model container
    :param models: ~gammalib.GModels, gammalib model container
    :param bmax: float, maximum latitude of models to include
    :param dist_sigma: float, minimum distance in sigma's from pre-existing source to add as new
    """

    # pre-filter known sources based on information in the hawc catalog
    sources = hawc_sources[hawc_sources['dist_tevcat'] > dist_sigma * hawc_sources['pos_err']]

    # prepare new gammalib model container
    newpt = 0
    newext = 0
    newmodels = gammalib.GModels()
    warning = ''
    n_pwn_del = 0

    # loop over sources and compare with input models
    for hsource in sources:
        # assume source is new
        new = 1
        # hawc source position as gammalib.GSkyDir
        hdir = gammalib.GSkyDir()
        ra = np.double(hsource['ra'])
        dec = np.double(hsource['dec'])
        hdir.radec_deg(ra, dec)
        # hawc source radius
        hradius = hsource['radius']
        ######################################### treat pointlike sources for HAWC ###########
        if hradius == 0.:
            # case of sources pointlike for HAWC
            # loop over gammalib container and determine closest neighbor
            # initialize at 1000
            dist_min = 1000.
            for source in models:
                dir = get_model_dir(source)
                dist = hdir.dist_deg(dir)
                dist_min = np.minimum(dist, dist_min)
            if dist_min < dist_sigma * hsource['pos_err']:
                # closeby source foud in container, source will not be added
                new = 0
            else:
                # source will be added to container, set spatial model
                spatial = gammalib.GModelSpatialPointSource(hdir)
                newpt += 1
        ######################################### treat extended sources for HAWC ############
        else:
            # "corrected" distance, i.e, centre distance - radius of source - error on position
            # initialize at 1000
            dist_corr = 1000.
            for source in models:
                dir = get_model_dir(source)
                dist = dir.dist_deg(dir)
                # get source radius according to model used
                radius = get_model_radius(source)
                # use as radius the max between HAWC and model
                radius = np.maximum(hradius,radius)
                # subtract from distance the max radius
                dist -= radius + hsource['pos_err']
                dist_corr = np.minimum(dist,dist_corr)
            if dist_corr < 0.:
                # overlapping source foud in container, source will not be added
                new = 0
            else:
                # source will be added to container, set spatial model
                spatial = gammalib.GModelSpatialRadialDisk(hdir,hradius)
        ######################################### spectra #####################################
        if new == 1 and np.abs(hdir.b_deg()) < bmax:
            # spectral model
            # identify spectral model in table 3
            spectrum = hawc_spectra[(hawc_spectra['name'] == hsource['name']) &
                                    (hawc_spectra['radius'] == hsource['radius'])][0]
            # flux scaling factor
            # 1e-15 ph TeV-1 -> ph MeV-1
            spectral = gammalib.GModelSpectralPlaw(1.e-21 * spectrum['F7'], spectrum['index'], eref)
            if spectrum['index'] > -2.4:
                warning += '2HWC source {} has hard index {}'.format(hawc_spectra['name'],spectrum['index'])
            # assemble model and append to container
            model = gammalib.GModelSky(spatial, spectral)
            model.name(hsource['name'])
            newmodels.append(model)
            # delete synthetic PWN
            print('PWN')
            rname, pwn_dict, distx, disty, radr, frlog = find_source_to_delete(pwn_dict,hdir.l_deg(),hdir.b_deg(),
                                                                               get_model_radius(model),flux_Crab(model,1.,1000.),radmin=radmin)
            pwn_models.remove(rname)
            n_pwn_del += 1
            pwn_distx.append(distx)
            pwn_disty.append(disty)
            pwn_radr.append(radr)
            pwn_frlog.append(frlog)
            print('-----')
        else:
            # source already present, skip
            pass

    # add new models to original container
    for model in newmodels:
        models.append(model)

    return models, newpt, newext, warning, pwn_models, n_pwn_del, pwn_distx, pwn_disty, pwn_radr, pwn_frlog

