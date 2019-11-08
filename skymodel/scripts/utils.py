import gammalib
import numpy as np
import pdb

def format_ax(ax):
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.grid()
    return

def flux_Crab(model,emin,emax):
    emin = gammalib.GEnergy(emin,'TeV')
    emax = gammalib.GEnergy(emax,'TeV')
    flux = model.spectral().flux(emin, emax)
    # convert to Crab units (Meyer model, flux > 1 TeV in ph cm-2 s-1)
    flux /= 2.0744340476909142e-11
    return flux

def dist_from_gammalib(models):
    #extracts flux, longitude and latitude distribution from gammalib model container
    fluxes = []
    lons = []
    lats = []
    names = []
    # energy boundaries for flux calculation

    for model in models:
        src_dir = get_model_dir(model)
        lons.append(src_dir.l_deg())
        lats.append(src_dir.b_deg())
        flux = flux_Crab(model,1.,1000.)
        fluxes.append(flux)
        names.append(model.name())

    return lons, lats, fluxes, names

def get_model_dir(model):
    """
    extracts sky direction for either analytical model or DiffuseMap
    :param model: ~gammalib.GModelSky
    :return: src_dir: ~gammalib.GSkyDir
    """
    if model.spatial().type() == 'DiffuseMap':
        # retrieve direction from map
        src_dir = model.spatial().region().centre()
    else:
        # retrieve direction from analytical model
        src_dir = model.spatial().dir()
    return src_dir

def get_model_radius(model):
    """
    extracts radius for extended source
    :param model: ~gammalib.GModelSky
    :return: radius: float
    """
    if model.spatial().type() == 'PointSource':
        radius = 0.
    elif model.spatial().type() == 'RadialGaussian':
        radius = 2 * model['Sigma'].value()
    elif model.spatial().type() == 'EllipticalGaussian':
        radius = 2 * model['MajorRadius'].value()
    elif model.spatial().type() == 'RadialShell':
        radius = model['Radius'].value() + model['Width'].value()
    elif model.spatial().type() == 'RadialDisk':
        radius = model['Radius'].value()
    elif model.spatial().type() == 'DiffuseMap':
        radius = model.spatial().region().radius()
    else:
        print('model {} has spatial model type {} which is not implemented'.format(
            model.name(), model.spatial().type()))
    return radius

def find_source_to_delete(d,lon,lat,flux):

    # calculate distance
    # use flat approx, only valid for small distances
    # but it does not matter here
    distx = np.abs(d['GLON'] - lon)
    distx[distx > 180] = 360 - distx[distx > 180]
    disty = d['GLAT'] - lat
    dist = np.sqrt(distx**2 + disty**2)

    # calculate flux_ratio log
    frlog = np.log10(d['flux']/flux)

    # figure of merit to decide what source to eliminate
    fom = np.sqrt((dist/180)**2 + frlog**2)

    # eliminate closer source = minimum fom
    s = np.where(fom == np.min(fom))
    name = d['name'][s][0]

    # create new source dictionary without source to be removed
    for key in d.keys():
        if key == 'name':
            pass
        else:
            d[key] = d[key][d['name'] != name]
    d['name'] = d['name'][d['name'] != name]

    # prints for checking algorithm works correctly
    print('name {}, fom {}, dist{}, distx {}, disty {}, frlog {}'.format(name, fom[s],dist[s], distx[s], disty[s], frlog[s]))

    return name, d



