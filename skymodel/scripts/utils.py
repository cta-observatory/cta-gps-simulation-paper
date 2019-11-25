import gammalib
import numpy as np
import pdb
import matplotlib.pyplot as plt

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
    # convert to Crab units
    # Set Crab TeV spectral model based on a power law
    crab = gammalib.GModelSpectralPlaw(5.7e-16, -2.48, gammalib.GEnergy(0.3, 'TeV'))
    # calculate crab flux over the same energy range
    crab_flux = crab.flux(emin, emax)
    # remormalize crab flux so that it matches the Meyer model > 1 TeV (consistent with gamma-cat)
    crab_flux_1TeV = crab.flux(gammalib.GEnergy(1.,'TeV'), gammalib.GEnergy(1000.,'TeV'))
    # (Meyer model, flux > 1 TeV in ph cm-2 s-1)
    crab_flux *= 2.0744340476909142e-11 / crab_flux_1TeV
    flux /= crab_flux
    return flux

def dist_from_gammalib(models,emin=1.,emax=1000):
    #extracts flux, longitude and latitude distribution from gammalib model container
    fluxes = []
    lons = []
    lats = []
    names = []
    radii = []

    for model in models:
        src_dir = get_model_dir(model)
        lons.append(src_dir.l_deg())
        lats.append(src_dir.b_deg())
        rad = get_model_radius(model)
        radii.append(rad)
        flux = flux_Crab(model,emin,emax)
        fluxes.append(flux)
        names.append(model.name())

    return lons, lats, radii, fluxes, names

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

def delete_source_fom(distx,disty,frlog):
    fom = np.sqrt((distx / 180) ** 2 + (disty / 10) ** 2 + (frlog / 0.3) ** 2)
    return fom

def find_source_to_delete(d,lon,lat,flux):

    # calculate distance
    # use flat approx, only valid for small distances
    # but it does not matter here
    distx = np.abs(d['GLON'] - lon)
    distx[distx > 180] = 360 - distx[distx > 180]
    disty = d['GLAT'] - lat
    #dist = np.sqrt(distx**2 + disty**2)

    # calculate flux_ratio log
    frlog = np.log10(d['flux']/flux)

    # figure of merit to decide which source to eliminate
    fom = delete_source_fom(distx,disty,frlog)

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

    # # prints for checking algorithm works correctly
    # print('name {}, fom {}, distx {}, disty {}, frlog {}'.format(name, fom[s], distx[s], disty[s], frlog[s]))

    return name, d, distx[s], disty[s], frlog[s]

def get_syn_model(filename,fmin,emin=1.,emax=1000.):
    """
    Load synthetic population already given in gammalib format
    :param filename: str, name of XML file
    :param fmin: float, minimum flux (Crab)
    :param emin: float, minimum energy for flux computation (TeV)
    :param emax: float, maximum energy for flux computation (TeV)
    :return: ~gammalib.GModels, dic, output models and dictionary of their properties
    """

    # load input models models
    models = gammalib.GModels(filename)

    # parameter distribution for the requested energy range
    lons, lats, radii, fluxes, names = dist_from_gammalib(models,emin=emin,emax=emax)

    # output models
    outmodels = gammalib.GModels()
    for s, model in enumerate(models):
        if fluxes[s] > fmin:
            outmodels.append(model)
        else:
            pass

    # regenerate distributions for standard energy range
    lons, lats, radii, fluxes, names = dist_from_gammalib(outmodels,emin=1.,emax=1000.)

    # create dictionary with source parameters
    d = {'name'   : np.array(names),
         'GLON'   : np.array(lons),
         'GLAT'   : np.array(lats),
         'radius' : np.array(radii),
         'flux'   : np.array(fluxes)}

    return outmodels, d

def plot_del_sources(distx,disty,frlog,namestr,namefull):
    fig4 = plt.figure('Deleted {} XY'.format(namefull))
    ax4 = plt.subplot()
    ax4.set_xlabel("Longitude distance (deg)", fontsize=14)
    ax4.set_ylabel('Latitude distance (deg)', fontsize=14)
    format_ax(ax4)
    ax4.scatter(distx, disty)
    # add FOM contours
    xmin, xmax = ax4.get_xlim()
    ymin, ymax = ax4.get_ylim()
    xv, yv = np.meshgrid(np.linspace(xmin, xmax, 50),
                         np.linspace(ymin, ymax, 50))
    fom = delete_source_fom(xv, yv, 0.3)
    cs = ax4.contour(xv, yv, fom)
    ax4.clabel(cs, inline=1, fontsize=10)
    fig4.savefig('{}XY.png'.format(namestr), dpi=300)

    fig5 = plt.figure('Deleted {} Flux-X'.format(namefull))
    ax5 = plt.subplot()
    ax5.set_xlabel("Longitude distance (deg)", fontsize=14)
    ax5.set_ylabel('Log10(flux ratio)', fontsize=14)
    format_ax(ax5)
    ax5.scatter(distx, frlog)
    # add FOM contours
    xmin, xmax = ax5.get_xlim()
    ymin, ymax = ax5.get_ylim()
    xv, yv = np.meshgrid(np.linspace(xmin, xmax, 50),
                         np.linspace(ymin, ymax, 50))
    fom = delete_source_fom(xv, 2, yv)
    cs = ax5.contour(xv, yv, fom)
    ax5.clabel(cs, inline=1, fontsize=10)
    fig5.savefig('{}Flux-X.png'.format(namestr), dpi=300)

    return





