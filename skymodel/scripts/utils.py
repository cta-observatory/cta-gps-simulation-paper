import gammalib
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
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

def cube_flux(model,emin,emax):
    # get data
    filename = model.spatial().filename().file()
    hdulist = fits.open(filename)
    cube = hdulist[0].data
    energies = hdulist[1].data['Energy']
    binsize1 = np.abs(hdulist[0].header['CDELT1'])
    binsize2 = np.abs(hdulist[0].header['CDELT2'])
    # multiply by solid angle (only works for tangential projection or for reference latitude ~0)
    cube *= np.deg2rad(binsize1) * np.deg2rad(binsize2)
    fluxes = np.sum(cube, axis=(1, 2))
    # correct by spectral model
    for s, energy in enumerate(energies):
        fluxes[s] *= model.spectral().eval(gammalib.GEnergy(np.double(energy),'MeV'))
    # select energy range
    # emin, emax are in TeV and the energy vector in MeV
    # just select on bins, may be inaccurate close to threshold
    fluxes = fluxes[(energies >= 1.e6 * emin) & (energies <= 1.e6 * emax)]
    energies = energies[(energies >= 1.e6 * emin) & (energies <= 1.e6 * emax)]
    # integrate over energy in piece-wise power law approximation
    gamma = - (np.log(fluxes[1:]) - np.log(fluxes[:-1])) / (np.log(energies[1:]) - np.log(energies[:-1]))
    # integral flux in individual bins between two nodes
    int = fluxes[:-1] * energies[:-1] / (-gamma + 1) * (np.power(energies[1:]/energies[:-1],-gamma+1) - 1)
    # sum over bins
    int = np.sum(int)
    return int

def flux_Crab(model,Emin,Emax):
    emin = gammalib.GEnergy(Emin,'TeV')
    emax = gammalib.GEnergy(Emax,'TeV')
    # deal with special case of cubes
    if model.spatial().type() == 'DiffuseMapCube':
        # get cube flux
        flux = cube_flux(model,Emin,Emax)
    else:
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
        if model.classname() == 'GModelSky':
            src_dir = get_model_dir(model)
            lons.append(src_dir.l_deg())
            lats.append(src_dir.b_deg())
            rad = get_model_radius(model)
            radii.append(rad)
            flux = flux_Crab(model,emin,emax)
            fluxes.append(flux)
            names.append(model.name())
        else:
            pass

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
    elif model.spatial().type() == 'DiffuseMapCube':
        # region does not work, extract manually from the map
        # call the energies method to load the cube
        model.spatial().energies()
        # assume it is center of the map
        sl = model.spatial().cube().extract(0)
        ctr_pix = gammalib.GSkyPixel((sl.nx() - 1) / 2, (sl.ny() - 1) / 2)
        src_dir = sl.pix2dir(ctr_pix)
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
    try:
        if model.spatial().type() == 'DiffuseMapCube':
            # region not implemented in gammalib, extract manually extent of the map
            # call the energies method to load the cube
            model.spatial().energies()
            # half size along x direction
            slice = model.spatial().cube().extract(0)
            ctr_pix = gammalib.GSkyPixel((slice.nx() - 1) / 2, (slice.ny() - 1) / 2)
            ctr_dir = slice.pix2dir(ctr_pix)
            bor_pix = gammalib.GSkyPixel(0., (slice.ny() - 1) / 2)
            bor_dir = slice.pix2dir(bor_pix)
            radius = bor_dir.dist_deg(ctr_dir)
        else:
            circle = gammalib.GSkyRegionCircle(model.spatial().region())
            radius = circle.radius()
    except:
        print('Cannot extract radius for model {} with spatial model type {}'.format(
            model.name(), model.spatial().type()))
    return radius

def delete_source_fom(distx,disty,radr, frlog):
    fom = np.sqrt((distx / 180) ** 2 + (disty / 10) ** 2 + (radr / 1.) ** 2 + (frlog / 0.3) ** 2)
    return fom

def pop_source(d,name):
    # remove source from dictionary by name
    # boolean mask to identify sources to keep in final dictionary
    m = (d['name'] != name)
    # handle composites for which name is an array
    if len(np.shape(m)) > 1:
        m = np.product(m, axis=1) # both names must coincide to be the same composite
        m = m.astype('bool')
    for key in d.keys():
        if key == 'name':
            pass
        else:
            d[key] = d[key][m]
    d['name'] = d['name'][m]
    return d

def find_source_to_delete(d,lon,lat,rad,flux, radmin =0.05):


    # calculate distance in lon and lat
    distx = np.abs(d['GLON'] - lon)
    distx[distx > 180] = 360 - distx[distx > 180]
    disty = d['GLAT'] - lat
    #dist = np.sqrt(distx**2 + disty**2)

    # calculate radius relative difference
    # take into the fact that current instruments cannot resolve objects with radii < 0.05 deg
    # based on minimum measured size in HGPS
    if rad < radmin:
        rad = radmin
    radr = 0.5 * (np.maximum(d['radius'],radmin) - rad) / (np.maximum(d['radius'],radmin) + rad)

    # calculate flux_ratio log
    frlog = np.log10(d['flux']/flux)

    # figure of merit to decide which source to eliminate
    fom = delete_source_fom(distx,disty,radr,frlog)

    # eliminate closer source = minimum fom
    s = np.where(fom == np.min(fom))
    name = d['name'][s][0]
    d = pop_source(d,name)

    # # prints for checking algorithm works correctly
    # print('name {}, fom {}, distx {}, disty {}, radr {}, frlog {}'.format(name, fom[s], distx[s], disty[s], radr[s], frlog[s]))

    return name, d, distx[s], disty[s], radr[s], frlog[s]

def get_syn_model(filename,fmin,bmax,emin=1.,emax=1000.):
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
        if fluxes[s] > fmin and np.abs(lats[s]) < bmax:
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

def plot_del_sources(distx,disty,radr, frlog,namestr,namefull):
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
    fom = delete_source_fom(xv, yv, 0.2, 0.3)
    cs = ax4.contour(xv, yv, fom)
    ax4.clabel(cs, inline=1, fontsize=10)
    fig4.savefig('{}XY.png'.format(namestr), dpi=300)

    fig5 = plt.figure('Deleted {} Flux-radius'.format(namefull))
    ax5 = plt.subplot()
    ax5.set_xlabel("Relative radius difference", fontsize=14)
    ax5.set_ylabel('Log10(flux ratio)', fontsize=14)
    format_ax(ax5)
    ax5.scatter(radr, frlog)
    # add FOM contours
    xmin, xmax = ax5.get_xlim()
    ymin, ymax = ax5.get_ylim()
    xv, yv = np.meshgrid(np.linspace(xmin, xmax, 50),
                         np.linspace(ymin, ymax, 50))
    fom = delete_source_fom(20, 2, xv, yv)
    cs = ax5.contour(xv, yv, fom)
    ax5.clabel(cs, inline=1, fontsize=10)
    fig5.savefig('{}Flux-rad.png'.format(namestr), dpi=300)

    return

def set_composites(pwn_dict,snr_dict):

    # create arrays to store output quantities
    names = []
    lons = np.array([])
    lats = np.array([])
    radii = np.array([])
    fluxes = np.array([])

    # loop over SNRs
    for s, snrname in enumerate(snr_dict['name']):
        #calculate distance along x and y for all pwne
        distx = snr_dict['GLON'][s] - pwn_dict['GLON']
        disty = snr_dict['GLAT'][s] - pwn_dict['GLAT']
        # total distance (flat approx)
        dist = np.sqrt(distx**2 + disty**2)
        # calculate radius sum
        sumrad = snr_dict['radius'][s] + pwn_dict['radius']
        # if there is an overlapping source
        if np.any(dist < sumrad):
            distance = np.min(dist)
            pwnname = pwn_dict['name'][dist==distance][0]
            names.append([pwnname,snrname])
            # for GLON and GLAT take flux-weighted average of individual objects
            pwn_lon = pwn_dict['GLON'][dist==distance]
            pwn_lat = pwn_dict['GLAT'][dist==distance]
            pwn_flux = pwn_dict['flux'][dist == distance]
            lon = snr_dict['GLON'][s] * snr_dict['flux'][s] + pwn_lon * pwn_flux
            lon /= (snr_dict['flux'][s] + pwn_flux)
            lons = np.append(lons,lon)
            lat = snr_dict['GLAT'][s] * snr_dict['flux'][s] + pwn_lat * pwn_flux
            lat /= (snr_dict['flux'][s] + pwn_flux)
            lats = np.append(lats,lat)
            # radius is max encompassing two objects
            rad_pwn = pwn_dict['radius'][dist==distance][0]
            if distance < np.abs(snr_dict['radius'][s] - rad_pwn):
                # one of the objects is contained
                rad = np.maximum(snr_dict['radius'][s], rad_pwn)
            else:
                # partially overlapping
                rad = (distance + rad_pwn + snr_dict['radius'][s]) / 2
            radii = np.append(radii,rad)
            # flux is sum of two fluxes
            flux = pwn_dict['flux'][dist==distance] + snr_dict['flux'][s]
            fluxes = np.append(fluxes,flux)
            # # test print
            msg = '{}/{} distance: {} deg, radii: {}/{} deg'.format(pwnname,snrname, distance,rad_pwn,snr_dict['radius'][s])
            print(msg)

    # convert name to numpy array
    names = np.array(names)

    # check that the same pwn is not used twice
    m = np.zeros(len(names), dtype=bool)
    m[np.unique(names[:,0], return_index=True)[1]] = True
    if len(names[:,0][~m]) > 0:
        print('WARNING: the same PWN appears as part of two composites')
        print(names[:,0][~m])
    else:
        pass

    # create composite dictionary with source parameters
    d = {'name'   : np.array(names),
         'GLON'   : np.array(lons),
         'GLAT'   : np.array(lats),
         'radius' : np.array(radii),
         'flux'   : np.array(fluxes)}

    # remove composite members for snr and pwn dictionaries
    for name in d['name']:
        pwn_dict = pop_source(pwn_dict,name[0])
        snr_dict = pop_source(snr_dict,name[1])

    # return results
    return d, pwn_dict, snr_dict









