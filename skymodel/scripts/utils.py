import gammalib
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

def dist_from_gammalib(models):
    #extracts flux, longitude and latitude distribution from gammalib model container
    fluxes = []
    lons = []
    lats = []
    # energy boundaries for flux calculation
    emin = gammalib.GEnergy(1.,'TeV')
    emax = gammalib.GEnergy(1000,'TeV')
    for model in models:
        src_dir = get_model_dir(model)
        lons.append(src_dir.l_deg())
        lats.append(src_dir.b_deg())
        flux = model.spectral().flux(emin,emax)
        # convert to Crab units (Meyer model, flux > 1 TeV in ph cm-2 s-1)
        flux /= 2.0744340476909142e-11
        fluxes.append(flux)

    return lons, lats, fluxes

def get_model_dir(model):
    """
    extracts sky direction for either analytical model or DiffuseMap
    :param model: ~gammalib.GModelSky
    :return: src_dir: ~gammalib.GSkyDir
    """
    if model.spatial().type() == 'DiffuseMap':
        # # retrieve source direction from map, use mid pixel as approximation
        dmap = model.spatial().map()
        src_dir = dmap.inx2dir(int(dmap.npix() / 2))
    else:
        # retrieve direction from analytical model
        src_dir = model.spatial().dir()
    return src_dir