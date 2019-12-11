import gammalib
import numpy as np
from astropy.io import fits
import os
from utils import *

# convert binary population model from txt format provided by Guillaume
# to gammalib format

# some additional binary population parameters not in the files (see README in binpop)
# spectral index
gamma = 2.5
# maximum latitude
bmax = 2.
# threshold energy for lightcurve (MeV)
eref = 1.e6
# reference time for phasogram
t0 = gammalib.GTime()
t0.mjd(51544.5)


def phasogram_file(phases, fname, outdir):
    # derive normalization from flux
    norm = phases[:, 1] / np.max(phases[:, 1])

    # create fits columns and put them into table
    phase = fits.Column(name='PHASE', format='D', array=phases[:, 0])
    norm = fits.Column(name='NORM', format='D', array=norm)
    tbhdu = fits.BinTableHDU.from_columns([phase, norm])
    tbhdu.verify('fix')

    # write file
    tbhdu.writeto(fname, overwrite = True)

    return


def bin2gammalib(filename, flux_thresh, outdir):
    data = open(filename).readlines()

    # convert phasogram to numpy array
    phasogram = np.genfromtxt(data[14:])

    # only include model if peak in lightcurve is > threshold
    if np.max(phasogram[:, 2]) > flux_thresh:

        # final name
        num = int(data[0].split('#')[-1])
        name = 'bin{:03d}'.format(num)

        # spatial model
        # read Galactic longitude
        lon = float(data[1].split('   ')[-1])
        # Galactic latitude is random with Gaussian distribution with 95% containment at bmax
        lat = np.random.normal(0, bmax / 2)
        src_dir = gammalib.GSkyDir()
        src_dir.lb_deg(lon, lat)
        spatial = gammalib.GModelSpatialPointSource(src_dir)

        # spectral model
        # average flux in ph/cm2/s
        f0 = np.average(phasogram[:, 1])
        # convert to differential flux at eref in ph/cm2/s/MeV
        n0 = (gamma - 1) * f0 / eref
        spectral = gammalib.GModelSpectralPlaw(n0, -gamma, gammalib.GEnergy(np.double(eref), 'MeV'))

        # orbital phase model
        # create phasogram file
        fname = 'phasecurve_' + name + '.fits'
        phasogram_file(phasogram, outdir + fname, outdir)
        # period in days
        per = float(data[4].split(' ')[-1])
        # convert to frequency in s
        per *= 86400
        f0 = 1 / per
        temporal = gammalib.GModelTemporalPhaseCurve(gammalib.GFilename(fname),
                                                     t0, np.random.random(), # phase at t0 randomly drawn
                                                     f0, 0., 0., # orbit is stable, no derivatives
                                                     1., True) # normalized so that spectral model represents average flux

        # assemble final model
        model = gammalib.GModelSky(spatial, spectral, temporal)
        model.name(name)

    else:
        # empty model
        model = gammalib.GModelSky()

    return model

def get_binpop_models(datadir,flux_thresh,outdir):

    models = gammalib.GModels()

    # get file names
    files = []
    # r=root, d=directories, f = files
    for r, d, f in os.walk(datadir):
        for file in f:
            if '.txt' in file:
                files.append(os.path.join(r, file))

    # load models
    for file in files:
        model = bin2gammalib(file,flux_thresh,outdir)
        if model.size() > 0:
            models.append(model)
        else:
            pass

    # dictionary with parameters
    lons, lats, radii, fluxes, names = dist_from_gammalib(models)

    d = {'name' : np.array(names),
            'GLON' : np.array(lons),
            'radius' : np.array(radii),
            'GLAT' : np.array(lats),
            'flux' : np.array(fluxes)}

    return models, d
