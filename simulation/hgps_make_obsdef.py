#! /usr/bin/env python
# ==========================================================================
# Generate HGPS observation definition file
#
# Copyright (C) 2019 Juergen Knoedlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================
import os
import sys
import math
import gammalib
import matplotlib.pyplot as plt


# =============================================== #
# Set positions for dedicated survey observations #
# =============================================== #
def set_survey_pointings(lmin, lmax, b, dl=0.7):
    """
    """
    # Initialise positions
    pointings = []

    # Initialise longitude
    lon = lmin
    if lon > 180.0:
        lon -= 360.0

    # Loop over longitudes
    while lon < lmax:

        # Loop over latitudes
        for lat in b:

            # Append pointing
            pointings.append({'lon': lon, 'lat': lat})

        # Increment longitude
        lon += dl

    # Return pointings
    return pointings


# ============= #
# Set pointings #
# ============= #
def set_pointings():
    """
    """
    # Initialise pointings
    pointings = []

    # Set pointings
    pointings.extend(set_survey_pointings(244.5, 77.5, [-1.8, 1.0]))

    # Return pointings
    return pointings


# ==================================== #
# Setup H.E.S.S. Galactic plane survey #
# ==================================== #
def set_hgps(pointings, exposure=28.0*60.0, caldb='hess'):
    """
    Setup H.E.S.S. Galactic plane survey

    Parameters
    ----------
    pointings : list of dict
        HGPS pointings

    Returns
    -------
    obsdef : list of dict
        List of pointing definitions
    """
    # Initialise observation definition
    obsdef = []

    # Initialise start time in seconds (1-1-2021)
    tmin = 7671.0 * 86400.0

    # Set geographic latitude of H.E.S.S. array
    geolat = -23.27178

    # Loop over all pointings
    for pointing in pointings:

        # Compute Right Ascension and Declination of pointing
        pnt = gammalib.GSkyDir()
        pnt.lb_deg(pointing['lon'], pointing['lat'])
        ra  = pnt.ra_deg()
        dec = pnt.dec_deg()

        # Compute best possible zenith angle
        zenith = abs(dec - geolat)

        # Set IRF
        irf = 'dummy'

        # Set positions, start time and duration
        obs = {'lon': pointing['lon'], \
               'lat': pointing['lat'], \
               'tmin': tmin, \
               'duration': exposure,
               'zenith': zenith,
               'caldb': caldb,
               'irf': irf}

        # Update start time for next pointing; add 2 min for slew
        tmin += exposure + 120.0

        # Append observation
        obsdef.append(obs)

    # Return observation definition
    return obsdef


# ============================================== #
# Set map from observation definition dictionary #
# ============================================== #
def set_map(obsdef, radius=2.0):
    """
    Set map from observation definition dictionary

    Parameters
    ----------
    obsdef : list of dict
        List of pointing definitions
    """
    # Create sky map
    map = gammalib.GSkyMap('CAR', 'GAL', 340.0, 0.0, -0.1, 0.1, 950, 100)

    # Loop over observations
    for obs in obsdef:

        # Set sky region
        centre   = gammalib.GSkyDir()
        centre.lb_deg(obs['lon'], obs['lat'])
        circle   = gammalib.GSkyRegionCircle(centre, radius)
        region   = gammalib.GSkyRegionMap(circle)
        exposure = region.map().copy()
        exposure *= obs['duration']/3600.0

        # Add sky region
        map += exposure

    # Return map
    return map


# ====================================== #
# Plot observation definition dictionary #
# ====================================== #
def plot_obsdef(obsdef):
    """
    Plot observation definition dictionary

    Parameters
    ----------
    obsdef : list of dict
        List of pointing definitions
    """
    # Create figure
    fig = plt.figure('Map', (20,3))
    fig.subplots_adjust(left=0.05, right=0.98, top=0.98, bottom=0.02)
    sub = plt.subplot()

    # Create map from observation definition dictionary
    map = set_map(obsdef)

    # Create array from skymap
    array = []
    v_max = 0.0
    for iy in range(map.ny()):
        row = []
        for ix in range(map.nx()):
            index = ix+(map.ny()-iy-1)*map.nx()
            value = map[index]
            row.append(value)
        array.append(row)

    # Get skymap boundaries
    lon_min = map.pix2dir(gammalib.GSkyPixel(0.0,map.ny()/2)).l_deg()
    lon_max = map.pix2dir(gammalib.GSkyPixel(map.nx(),map.ny()/2)).l_deg()
    lat_min = map.pix2dir(gammalib.GSkyPixel(map.nx()/2,0)).b_deg()
    lat_max = map.pix2dir(gammalib.GSkyPixel(map.nx()/2,map.ny())).b_deg()
    aspect  = abs((lon_max-lon_min)/(lat_max-lat_min))

    # Show Aitoff projection
    c    = sub.imshow(array, extent=(lon_min,lon_max,lat_min,lat_max),
                      cmap=plt.get_cmap('jet'))
    cbar = plt.colorbar(c, orientation='horizontal', shrink=0.2)
    sub.set_xlabel('Longitude (deg)', fontsize=14)
    sub.set_ylabel('Latitude (deg)', fontsize=14)

    # Show map
    plt.show()

    # Return
    return


# ======================================= #
# Write observation definition dictionary #
# ======================================= #
def write_obsdef(filename, obsdef, idstart):
    """
    Write observation definition dictionary

    Parameters
    ----------
    filename : str
        Observation definition file name
    obsdef : list of dict
        List of pointing definitions
    idstart : int
        First identifier of observation definition file
    """
    # Open file
    f = open(filename, 'w')

    # Write header
    f.write('id,ra,dec,tmin,duration,caldb,irf,emin,emax\n')

    # Initialise identifier
    obsid = idstart

    # Loop over pointings
    for obs in obsdef:

        # If we have lon,lat then convert into RA,Dec
        if 'lon' in obs and 'lat' in obs:
            lon = obs['lon']
            lat = obs['lat']
            pnt = gammalib.GSkyDir()
            pnt.lb_deg(lon,lat)
            ra  = pnt.ra_deg()
            dec = pnt.dec_deg()
        else:
            ra  = obs['ra']
            dec = obs['dec']

        # Set site dependent energy thresholds
        if 'South' in obs['irf']:
            emin =   0.030
            emax = 160.0
        else:
            emin =  0.030
            emax = 50.0

        # Write information
        f.write('%6.6d,%8.4f,%8.4f,%.4f,%.4f,%s,%s,%.3f,%.1f\n' %
                (obsid, ra, dec, obs['tmin'], obs['duration'], obs['caldb'], \
                 obs['irf'], emin, emax))

        # Increment identifier
        obsid += 1

    # Close file
    f.close()

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Dump header
    print('********************************************')
    print('* Create HGPS observation definition files *')
    print('********************************************')

    # Create pointings
    pointings = set_pointings()

    # Create observation definition dictionary
    obsdef = set_hgps(pointings)

    # Plot observation definition dictionary
    plot_obsdef(obsdef)

    # Write observation definition ASCII file
    write_obsdef('hgps.dat', obsdef, 110000)

    # We are done
    print('... done.')
