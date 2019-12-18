#! /usr/bin/env python
# ==========================================================================
# Display pointings
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
import sys
import math
import gammalib
import cscripts
import matplotlib.pyplot as plt
from   matplotlib.colors import LogNorm
from matplotlib.image import NonUniformImage


# =============================== #
# Extract pointings from CSV file #
# =============================== #
def get_pointings(filename):
    """
    Extract pointings from CSV file

    Parameters
    ----------
    filename : str
        File name of observation definition CSV file

    Returns
    -------
    pnt : list of dict
        Pointings
    """
    # Initialise pointings
    pnt = []

    # Open CSV file
    csv = gammalib.GCsv(filename, ',')

    # Loop over rows
    for i in range(csv.nrows()):

        # Skip header
        if i == 0:
            continue

        # Extract information
        id       = csv[i,1]
        lon      = float(csv[i,2])
        lat      = float(csv[i,3])
        tstart   = float(csv[i,5])
        duration = float(csv[i,6])
        irf      = csv[i,8]
        zenith   = float(csv[i,11])

        # Convert direction
        dir = gammalib.GSkyDir()
        dir.lb_deg(lon, lat)

        # Derive attributes
        south = 'South' in irf

        # Create entry
        entry = {'id': id, 'l': lon, 'b': lat,
                 'ra': dir.ra_deg(), 'dec': dir.dec_deg(),
                 'time': tstart, 'duration': duration, \
                 'zenith': zenith, 'south': south}

        # Append pointing
        pnt.append(entry)

    # Return pointings
    return pnt


# ===================== #
# Get binary population #
# ===================== #
def get_binpop(pathname='github/cta-gps-simulation-paper/skymodel/binpop'):
    """
    Get binary population

    Parameters
    ----------
    pathname : str
        Path name of binary population

    Returns
    -------
    pop : list of dict
        Binary population
    """
    # Initialise population
    pop = []

    # Loop over all files
    for i in range(200):

        # Build filename
        filename = '%s/grbinary%d.txt' % (pathname, i+1)

        # Read file
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()

        # Extract longitude (deg) and peridicity (days)
        lon = float(lines[1].split()[3])
        per = float(lines[4].split()[3])

        # Create entry
        entry = {'lon': lon, 'period': per}

        # Append binary
        pop.append(entry)

    # Return population
    return pop


# ============== #
# Get longitudes #
# ============== #
def get_longitudes(pnt, lon, rad=3.0):
    """
    Selects longitudes from a list of pointings

    Parameters
    ----------
    pnt : list of dict
        Pointings
    lon : float
        Galactic longitude (deg)
    rad : float, optional
        Selection radius (deg)

    Returns
    -------
    pnt : list of dict
        Pointings
    """
    # Initialise list of pointings
    pointings = []

    # Loop over all pointings
    for p in pnt:

        # Compute longitude difference
        dl = p['l'] - lon
        if dl > 180.0:
            dl -= 360.0
        if dl < -180.0:
            dl += 360.0

        # Append pointing if close to requested longitude
        if abs(dl) < rad:
            pointings.append(p)

    # Return pointings
    return pointings


# ================= #
# Get periodicities #
# ================= #
def get_periodicities(pnt):
    """
    Get periodicities for pointings

    Parameters
    ----------
    pnt : list of dict
        Pointings

    Returns
    -------
    periodicities : list of float
        Periodicities
    """
    # Initialise periodicities
    periodicities = []

    # Loop over pointings
    for i1, p1 in enumerate(pnt):

        # Loop over other pointings
        for i2, p2 in enumerate(pnt):

            # Skip all p2 that are not after p1
            if i2 <= i1:
                continue

            # Compute time differences in days
            t = abs(p2['time'] - p1['time']) / 86400.0

            # Skip zero times
            if t == 0:
                print('WARNING: pointings "%s" and "%s" at (l,b)=(%.2f,%.2f) and '
                      '(l,b)=(%.2f,%.2f) have identical times' %
                      (p1['id'],p2['id'],p1['l'],p1['b'],p2['l'],p2['b']))
                continue

            # Set entry
            entry = {'periodicity': t, 'duration': p1['duration']+p2['duration']}

            # Append entry
            periodicities.append(entry)

    # Return periodicities
    return periodicities


# ============== #
# Plot pointings #
# ============== #
def plot_pointings(pnt, rad=3.0):
    """
    Plot pointings

    Parameters
    ----------
    pnt : list of dict
        Pointings
    rad : float, optional
        Selection radius (deg)
    """
    # Create figure
    plt.figure(figsize=(10,6))

    # Setup figure
    ax = plt.gca()
    ax.cla()
    ax.set_xlim((180, -180))

    # Initialise ymin and ymax
    ymin = 0.0
    ymax = 0.0

    # Loop over pointings
    for p in pnt:

        # Get longitude
        l = p['l']
        if l > 180.0:
            l = l - 360.0

        # Set color
        if p['south']:
            color='r'
        else:
            color='b'

        # Get y axis
        b = p['b'] + (p['time'] - 662774400.0) / 86400.0

        # Set y limits
        if b < ymin:
            ymin = b
        if b > ymax:
            ymax = b

        # Flag pointings with identical time in green
        for p2 in pnt:
            if p['id'] != p2['id'] and p['south'] == p2['south'] and \
               p['time'] == p2['time']:
                color='c'
                print('WARNING: pointings "%s" and "%s" at (l,b)=(%.2f,%.2f) and '
                      '(l,b)=(%.2f,%.2f) have identical times' %
                      (p['id'],p2['id'],p['l'],p['b'],p2['l'],p2['b']))
                break

        # Set circle
        circle = plt.Circle((l, b), rad, color=color, fill=False)

        # Add circle
        ax.add_artist(circle)

    # Set y limits
    ax.set_ylim((ymin-10, ymax+10))

    # Plot title and labels
    plt.xlabel('Galactic longitude (deg)', fontsize=14)
    plt.ylabel('Galactic latitude (deg) + Time (days)', fontsize=14)
    plt.title('GPS paper pointings')

    # Return
    return


# ================== #
# Plot periodicities #
# ================== #
def plot_periodicities(pnt, rad=3.0):
    """
    Plot periodicities

    Parameters
    ----------
    pnt : list of dict
        Pointings
    rad : float, optional
        Selection radius (deg)
    """
    # Set parameters of known binaries
    knowns = [{'name': 'LS 5039',           'period':    3.9, 'lon':  16.88},
              {'name': 'LSI+61 303',        'period':   26.5, 'lon': 135.68},
              {'name': 'PSR B1259-63',      'period': 1236.7, 'lon': 304.18},
              {'name': 'HESS J1832-093',    'period':  60.0,  'lon':  22.48}, # from 1DC model
              {'name': 'HESS J0632+057',    'period': 315.0,  'lon': 205.67},
              {'name': '1FGL J1018.6-5856', 'period': 16.5,   'lon': 284.36}]
    
    # Create figure
    plt.figure(figsize=(10,6))

    # Initialise array
    P_min    =    0.02           # 30 min in days
    P_max    = 3650.0            # 10 years in days
    n_logP   = 50
    n_lons   = 360
    logP_min = math.log10(P_min)
    logP_max = math.log10(P_max)
    dlogP    = (logP_max-logP_min) / float(n_logP)
    array    = []
    lons     = [l-180.0 for l in range(n_lons)]
    logPs    = [i*dlogP+logP_min for i in range(n_logP)]
    Ps       = [math.pow(10.0,logP) for logP in (logPs)]

    # Get binary population
    pop = get_binpop()

    # Loop over all longitudes
    for lon in lons:

        # Get pointings for longitude
        pnt_lon = get_longitudes(pnt, lon, rad=rad)

        # Get periodicities for longitude
        periodicities = get_periodicities(pnt_lon)

        # Initialise row
        row = [0.0 for i in range(n_logP)]

        # Fill periodicities
        for periodicity in periodicities:

            # Compute logP
            logP   = math.log10(periodicity['periodicity'])
            i_logP = int((logP - logP_min) / dlogP)

            # Fill array
            if i_logP >= 0 and i_logP < n_logP:
                row[i_logP] += periodicity['duration']

        # Append row
        array.append(row)

    # Rotate array
    vmin         = 1.0e30
    vmax         = 0.0
    array_rotate = []
    for iy in range(n_logP):
        row = []
        #for ix in range(n_lons-1,-1,-1):
        for ix in range(n_lons):
            value = array[ix][iy]
            if value > vmax:
                vmax = value
            if value > 0.0 and value < vmin:
                vmin = value
            row.append(value)
        array_rotate.append(row)

    # Show color plot
    c    = plt.pcolor(lons, Ps, array_rotate, norm=LogNorm(vmin=vmin, vmax=vmax),
                      cmap='jet')
    cbar = plt.colorbar(c, orientation='vertical', shrink=0.7)
    cbar.ax.set_ylabel('exposure (sec)')
    plt.yscale('log')
    plt.gca().invert_xaxis()
    plt.xlabel('Galactic longitude (deg)', fontsize=14)
    plt.ylabel('Periodicity (days)', fontsize=14)
    plt.title('GPS paper periodicities')

    # Overplot binary population
    for binary in pop:
        lon = binary['lon']
        if lon > 180.0:
            lon -= 360.0
        plt.plot([lon],[binary['period']],'g.', markersize=8)

    # Overplot known binaries
    for known in knowns:
        lon = known['lon']
        if lon > 180.0:
            lon -= 360.0
        period = known['period']*1.4
        #plt.plot([lon],[known['period']],'r*', markersize=14)
        plt.plot([lon+0.4],[known['period']*1.04],'w*', markersize=14)
        plt.plot([lon+0.4],[known['period']*0.96],'w*', markersize=14)
        plt.plot([lon-0.4],[known['period']*1.04],'w*', markersize=14)
        plt.plot([lon-0.4],[known['period']*0.96],'w*', markersize=14)
        plt.plot([lon],[known['period']],'r*', markersize=14)
        plt.text(lon+0.4,period*1.04,known['name'], horizontalalignment='center', color='white')
        plt.text(lon+0.4,period*0.96,known['name'], horizontalalignment='center', color='white')
        plt.text(lon-0.4,period*1.04,known['name'], horizontalalignment='center', color='white')
        plt.text(lon-0.4,period*0.96,known['name'], horizontalalignment='center', color='white')
        plt.text(lon,period,known['name'], horizontalalignment='center')

    # Return
    return


# ============ #
# Plot zeniths #
# ============ #
def plot_zeniths(pnt, rad=3.0):
    """
    Plot zenith angles

    Parameters
    ----------
    pnt : list of dict
        Pointings
    rad : float, optional
        Selection radius (deg)
    """
    # Create figure
    plt.figure(figsize=(10,6))

    # Initialise array
    zenith_min =  0.0
    zenith_max = 50.0
    n_zeniths  = 20
    n_lons     = 360
    dzenith    = (zenith_max-zenith_min) / float(n_zeniths)
    array     = []
    lons      = [l-180.0 for l in range(n_lons)]
    zeniths   = [i*dzenith+zenith_min for i in range(n_zeniths)]

    # Initialise best zenith angle array
    best_zeniths = []

    # Fill best zenith angles
    for p in pnt:

        # Get longitude
        l = p['l']
        if l > 180.0:
            l = l - 360.0

        # Compute best zenith angle
        if p['south']:
            zenith = abs(p['dec'] + 24.58)
        else:
            zenith = abs(p['dec'] - 28.7569)

        # Create entry
        entry = {'lon': l, 'zenith': zenith}

        # Append best zenith angle
        best_zeniths.append(entry)

    # Sort list
    best_zeniths = sorted(best_zeniths, key=lambda k: k['lon'])

    # Extract arrays for plotting
    best_l = [k['lon']    for k in best_zeniths]
    best_z = [k['zenith'] for k in best_zeniths]

    # Loop over all longitudes
    for lon in lons:

        # Get pointings for longitude
        pnt_lon = get_longitudes(pnt, lon, rad=rad)

        # Initialise row
        row = [0.0 for i in range(n_zeniths)]

        # Fill zenith angles
        for p in pnt_lon:
            iz = int((p['zenith']-zenith_min)/dzenith)
            if iz >= 0 and iz < n_zeniths:
                row[iz] += p['duration']

        # Append row
        array.append(row)

    # Rotate array
    vmin         = 1.0e30
    vmax         = 0.0
    array_rotate = []
    for iy in range(n_zeniths-1,-1,-1):
        row = []
        for ix in range(n_lons-1,-1,-1):
            value = array[ix][iy]
            if value > vmax:
                vmax = value
            if value > 0.0 and value < vmin:
                vmin = value
            row.append(value)
        array_rotate.append(row)

    # Show color plot
    aspect = float((lons[n_lons-1]-lons[0])/(zenith_max-zenith_min))
    c      = plt.imshow(array_rotate, extent=(lons[n_lons-1],lons[0],zenith_min,zenith_max),
                        cmap='jet',
                        interpolation='nearest',
                        norm=LogNorm(vmin=vmin, vmax=vmax),
                        aspect=aspect)
    cbar   = plt.colorbar(c, orientation='vertical', shrink=0.7)
    cbar.ax.set_ylabel('exposure (sec)')
    plt.plot(best_l, best_z, 'k-')
    plt.xlabel('Galactic longitude (deg)', fontsize=14)
    plt.ylabel('Zenith angle (deg)', fontsize=14)
    plt.title('GPS paper zenith angles')

    # Return
    return


# ============== #
# Show pointings #
# ============== #
def show_pointings():
    """
    Show pointings
    """
    # Set filename
    filename = 'GPS_pointings.dat'

    # Get pointings
    pnt = get_pointings(filename)

    # Plot pointings
    plot_pointings(pnt)

    # Plot periodicities
    plot_periodicities(pnt)

    # Plot zeniths
    plot_zeniths(pnt)

    # Show plot
    plt.show()

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show pointings
    show_pointings()
