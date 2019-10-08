# ==========================================================================
# this script is partially based on
# make_pointings.py
# distributed as part of ctools
# http://cta.irap.omp.eu/ctools/
# under the terms of the GNU General Public License
# <http://www.gnu.org/licenses/>
# ==========================================================================

import gammalib
import math

nrows = 2
step = 2.25  # deg
obs_time = 0.5  # h
caldb = 'prod3b-v2'
rad = 5.  # FoV cut
idstart = 0
t0 = 7671.0 * 86400.0  # s, corresponds to (1-1-2021)
emin = 0.013  # TeV
emax = 199.  # TeV

# characteristics of GPS KSP from Science with CTA
GPS = [{'name': 'inner region', 'lmin': -60., 'lmax': 60., 'duration': 780., 'array': 'South'},
       {'name': 'Vela, Carina', 'lmin': 240., 'lmax': 300., 'duration': 180.,
        'array': 'South'},
       {'name': 'other South', 'lmin': 210., 'lmax': 240., 'duration': 60., 'array': 'South'},
       {'name': 'Cygnus, Perseus', 'lmin': 60., 'lmax': 150., 'duration': 450.,
        'array': 'North'},
       {'name': 'anticentre', 'lmin': 150., 'lmax': 210., 'duration': 150., 'array': 'North'}]


def get_positions(xmin, xmax, nrows, step, fixed_xstep=False):
    positions = []

    # Determine ystep and ymin
    ystep = step * math.sqrt(0.75)
    ymin = - ystep * (nrows - 1) / 2

    y = ymin

    for row in range(nrows):
        # Choose between fixed xstep or equilateral trangles
        if fixed_xstep:
            xstep = fixed_xstep
        else:
            xstep = step
        # Determine xstep and number of pointings in row. The xstep increases
        # with the cosine of the latitude so that the distance of the step
        # is invariant of latitude.
        xstep = xstep / math.cos(gammalib.deg2rad * y)
        nx = int((xmax - xmin) / xstep + 1.5)

        # Set x offset. For every second row the x position is displace by
        # half a step size.
        if row % 2 == 0:
            x = xmin
        else:
            x = xmin + 0.5 * xstep

        # Append pointings
        for pnt in range(nx):
            positions.append({'x': x, 'y': y})
            x += xstep
            if x >= 360.0:
                x -= 360.0

        # Increment y
        y += ystep

    # Return positions
    return positions


# ============== #
# Set start time #
# ============== #
def set_tmin_for_next_pointing(tmin, duration):
    """
    Set start time for next pointing

    Adds duration and 2 min for slew, and takes care of day and night. It is
    assumed that the night lasts from [0,10] and the day lasts from [10,24]

    Parameters
    ----------
    tmin : float
        Start time of current pointing (sec)
    duration : float
        Duration of current pointing (sec)

    Returns
    -------
    tmin : float
        Start time of new pointing
    """
    # Update start time, add 2 min for slew
    tmin += duration + 120.0

    # If we have more than 10 hours in day then go to next day
    hours_in_day = tmin / 3600.0 % 24.0
    if hours_in_day > 10.0:
        days = int(tmin / 86400.0)
        days += 1
        tmin = float(days) * 86400.0

    # Dump result
    # print(tmin, hours_in_day)

    # Return start time
    return tmin


def set_patch(tmin, xmin, xmax, duration, array):
    obsdef = []

    # set as many observations as allowed by total duration time and observation length
    n_pnt = int(duration / obs_time) - 2  # - 2 to allow for margin in get_positions
    xstep = 2 * (xmax - xmin) / n_pnt
    pointings = get_positions(xmin, xmax, nrows, step, fixed_xstep=xstep)

    # Set observations
    total_duration = 0
    for pnt in pointings:
        # Set positions, start time and duration
        obs = {'lon': pnt['x'], 'lat': pnt['y'], 'rad': rad, \
               'tmin': tmin, 'duration': obs_time * 3600.}

        # Update start time for next pointing
        tmin = set_tmin_for_next_pointing(tmin, obs_time * 3600.)

        # Find Dec for pointing
        pnt_dir = gammalib.GSkyDir()
        pnt_dir.lb_deg(pnt['x'], pnt['y'])
        dec = pnt_dir.dec_deg()

        # Automatic IRF setting
        # Set geographic latitude of array
        if array == 'North':
            geolat = +28.7569
        else:
            geolat = -24.58
        # Compute best possible zenith angle
        zenith = abs(dec - geolat)
        # Assign IRF zenith
        if zenith < 30.0:
            irfz = 20
        elif zenith < 50.0:
            irfz = 40
        else:
            irfz = 60
        irf = '{}_z{}_50h'.format(array, irfz)

        # Set IRF information
        obs['caldb'] = caldb
        obs['irf'] = irf

        # Append observation
        obsdef.append(obs)

        # Add to total duration
        total_duration += obs_time

    # check that total duration is preserved
    if total_duration == duration:
        pass
    else:
        print(
            'WARNING: total duration of {} h and requested observation of {} do not match'.format(
                total_duration, duration))

    return obsdef


def write_obsdef(filename, obsdef, idstart, emin, emax):
    """
    Write observation definition file

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
    f.write('id,ra,dec,rad,tmin,duration,caldb,irf,emin,emax\n')

    # Initialise identifier
    obsid = idstart

    # Loop over pointings
    for obs in obsdef:
        # We have lon,lat then convert into RA,Dec
        lon = obs['lon']
        lat = obs['lat']
        pnt = gammalib.GSkyDir()
        pnt.lb_deg(lon, lat)
        ra = pnt.ra_deg()
        dec = pnt.dec_deg()

        # Write information
        f.write('%6.6d,%8.4f,%8.4f,%8.4f,%.4f,%.4f,%s,%s,%.3f,%.1f\n' %
                (obsid, ra, dec, obs['rad'], obs['tmin'], obs['duration'], obs['caldb'], \
                 obs['irf'], emin, emax))

        # Increment identifier
        obsid += 1

    # Close file
    f.close()


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    obsdef = []

    # Southern survey
    tmin = t0
    for patch in GPS:
        if patch['array'] == 'South':
            obsdef.extend(
                set_patch(tmin, patch['lmin'], patch['lmax'], patch['duration'],
                          patch['array']))
            # Update start time
            tmin = set_tmin_for_next_pointing(obsdef[-1]['tmin'],
                                              obsdef[-1]['duration'])

    # Northern survey
    tmin = t0
    for patch in GPS:
        if patch['array'] == 'North':
            obsdef.extend(
                set_patch(tmin, patch['lmin'], patch['lmax'], patch['duration'],
                          patch['array']))
            # Update start time
            tmin = set_tmin_for_next_pointing(obsdef[-1]['tmin'],
                                              obsdef[-1]['duration'])

    write_obsdef('GPS_pointings.dat', obsdef, idstart, emin, emax)
