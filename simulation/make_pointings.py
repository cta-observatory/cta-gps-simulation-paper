#! /usr/bin/env python
# ==========================================================================
# Generates schedule pointing sequence for GPS KSP
#
# This script is partially based on
#
# make_pointings.py
#
# distributed as part of ctools
# http://cta.irap.omp.eu/ctools/
# under the terms of the GNU General Public License
# <http://www.gnu.org/licenses/>
# ==========================================================================
import copy
import math
import gammalib


# ========== #
# Parameters #
# ========== #
nrows    = 2                 # Number of rows
step     = 2.25              # Stepsize (deg)
obs_time = 0.5               # Pointing duration (h)
caldb    = 'prod3b-v2'       # Calibration database
rad      = 5.0               # FoV cut
idstart  = 0                 # Initial ID
t0       = 7671.0 * 86400.0  # s, corresponds to (1-1-2021)
emin     =   0.030           # TeV
emax     = 199.0             # TeV


# ================================================ #
# Characteristics of GPS KSP from Science with CTA #
# ================================================ #
GPS_STP = [{'name': 'Inner region',    'lmin': -60.0, 'lmax': 60.0,  'duration': 300.0,
           'array': 'South',  'invert': False, 'maxzenith': 50},
           {'name': 'Cygnus Perseus',  'lmin': 60.0,  'lmax': 150.0, 'duration': 180.0,
            'array': 'North', 'invert': False, 'maxzenith': 40}]
GPS_LTP = [{'name': 'Inner region',    'lmin': -60.0, 'lmax': 60.0,  'duration': 480.0,
            'array': 'South', 'invert': False, 'maxzenith': 50},
           {'name': 'Vela Carina',     'lmin': 240.0, 'lmax': 300.0, 'duration': 180.0,
            'array': 'South', 'invert': False, 'maxzenith': 40},
           {'name': 'Outer South',     'lmin': 210.0, 'lmax': 240.0, 'duration':  60.0,
            'array': 'South', 'invert': False, 'maxzenith': 40},
           {'name': 'Cygnus Perseus',  'lmin': 60.0,  'lmax': 150.0, 'duration': 270.0,
            'array': 'North', 'invert': False, 'maxzenith': 40},
           {'name': 'Anticentre',      'lmin': 150.0, 'lmax': 210.0, 'duration': 150.0,
            'array': 'North', 'invert': False, 'maxzenith': 40}]


# ====================== #
# Get pointing positions #
# ====================== #
def get_positions(xmin, xmax, nrows, step, fixed_xstep=False, invert=False):
    """
    Get pointing positions

    Parameters
    ----------
    xmin : float
        Longitude minimum (deg)
    xmax : float
        Longitude maximum (deg)
    nrows : integer
        Number of rows
    step : float
        Longitude stepsize (deg)
    fixed_xstep : bool, optional
        Fix longitude step size
    invert : bool, optional
        Invert longitude stepping

    Returns
    -------
    positions : list of tuple of float
        Pointing positions (x,y) in degrees
    """
    # Initialise positions
    positions = []

    # Determine ystep and ymin
    ystep = step * math.sqrt(0.75)
    ymin  = - ystep * (nrows - 1) / 2.0

    # Initialise y
    y = ymin

    # Loop over rows
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
        nx    = int((xmax - xmin) / xstep + 0.5)

        # Set x offset. For every second row the x position is displaced by
        # half a step size.
        if invert:
            if row % 2 == 0:
                x = xmax
            else:
                x = xmax - 0.5 * xstep
        else:
            if row % 2 == 0:
                x = xmin
            else:
                x = xmin + 0.5 * xstep

        # Append pointings
        for pnt in range(nx):

            # Check that we do not exceed maximum
            if x < xmin:
                print('WARNING: longitude %.3f < minimum %.3f' % (x, xmin))
            if x > xmax:
                print('WARNING: longitude %.3f > maximum %.3f' % (x, xmax))

            # Append pointing
            positions.append({'x': x, 'y': y})
            
            # Increment or decrement x
            if invert:
                x -= xstep
                if x < 0.0:
                    x += 360.0
            else:
                x += xstep
                if x >= 360.0:
                    x -= 360.0

        # Increment y
        y += ystep

    # Return positions
    return positions


# ======================== #
# Compute Sun zenith angle #
# ======================== #
def zenith_sun(time, array='South'):
    """
    Compute zenith angle of the Sun at a given time

    Parameters
    ----------
    time : `~gammalib.GTime()`
        Time for which the Sun zenith angle is to be computed
    array : string, optional
        Array site

    Returns
    -------
    zenith : float
        Sun zenith angle in (degrees)
    """
    # Compute Right Ascension and Declination of Sun in degrees and convert
    # the Declination into radians
    sun = gammalib.GSkyDir()
    sun.sun(time)

    # Get zenith angle
    zenith = zenith_dir(time, sun, array=array)

    # Return zenith angle
    return zenith


# ========================= #
# Compute Moon zenith angle #
# ========================= #
def zenith_moon(time, array='South'):
    """
    Compute zenith angle of the Moon at a given time

    Parameters
    ----------
    time : `~gammalib.GTime()`
        Time for which the Moon zenith angle is to be computed
    array : string, optional
        Array site

    Returns
    -------
    zenith : float
        Moon zenith angle in (degrees)
    """
    # Compute Right Ascension and Declination of Moon in degrees and convert
    # the Declination into radians
    moon = gammalib.GSkyDir()
    moon.moon(time)

    # Get zenith angle
    zenith = zenith_dir(time, moon, array=array)

    # Return zenith angle
    return zenith


# ================================== #
# Compute sky direction zenith angle #
# ================================== #
def zenith_dir(time, dir, array='South'):
    """
    Compute zenith angle of the Moon at a given time

    Parameters
    ----------
    time : `~gammalib.GTime()`
        Time for which the Moon zenith angle is to be computed
    dir : `~gammalib.GSkyDir()`
        Sky direction
    array : string, optional
        Array site

    Returns
    -------
    zenith : float
        Zenith angle in (degrees)
    """
    # Set geographic longitude and latitude
    if array == 'South':
        geolon =  79.4041
        geolat = -24.6272
    else:
        geolon =  17.8920
        geolat = +28.7622

    # Compute local apparent siderial time which is equal to the Right
    # Ascension that goes at a given time through the meridan
    last_deg = time.last(geolon) * 15.0

    # Compute hour angle of the Moon and it's cosine
    h     = dir.ra_deg() - last_deg
    cos_h = math.cos(h * gammalib.deg2rad)

    # Compute some sines and cosines
    sin_geolat = math.sin(geolat * gammalib.deg2rad)
    cos_geolat = math.cos(geolat * gammalib.deg2rad)
    sin_dirdec = math.sin(dir.dec())
    cos_dirdec = math.cos(dir.dec())

    # Compute cosine of zenith angle
    cos_zenith = sin_geolat * sin_dirdec + cos_geolat * cos_dirdec * cos_h

    # Compute zenith angle
    zenith = math.acos(cos_zenith) * gammalib.rad2deg

    # Return zenith angle
    return zenith


# ============== #
# Set start time #
# ============== #
def set_tmin_for_next_pointing(tmin, duration, array='South',
                               sunzenith=105.0, moonzenith=90.0,
                               skip_days=6.0):
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
    array : str, optional
        Array
    sunzenith : float, optional
        Sun zenith angle constraint (deg)
    moonzenith : float, optional
        Moon zenith angle constraint (deg)
    skip_days : float, optional
        Number of days to skip after Sun rise

    Returns
    -------
    tmin : float
        Start time of new pointing
    """
    # Set time
    tref = gammalib.GTimeReference(51544.5, 's', 'TT', 'LOCAL')
    time = gammalib.GTime(tmin, tref)

    # Add duration and 2 min for slew
    time += duration + 120.0

    # If the Sun or the Moon is above zenith angle limit then increment
    # time until the Sun and the Moon is again below the limit
    add_days = skip_days
    while True:
        sun  = zenith_sun(time, array=array)
        moon = zenith_moon(time, array=array)
        if sun < sunzenith and add_days > 0.0:
            time    += 86400.0 * add_days
            add_days = 0.0
        if sun < sunzenith or moon < moonzenith:
            time += 240.0  # Add 4 min and check again
        else:
            break

    # Convert back to seconds
    tmin = time.convert(tref)

    # Return start time
    return tmin


# ========================== #
# Set pointing pattern patch #
# ========================== #
def set_patch(name, tmin, xmin, xmax, duration, array, invert=False,
              maxzenith=40):
    """
    Set pointing pattern patch

    Parameters
    ----------
    name : string
        Name of pointing pattern patch
    tmin : float
        Start time (s)
    xmin : float
        Longitude minimum (deg)
    xmax : float
        Longitude maximum (deg)
    duration : float
        Pointing duration (s)
    array : string
        Array site ('South' or 'North')
    invert : bool, optional
        Invert longitude stepping
    maxzenith : float, optional
        Maximum zenith angle

    Returns
    -------
    obsdef : list of dict
        Observation definition
    """
    # Initialise observation definition
    obsdef = []

    # Set as many observations as allowed by total duration time and
    # observation length
    n_pnt = int(duration / obs_time)

    # If there are pointings then get positions
    if n_pnt > 0:
        xstep     = 2 * (xmax - xmin) / n_pnt
        pointings = get_positions(xmin, xmax, nrows, step, fixed_xstep=xstep,
                                  invert=invert)
    else:
        pointings = []

    # Set observations
    total_duration = 0
    for pnt in pointings:

        # Set positions, start time and duration
        obs = {'name': name, 'lon': pnt['x'], 'lat': pnt['y'], 'rad': rad, \
               'tmin': tmin, 'duration': obs_time * 3600.0, \
               'scheduled': False}

        # Update start time for next pointing
        tmin = set_tmin_for_next_pointing(tmin, obs_time * 3600.0)

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
        zenith             = abs(dec - geolat)
        obs['best_zenith'] = zenith

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
        obs['irf']   = irf

        # Set maximum zenith angle
        obs['maxzenith'] = maxzenith

        # Kluge for 'Inner region' negative longitudes
        if 'Inner region' in name:
            if pnt['x'] < 48.0 or pnt['x'] > 100.0:
                obs['maxzenith'] = 40.0

        # Append observation
        obsdef.append(obs)

        # Add to total duration
        total_duration += obs_time

    # Check that total duration is preserved
    if total_duration == duration:
        pass
    else:
        print(
            'WARNING: total duration of {} h and requested observation of {} do not match'.format(
                total_duration, duration))

    # Return
    return obsdef


# ================================= #
# Write observation definition file #
# ================================= #
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
    emin : float
        Minimum energy (TeV)
    emax : float
        Maximum energy (TeV)
    """
    # Open file
    f = open(filename, 'w')

    # Write header
    f.write('name,id,lon,lat,rad,tmin,duration,caldb,irf,emin,emax,zenith\n')

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
        f.write('%s,%6.6d,%8.4f,%8.4f,%8.4f,%.4f,%.4f,%s,%s,%.3f,%.1f,%.2f\n' %
                (obs['name'], obsid, obs['lon'], obs['lat'], obs['rad'],
                 obs['tmin'], obs['duration'],
                 obs['caldb'], obs['irf'],
                 emin, emax, obs['zenith']))

        # Increment identifier
        obsid += 1

    # Close file
    f.close()


# ===================== #
# Schedule observations #
# ===================== #
def schedule_observations_v1(obsdef, dl=1.0):
    """
    Schedule observations

    Parameters
    ----------
    obsdef : list of dict
        Observation definition
    dl : float, optional
        Longitude step size (deg)

    Returns
    -------
    obsdef : list of dict
        Scheduled observation definition
    """
    # Initialise statistics
    nobs     = len(obsdef)
    last_lon = 1000.0
    last_lat = 0.0

    # Initialise times for North and South
    time_south = t0
    time_north = t0

    # Initialise scheduled observation definitions
    scheduled_obsdef = []

    # Schedule observations until all are treated
    while len(scheduled_obsdef) < nobs:

        # Initialse number of observations found in one path
        n_obs_path = 0

        # Find next valid observation
        for obs in obsdef:
    
            # Consider next valid observation
            if obs['duration'] > 0.0:

                # If latitude is comparable then skip observation. This leads
                # to a toggle between the positive and negative latitudes of
                # the two-row pattern
                if (abs(obs['lat']-last_lat) < 0.1):
                    continue

                # If longitude difference is smaller than minimum longitude
                # step then skip observation
                if last_lon != 1000.0:
                    lon_step = last_lon - obs['lon']
                    if lon_step < dl:
                        continue

                # Set time
                if 'South' in obs['irf']:
                    obs['tmin'] = time_south
                else:
                    obs['tmin'] = time_north

                # Schedule observation
                scheduled_obsdef.append(obs)

                # Update last pointing
                last_lon    = obs['lon']
                last_lat    = obs['lat']
                n_obs_path += 1

                # Update start time for next pointing
                if 'South' in obs['irf']:
                    time_south = set_tmin_for_next_pointing(time_south,
                                                            obs_time*3600.0,
                                                            array='South')
                else:
                    time_north = set_tmin_for_next_pointing(time_north,
                                                            obs_time*3600.0,
                                                            array='North')

                # Debug
                #print(obs['lon'],obs['lat'])

                # Remove observation from list of valid observations
                obs['duration'] = 0.0

                # Break to restart from the beginning of the list
                break

        # If no observation was found then reinitialise some stuff
        if n_obs_path == 0:
            last_lon = 1000.0

    # Return scheduled observation
    return scheduled_obsdef


# ===================== #
# Schedule observations #
# ===================== #
def schedule_observations_v2(obsdef, tstart, dl_min=1.0, dl_max=5.0,
                             array='South', skip_days=6.0):
    """
    Schedule observations

    Parameters
    ----------
    obsdef : list of dict
        Observation definition
    tstart : float
        Start time (s)
    dl_min : float, optional
        Minimum longitude step (deg)
    dl_max : float, optional
        Maximum longitude step (deg)
    array : string
        Array site ('South' or 'North')
    skip_days : float, optional
        Number of days to skip after Sun rise

    Returns
    -------
    obsdef : list of dict
        Scheduled observation definition
    """
    # Initialise statistics
    nobs         = len(obsdef)
    last_lon     = 1000.0
    last_lat     = 0.0
    lon_min_step = dl_min
    lon_max_step = dl_max
    n_blocked    = 0

    # Initialise time
    time = set_tmin_for_next_pointing(tstart, 0.0, array=array, skip_days=skip_days)

    # Initialise scheduled observation definitions
    scheduled_obsdef = []

    # Schedule observations until all are treated
    while len(scheduled_obsdef) < nobs:

        # Find next valid observation with smallest zenith angle
        max_zenith =   0.0
        min_zenith = 180.0
        for obs in obsdef:
    
            # Consider next non-scheduled observation
            if not obs['scheduled']:

                # If latitude is comparable then skip observation. This leads
                # to a toggle between the positive and negative latitudes of
                # the two-row pattern
                if abs(obs['lat']-last_lat) < 0.1:
                    continue

                # If longitude difference is smaller than minimum longitude
                # step then skip observation
                if last_lon != 1000.0:
                    if obs['lon'] - last_lon < lon_min_step:
                        continue

                # If longitude difference is larger than maximum longitude
                # step the skip step
                if last_lon != 1000.0:
                    if obs['lon'] - last_lon > lon_max_step:
                        continue

                # Compute zenith angle
                dir    = gammalib.GSkyDir()
                dir.lb_deg(obs['lon'], obs['lat'])
                tref   = gammalib.GTimeReference(51544.5, 's', 'TT', 'LOCAL')
                gtime  = gammalib.GTime(time, tref)
                zenith = zenith_dir(gtime, dir, array=array)
                if zenith > obs['maxzenith']:
                    continue
                if zenith < min_zenith:
                    min_zenith = zenith
                    max_zenith = obs['maxzenith']
                    min_obs    = obs

        # If an observation with an acceptable zenith angle was found then
        # use it
        if min_zenith < max_zenith:

            # Set time and duration
            min_obs['tmin']     = time

            # Set zenith angle
            min_obs['zenith'] = min_zenith

            # Assign IRF zenith
            if min_zenith < 30.0:
                irfz = 20
            elif min_zenith < 50.0:
                irfz = 40
            else:
                irfz = 60
            irf = '{}_z{}_50h'.format(array, irfz)

            # Set IRF information
            min_obs['irf'] = irf

            # Schedule observation
            scheduled_obsdef.append(min_obs)

            # Update last pointing
            last_lon = min_obs['lon']
            last_lat = min_obs['lat']

            # Dump
            gtime  = gammalib.GTime(time, tref)
            print('%8.3f %8.3f %8.3f  %4d/%4d  %s %s' %
                  (min_obs['lon'], min_obs['lat'], min_zenith,
                   len(scheduled_obsdef),nobs,
                   gtime.utc(), min_obs['name']))

            # Signal that observation was scheduled
            min_obs['scheduled'] = True

        # ... otherwise remove the longitude and latitude constraints
        else:
            # Remove constraints
            last_lon = 1000.0
            last_lat =    0.0

            # Log number of blocks
            n_blocked += 1

            # If blocked for more than 1000 iterations signal
            if n_blocked > 1000:
                print('WARNING: Scheduling blocked' % (n_blocked, min_zenith, max_zenith))
                for obs in obsdef:
                    if not obs['scheduled']:
                        print('- %s %.2f %.2f %.2f' % (obs['name'], obs['lon'],
                              obs['lat'], obs['maxzenith']))
                break

        # Increment time for next pointing
        time = set_tmin_for_next_pointing(time, obs_time*3600.0, array=array,
                                          skip_days=skip_days)

    # Return scheduled observation
    return scheduled_obsdef


# ===================== #
# Schedule observations #
# ===================== #
def schedule_observations_v3(obsdef, tstart, dl_min=1.0, dl_max=5.0,
                             array='South', skip_days=6.0, dzenith_max=5.0):
    """
    Schedule observations

    Parameters
    ----------
    obsdef : list of dict
        Observation definition
    tstart : float
        Start time (s)
    dl_min : float, optional
        Minimum longitude step (deg)
    dl_max : float, optional
        Maximum longitude step (deg)
    array : string
        Array site ('South' or 'North')
    skip_days : float, optional
        Number of days to skip after Sun rise

    Returns
    -------
    obsdef : list of dict
        Scheduled observation definition
    """
    # Initialise statistics
    nobs         = len(obsdef)
    last_lon     = 1000.0
    last_lat     = 0.0
    lon_min_step = dl_min
    lon_max_step = dl_max
    n_blocked    = 0

    # Initialise time
    time = set_tmin_for_next_pointing(tstart, 0.0, array=array, skip_days=skip_days)

    # Initialise scheduled observation definitions
    scheduled_obsdef = []

    # Schedule observations until all are treated
    while len(scheduled_obsdef) < nobs:

        # Find next valid observation with smallest zenith angle
        max_zenith  =   0.0
        min_zenith  = 180.0
        min_dzenith = 180.0
        for obs in obsdef:
    
            # Consider next non-scheduled observation
            if not obs['scheduled']:

                # If latitude is comparable then skip observation. This leads
                # to a toggle between the positive and negative latitudes of
                # the two-row pattern
                if abs(obs['lat']-last_lat) < 0.1:
                    continue

                # If longitude difference is smaller than minimum longitude
                # step then skip observation
                if last_lon != 1000.0:
                    if obs['lon'] - last_lon < lon_min_step:
                        continue

                # If longitude difference is larger than maximum longitude
                # step the skip step
                if last_lon != 1000.0:
                    if obs['lon'] - last_lon > lon_max_step:
                        continue

                # Compute zenith angle
                dir    = gammalib.GSkyDir()
                dir.lb_deg(obs['lon'], obs['lat'])
                tref   = gammalib.GTimeReference(51544.5, 's', 'TT', 'LOCAL')
                gtime  = gammalib.GTime(time, tref)
                zenith = zenith_dir(gtime, dir, array=array)

                # Compute zenith distance from culmination
                dzenith = zenith - obs['best_zenith']

                # If zenith distance from culmination is worse by dzenith_max
                # deg then skip pointing
                if dzenith > dzenith_max:
                    continue

                # Keep pointing with smallest distance from culmination
                if dzenith < min_dzenith:
                    min_dzenith = dzenith
                    min_zenith  = zenith
                    max_zenith  = obs['maxzenith']
                    min_obs     = obs

        # If an observation with an acceptable zenith angle was found then
        # use it
        if min_dzenith < dzenith_max:

            # Set time and duration
            min_obs['tmin'] = time

            # Set zenith angle
            min_obs['zenith'] = min_zenith

            # Assign IRF zenith
            if min_zenith < 30.0:
                irfz = 20
            elif min_zenith < 50.0:
                irfz = 40
            else:
                irfz = 60
            irf = '{}_z{}_50h'.format(array, irfz)

            # Set IRF information
            min_obs['irf'] = irf

            # Schedule observation
            scheduled_obsdef.append(min_obs)

            # Update last pointing
            last_lon = min_obs['lon']
            last_lat = min_obs['lat']

            # Dump
            gtime  = gammalib.GTime(time, tref)
            print('%8.3f %8.3f %8.3f  %4d/%4d  %s %s' %
                  (min_obs['lon'], min_obs['lat'], min_zenith,
                   len(scheduled_obsdef),nobs,
                   gtime.utc(), min_obs['name']))

            # Signal that observation was scheduled
            min_obs['scheduled'] = True

        # ... otherwise remove the longitude and latitude constraints
        else:
            # Remove constraints
            last_lon = 1000.0
            last_lat =    0.0

            # Log number of blocks
            n_blocked += 1

            # If blocked for more than 5000 iterations signal
            if n_blocked > 5000:
                print('WARNING: Scheduling blocked')
                for obs in obsdef:
                    if not obs['scheduled']:
                        print('- %s %.2f %.2f %.2f %.2f' % (obs['name'], obs['lon'],
                              obs['lat'], obs['best_zenith'], obs['maxzenith']))
                break

        # Increment time for next pointing
        time = set_tmin_for_next_pointing(time, obs_time*3600.0, array=array,
                                          skip_days=skip_days)

    # Return scheduled observation
    return scheduled_obsdef


# ================================ #
# Make pointings (STP/LTP version) #
# ================================ #
def make_pointings_stp_ltp():
    """
    """
    # Initialise observation definitions
    obsdef           = []
    obsdef_south_stp = []
    obsdef_south_ltp = []
    obsdef_north_stp = []
    obsdef_north_ltp = []

    # Southern survey
    print('Generate South observations:')
    print('============================')
    tmin = t0
    for patch in GPS_STP:
        if patch['array'] == 'South':

            # Set patch
            obs = set_patch(patch['name'], tmin, patch['lmin'], patch['lmax'],
                            patch['duration'],
                            patch['array'],
                            invert=patch['invert'],
                            maxzenith=patch['maxzenith'])

            # Add patch
            obsdef_south_stp.extend(obs)

            # Update start time
            tmin = set_tmin_for_next_pointing(obsdef_south_stp[-1]['tmin'],
                                              obsdef_south_stp[-1]['duration'])
    for patch in GPS_LTP:
        if patch['array'] == 'South':

            # Set patch
            obs = set_patch(patch['name'], tmin, patch['lmin'], patch['lmax'],
                            patch['duration'],
                            patch['array'],
                            invert=patch['invert'],
                            maxzenith=patch['maxzenith'])

            # Add patch
            obsdef_south_ltp.extend(obs)

            # Update start time
            tmin = set_tmin_for_next_pointing(obsdef_south_ltp[-1]['tmin'],
                                              obsdef_south_ltp[-1]['duration'])

    # Northern survey
    print('Generate North observations:')
    print('============================')
    tmin = t0
    for patch in GPS_STP:
        if patch['array'] == 'North':

            # Set patch
            obs = set_patch(patch['name'], tmin, patch['lmin'], patch['lmax'],
                            patch['duration'],
                            patch['array'],
                            invert=patch['invert'],
                            maxzenith=patch['maxzenith'])

            # Add patch
            obsdef_north_stp.extend(obs)

            # Update start time
            tmin = set_tmin_for_next_pointing(obsdef_north_stp[-1]['tmin'],
                                              obsdef_north_stp[-1]['duration'])
    for patch in GPS_LTP:
        if patch['array'] == 'North':

            # Set patch
            obs = set_patch(patch['name'], tmin, patch['lmin'], patch['lmax'],
                            patch['duration'],
                            patch['array'],
                            invert=patch['invert'],
                            maxzenith=patch['maxzenith'])

            # Add patch
            obsdef_north_ltp.extend(obs)

            # Update start time
            tmin = set_tmin_for_next_pointing(obsdef_north_ltp[-1]['tmin'],
                                              obsdef_north_ltp[-1]['duration'])

    # Compute start times
    t0_stp = t0
    t0_ltp = t0 + 2.0*365.25*86400.0

    # Schedule observations
    print('Schedule South observations:')
    print('============================')
    obsdef_south_stp = schedule_observations_v2(obsdef_south_stp, t0_stp,
                                                array='South',
                                                dl_min=2.0, dl_max=10.0,
                                                skip_days=3.0)
    obsdef_south_ltp = schedule_observations_v2(obsdef_south_ltp, t0_ltp,
                                                array='South',
                                                dl_min=2.0, dl_max=10.0,
                                                skip_days=7.0)
    obsdef.extend(obsdef_south_stp)
    obsdef.extend(obsdef_south_ltp)

    print('Schedule North observations:')
    print('============================')
    obsdef_north_stp = schedule_observations_v2(obsdef_north_stp, t0_stp,
                                                array='North',
                                                dl_min=2.0, dl_max=10.0,
                                                skip_days=4.0)
    obsdef_north_ltp = schedule_observations_v2(obsdef_north_ltp, t0_ltp,
                                                array='North',
                                                dl_min=2.0, dl_max=10.0,
                                                skip_days=10.0)
    obsdef.extend(obsdef_north_stp)
    obsdef.extend(obsdef_north_ltp)

    # Return obsdef
    return obsdef


# ========================== #
# Make pointings (Y version) #
# ========================== #
def make_pointings_y(shift=False, noskip=False):
    """
    """
    # Setup yearly time fraction
    time_fraction = [0.6, 0.4, 0.15, 0.2, 0.2, 0.15, 0.1, 0.1, 0.05, 0.05]

    # Setup number of days to skip
    if noskip:
        n_skip_south = [0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0]
        n_skip_north = [0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0]
    elif shift:
        n_skip_south = [1.5, 3.0,  3.0,  1.5,  2.5,  5.0,  7.0, 10.0, 12.0, 14.0]
        n_skip_north = [4.0, 5.0,  5.0,  4.0,  4.0,  5.0,  6.0,  8.0, 12.0, 13.0]
    else:
        n_skip_south = [3.5, 3.5,  6.5,  6.5,  6.5,  6.5,  6.5,  6.5,  6.5,  6.5]
        n_skip_north = [6.0, 6.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]

    # Assign yearly bunches
    gps_yearly = []
    for i in range(2):
        gps_yearly.append(copy.deepcopy(GPS_STP))
    for i in range(8):
        gps_yearly.append(copy.deepcopy(GPS_LTP))

    # Assign yearly observing times
    for year in range(10):
        if year < 2:
            for patch in gps_yearly[year]:
                patch['duration'] *= time_fraction[year]
        else:
            for patch in gps_yearly[year]:
                patch['duration'] *= time_fraction[year]

    # Check time fractions
    fraction_stp = 0.0
    fraction_ltp = 0.0
    for i in range(2):
        fraction_stp += time_fraction[i]
    for i in range(8):
        fraction_ltp += time_fraction[i+2]
    print('Time fraction STP ............: %.4f' % fraction_stp)
    print('Time fraction LTP ............: %.4f' % fraction_ltp)

    # Compute yearly durations and total
    south_time = 0.0
    north_time = 0.0
    for year in range(10):
        for patch in gps_yearly[year]:
            if patch['array'] == 'South':
                south_time += patch['duration']
            else:
                north_time += patch['duration']
    print('Total observing time South ...: %.2f' % south_time)
    print('Total observing time North ...: %.2f' % north_time)

    # Initialise observation definitions
    obsdef              = []
    obsdef_south_yearly = []
    obsdef_north_yearly = []

    # Southern survey
    print('Generate South observations:')
    print('============================')
    tmin = t0
    for year in range(10):

        # Initialise yearly list
        obsdef_yearly = []

        # Loop over all patched
        for patch in gps_yearly[year]:
            if patch['array'] == 'South':

                # Set patch
                obs = set_patch(patch['name'], tmin, patch['lmin'], patch['lmax'],
                                patch['duration'],
                                patch['array'],
                                invert=patch['invert'],
                                maxzenith=patch['maxzenith'])

                # Add patch
                obsdef_yearly.extend(obs)

                # Update start time
                tmin = set_tmin_for_next_pointing(obsdef_yearly[-1]['tmin'],
                                                  obsdef_yearly[-1]['duration'])

        # Add observation definition to yearly list
        obsdef_south_yearly.append(obsdef_yearly)

    # Northern survey
    print('Generate North observations:')
    print('============================')
    tmin = t0
    for year in range(10):

        # Initialise yearly list
        obsdef_yearly = []

        # Loop over all patched
        for patch in gps_yearly[year]:
            if patch['array'] == 'North':

                # Set patch
                obs = set_patch(patch['name'], tmin, patch['lmin'], patch['lmax'],
                                patch['duration'],
                                patch['array'],
                                invert=patch['invert'],
                                maxzenith=patch['maxzenith'])

                # Add patch
                obsdef_yearly.extend(obs)

                # Update start time
                tmin = set_tmin_for_next_pointing(obsdef_yearly[-1]['tmin'],
                                                  obsdef_yearly[-1]['duration'])

        # Add observation definition to yearly list
        obsdef_north_yearly.append(obsdef_yearly)

    # Compute start times
    t0_stp = t0
    t0_ltp = t0 + 2.0*365.25*86400.0

    # Schedule observations
    print('Schedule South observations:')
    print('============================')
    for year in range(10):

        # Set start time for year
        t0_year = t0 + year*365.25*86400.0

        # Shift time
        if shift:
            if year < 2:
                t0_year += (1-year)*29.5*86400.0
            else:
                t0_year += (10-year)*29.5*86400.0*2.0 - 365.25*86400.0

        # Set days to skip
        skip_days = n_skip_south[year]

        # Schedule observations
        obsdef_scheduled = schedule_observations_v3(obsdef_south_yearly[year],
                                                    t0_year,
                                                    array='South',
                                                    dl_min=2.0,
                                                    dl_max=10.0,
                                                    skip_days=skip_days)

        # Add scheduled observations
        obsdef.extend(obsdef_scheduled)

    print('Schedule North observations:')
    print('============================')
    for year in range(10):

        # Set start time for year
        t0_year = t0 + year*365.25*86400.0

        # Shift time
        if shift:
            if year < 2:
                t0_year += (1-year)*29.5*86400.0
            else:
                t0_year += (10-year)*29.5*86400.0*2.0 - 365.25*86400.0

        # Set days to skip
        skip_days = n_skip_north[year]

        # Schedule observations
        obsdef_scheduled = schedule_observations_v3(obsdef_north_yearly[year],
                                                    t0_year,
                                                    array='North',
                                                    dl_min=2.0,
                                                    dl_max=10.0,
                                                    skip_days=skip_days)

        # Add scheduled observations
        obsdef.extend(obsdef_scheduled)

    # Return obsdef
    return obsdef


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Make pointings
    #obsdef = make_pointings_stp_ltp()
    obsdef = make_pointings_y(shift=True, noskip=False)

    # Write observation definition file
    write_obsdef('GPS_pointings.dat', obsdef, idstart, emin, emax)
