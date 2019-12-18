#! /usr/bin/env python
# ==========================================================================
# Simulate binary for Galactic Plane survey
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
import shutil
import glob
import gammalib
import ctools
import cscripts
from multiprocessing import Pool


# ======================= #
# Enter working directory #
# ======================= #
def enter_wd(name, create=False):
    """
    Enter working directory

    Parameters
    ----------
    name : string
        Source name
    create : bool, optional
        Create directory tree

    Returns
    -------
    cwd : string
        Current working directory
    """
    # Get current working directory
    cwd = os.getcwd()

    # Create working directory filename
    fname1   = gammalib.replace_segment(gammalib.tolower(name), ' ', '')
    fname2   = gammalib.replace_segment(fname1, '+', 'p')
    filename = gammalib.replace_segment(fname2, '-', 'm')
    dirname  = 'gps/%s' % (filename)

    # Create directory tree
    if create:
        try:
            os.makedirs('gps')
        except:
            pass
        try:
            os.makedirs(dirname)
        except:
            pass

    # Step into working directory
    os.chdir(dirname)

    # Return
    return cwd


# ======================= #
# Exist working directory #
# ======================= #
def exit_wd(cwd):
    """
    Enter working directory

    Parameters
    ----------
    cwd : string
        Current working directory
    """
    # Step out of working directory
    os.chdir(cwd)

    # Return
    return


# ================================== #
# Create observation definition file #
# ================================== #
def create_obsdef(pntfile, name, lon, lat, rad=3.0):
    """
    Create observation definition file for binary

    Parameters
    ----------
    pntfile : string
        Pointing file
    name : string
        Source name
    lon : float
        Galactic longitude of binary (deg)
    lat : float
        Galactic latitude of binary (deg)
    rad : float, optional
        Maximum offset angle (deg)
    """
    # Enter working directory
    cwd = enter_wd(name, True)

    # Run make_pointings.py if pointings file does not exist
    if not os.path.isfile(pntfile):
        print('Create GPS pointings')
        os.system(os.path.expandvars('${CTAGITROOT}/simulation/gps/make_pointings.py'))
        #os.system('cp ../../gps.dat .')

    # Run csobsdef if observation definition file does not exist
    if not os.path.isfile('gps_obs.xml'):
        print('Create observation definition file')
        obsdef = cscripts.csobsdef()
        obsdef['inpnt']   = pntfile
        obsdef['outobs']  = 'gps_obs.xml'
        obsdef['rad']     = 5.0
        obsdef['logfile'] = 'gps_csobsdef.log'
        obsdef.logFileOpen()
        obsdef.execute()

    # Run csobsselect if observation definition file does not exist
    if not os.path.isfile('gps_obs_selected.xml'):
        print('Select pointings from observation definition file')
        obsselect = cscripts.csobsselect()
        obsselect['inobs']     = 'gps_obs.xml'
        obsselect['outobs']    = 'gps_obs_selected.xml'
        obsselect['pntselect'] = 'CIRCLE'
        obsselect['coordsys']  = 'GAL'
        obsselect['glon']      = lon
        obsselect['glat']      = lat
        obsselect['rad']       = rad
        obsselect['tmin']      = 'NONE'
        obsselect['tmax']      = 'NONE'
        obsselect['logfile']   = 'gps_csobsselect.log'
        obsselect.logFileOpen()
        obsselect.execute()

    # Exit working directory
    exit_wd(cwd)

    # Return
    return


# ============ #
# Setup models #
# ============ #
def setup_models(basedir, name, lc=True):
    """
    Setup model container for simulation

    Parameters
    ----------
    basedir : string
        Base directory
    name : string
        Name of source component

    Returns
    -------
    models : `~gammalib.GModels()`
        Model container
    """
    # Initialise model container
    models = gammalib.GModels()

    # Extract binary component
    binaries = gammalib.GModels(basedir+'/1dc/models/model_galactic_binaries.xml')
    binary   = binaries[name]

    # Optionally remove lightcurve
    if not lc:
        binary.temporal(gammalib.GModelTemporalConst())

    # Append binary to model container
    models.append(binary)

    # Append background model to container
    models.extend(gammalib.GModels(basedir+'/1dc/models/model_bkg.xml'))

    # Return model container
    return models


# ==================== #
# Simulate observation #
# ==================== #
def simulate_observation(obs, name, models, suffix, edisp=True, deadc=0.98,
                         debug=False):
    """
    Simulate events

    Parameters
    ----------
    obs : `~gammalib.GCTAObservation`
        CTA observation
    name : string
        Source name
    models : `~gammalib.GModels`
        Models
    suffix : string
        Suffix
    edisp : bool, optional
        Enable energy dispersion in simulation
    deadc : float, optional
        Deadtime correction factor
    debug : bool, optional
        Debug simulation
    """
    # Create directories
    try:
        os.makedirs('data')
    except:
        pass
    try:
        os.makedirs('log')
    except:
        pass

    # Set output filenames
    outevents = 'data/gps_baseline_%s%s.fits' % (obs.id(), suffix)
    logfile   = 'log/gps_ctobssim_%s%s.log'   % (obs.id(), suffix)

    # Continue only if output file does not exist
    if not os.path.isfile(outevents+'.gz'):

        # Dump processing
        print('Create %s' % outevents)

        # Create empty observation container
        container = gammalib.GObservations()

        # Append curent observation to container
        container.append(obs)

        # Attach models to container
        container.models(models)

        # Simulate events
        sim = ctools.ctobssim(container)
        sim['seed']      = int(obs.id())
        sim['edisp']     = edisp
        sim['outevents'] = outevents
        sim['deadc']     = deadc
        sim['logfile']   = logfile
        sim.logFileOpen()
        sim.execute()

        # Compress FITS file
        os.system('gzip %s' % outevents)

    # ... otherwise notify skipping
    else:
        print('Skip %s (exists already)' % outevents)

    # Return
    return


# ===================== #
# Simulate observations #
# ===================== #
def simulate_observations(name, processes=1, lc=True):
    """
    Simulate observations for binary

    Parameters
    ----------
    name : string
        Binary name
    processes : int, optional
        Number of parallel processes
    lc : bool, optional
        Use light curve
    """
    # Dump header
    print('Simulate observations')

    # Enter working directory
    cwd = enter_wd(name)

    # Setup models
    models = setup_models('/project-data/cta/data', name, lc=lc)

    # Save models
    if lc:
        models.save('gps_models.xml')
    else:
        models.save('gps_models_const.xml')

    # Load observations
    observations = gammalib.GObservations('gps_obs_selected.xml')

    # Append models to observation container
    observations.models(models)

    # Dump header
    print('Model setup completed')

    # Create a pool of processes (need to do this after loading the observations)
    if processes > 1:
        pool = Pool(processes=processes)

    # Loop over all observations
    for obs in observations:

        # Define process arguments
        if lc:
            args = (obs, name, models, '')
        else:
            args = (obs, name, models, '_const')

        # Start processes
        if processes > 1:
            pool.apply_async(simulate_observation, args)
            print('Schedule observation "%s"' % (obs.id()))
        else:
            simulate_observation(*args)

    # Wait for all processes to finish
    if processes > 1:
        pool.close()
        pool.join()

    # Exit working directory
    exit_wd(cwd)

    # Return
    return


# ============================================ #
# Collect data and write observation container #
# ============================================ #
def collect_data(name, lc=True):
    """
    Collect data and write observation container

    Parameters
    ----------
    name : string
        Source name
    lc : bool, optional
        Use light curve
    """
    # Dump header
    print('Generate observation definition file')

    # Enter working directory
    cwd = enter_wd(name)

    # Set filenames
    if lc:
        suffix = ''
    else:
        suffix = '_const'
    outxml = 'gps_obs_simulated%s.xml' % (suffix)

    # Continue only if XML file does not exist
    if not os.path.isfile(outxml):

        # Create empty XML file
        xmlout = gammalib.GXml()
        list   = xmlout.append('observation_list title="observation list"')

        # Load XML file
        xml = gammalib.GXml('gps_obs_selected.xml')

        # Get observation list
        obslist = xml.element('observation_list', 0)

        # Get number of observations
        nobs = obslist.elements('observation')

        # Loop over all observations
        for iobs in range(nobs):

            # Get current observation
            obs = obslist.element('observation', iobs)

            # Get observation ID
            id = obs.attribute('id')

            # Search for a corresponding event file
            files = glob.glob(os.path.expandvars('data/*%s%s.fits.gz' % (id, suffix)))

            # If we found one file then get the relative filename and append a
            # "EventList" parameter to the XML file
            if len(files) == 1:

                # Get event filename
                eventfile = 'data/'+os.path.basename(files[0])

                # Append XML element
                obs.append('parameter name="EventList" file="%s"' % eventfile)

                # Copy element to output file
                list.append(obs)

            # ... otherwise notify that an event file is missing
            else:
                print('No event file for observation %s' % id)

        # Save XML file
        xmlout.save(outxml)

    # Exit working directory
    exit_wd(cwd)

    # Return
    return


# ================ #
# Generate sky map #
# ================ #
def generate_skymap(name, lon, lat, emin=0.030, emax=199.0):
    """
    Generate sky map

    Parameters
    ----------
    name : string
        Source name
    lon : float
        Galactic longitude of binary (deg)
    lat : float
        Galactic latitude of binary (deg)
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    """
    # Dump header
    print('Generate sky map')

    # Enter working directory
    cwd = enter_wd(name)

    # Continue only if skymap does not exist
    if not os.path.isfile('skymap.fits'):

        # Setup task parameters
        skymap = ctools.ctskymap()
        skymap['inobs']       = 'gps_obs_simulated.xml'
        skymap['emin']        = emin
        skymap['emax']        = emax
        skymap['nxpix']       = 200
        skymap['nypix']       = 200
        skymap['binsz']       = 0.02
        skymap['coordsys']    = 'GAL'
        skymap['proj']        = 'CAR'
        skymap['xref']        = lon
        skymap['yref']        = lat
        skymap['bkgsubtract'] = 'RING'
        skymap['roiradius']   = 0.02
        skymap['inradius']    = 0.60
        skymap['outradius']   = 0.80
        skymap['iterations']  = 3
        skymap['threshold']   = 5.0
        skymap['inexclusion'] = 'NONE'
        skymap['outmap']      = 'skymap.fits'
        skymap['logfile']     = 'ctskymap.log'
        skymap.logFileOpen()

        # Generate sky map
        skymap.execute()

    # Exit working directory
    exit_wd(cwd)

    # Return
    return


# ============ #
# Phase events #
# ============ #
def phase_events(name, lon, lat, rad=0.2, lc=True):
    """
    Phase events

    Parameters
    ----------
    name : string
        Source name
    lon : float
        Galactic longitude of binary (deg)
    lat : float
        Galactic latitude of binary (deg)
    rad : float, optional
        Angular cut (deg)
    lc : bool, optional
        Use light curve
    """
    # Dump header
    print('Phase events')

    # Enter working directory
    cwd = enter_wd(name)

    # Set filenames
    if lc:
        suffix = ''
    else:
        suffix = '_const'
    inxml   = 'gps_obs_simulated%s.xml'       % (suffix)
    outxml  = 'gps_obs_phased%s.xml'          % (suffix)
    prefix  = 'phased%s/'                     % (suffix)
    outxml2 = 'gps_obs_selected_phased%s.xml' % (suffix)
    prefix2 = 'selected_phased%s/'            % (suffix)

    # Continue only if phased events do not exist
    if not os.path.isfile(outxml):

        # Create directory
        try:
            os.makedirs('phased%s' % suffix)
        except:
            pass

        # Setup task parameters
        phase = ctools.ctphase()
        phase['inobs']   = inxml
        phase['outobs']  = outxml
        phase['prefix']  = prefix
        phase['inmodel'] = 'gps_models.xml'
        phase['srcname'] = name
        phase.logFileOpen()

        # Phase events
        phase.execute()

    # Continue only if selected phased events do not exist
    if not os.path.isfile(outxml2):

        # Create directory
        try:
            os.makedirs('selected_phased%s' % suffix)
        except:
            pass

        # Set source position
        dir = gammalib.GSkyDir()
        dir.lb_deg(lon, lat)

        # Setup task parameters
        select = ctools.ctselect()
        select['inobs']  = outxml
        select['outobs'] = outxml2
        select['prefix'] = prefix2
        select['usepnt'] = False
        select['ra']     = dir.ra_deg()
        select['dec']    = dir.dec_deg()
        select['rad']    = rad
        select['emin']   = 'NONE'
        select['emax']   = 'NONE'
        select['tmin']   = 'NONE'
        select['tmax']   = 'NONE'
        select.logFileOpen()

        # Select events
        select.execute()

    # Exit working directory
    exit_wd(cwd)

    # Return
    return


# ============== #
# Check binaries #
# ============== #
def check_binary_simulation():
    """
    Check binaries
    """
    # Set binaries
    binaries = [{'name': 'LS 5039',          'lon':  16.88,   'lat': -1.29},
                {'name': 'LS I 61+303',      'lon': 135.6753, 'lat':  1.0861},
                {'name': 'PSR B1259-63',     'lon': 304.1836, 'lat': -0.9916},
                {'name': 'HESS J1832-093',   'lon':  22.4768, 'lat': -0.1539},
                {'name': 'HESS J0632+057',   'lon': 205.6665, 'lat': -1.4395},
                {'name': 'HESS J1018-589 A', 'lon': 284.3614, 'lat': -1.6784}]

    # Loop over binaries
    for binary in binaries:

        # Extract parameters
        name = binary['name']
        lon  = binary['lon']
        lat  = binary['lat']

        # Create observation definition file
        create_obsdef('GPS_pointings.dat', name, lon, lat)
        #create_obsdef('gps.dat', name, lon, lat)

        # Simulate observations
        simulate_observations(name, processes=4)
        simulate_observations(name, processes=4, lc=False)

        # Collect data
        collect_data(name)
        collect_data(name, lc=False)

        # Generate sky map
        generate_skymap(name, lon, lat)

        # Phase-fold events
        phase_events(name, lon, lat)
        phase_events(name, lon, lat, lc=False)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Dump header
    print('**********************************')
    print('* Check binary in GPS simulation *')
    print('**********************************')

    # Check binaries
    check_binary_simulation()

    # We are done
    print('... done.')
