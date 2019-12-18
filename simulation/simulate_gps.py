#! /usr/bin/env python
# ==========================================================================
# Simulate Galactic Plane survey
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


# ================ #
# Create directory #
# ================ #
def makedirs(dirname):
    """
    Create directory

    Parameters
    ----------
    dirname : str
        Directory name
    """
    #
    try:
        os.makedirs(dirname)
    except:
        pass

    # Return
    return


# ===================== #
# Create directory tree #
# ===================== #
def create_directory_tree(delete=False):
    """
    Create directory tree

    Parameters
    ----------
    delete : bool, optional
        Delete an existing output directory
    """
    # Optionally delete old directory tree
    if delete:
        os.system('rm -rf gps')

    # Create directory tree
    makedirs('gps')
    makedirs('gps/data')
    makedirs('gps/log')
    makedirs('gps/obs')
    makedirs('gps/models')

    # Copy model components
    for file in glob.glob(os.path.expandvars('github/cta-gps-simulation-paper/skymodel/output/*.xml')):
        if not os.path.isfile('gps/models/%s' % os.path.basename(file)):
            print('Copy %s' % file)
            shutil.copy(file, 'gps/models')
    for file in glob.glob(os.path.expandvars('github/cta-gps-simulation-paper/skymodel/output/*.fits')):
        if not os.path.isfile('gps/models/%s' % os.path.basename(file)):
            print('Copy %s' % file)
            shutil.copy(file, 'gps/models')
    for file in glob.glob(os.path.expandvars('github/cta-gps-simulation-paper/skymodel/output/*.txt')):
        if not os.path.isfile('gps/models/%s' % os.path.basename(file)):
            print('Copy %s' % file)
            shutil.copy(file, 'gps/models')
    for file in glob.glob(os.path.expandvars('github/cta-gps-simulation-paper/skymodel/iem/model_iem.xml')):
        if not os.path.isfile('gps/models/%s' % os.path.basename(file)):
            print('Copy %s' % file)
            shutil.copy(file, 'gps/models')
    for file in glob.glob(os.path.expandvars('github/cta-gps-simulation-paper/skymodel/iem/IEM_base.fits.gz')):
        if not os.path.isfile('gps/models/%s' % os.path.basename(file)):
            print('Copy %s' % file)
            shutil.copy(file, 'gps/models')

    # Return
    return


# ================================== #
# Create observation definition file #
# ================================== #
def create_obsdef(event_class):
    """
    Create observation definition file for binary

    Parameters
    ----------
    pntfile : string
        Pointing file
    event_class : dict
        Event type dictionary
    """
    # Get current working directory
    cwd = os.getcwd()

    # Step into data directory
    os.chdir('gps/obs')

    # Set filenames
    pntfile = 'GPS_pointings.dat'
    outobs  = 'gps_obs_%s_obsdef.xml' % (event_class['name'])
    logfile = '../log/gps_obs_%s.log' % (event_class['name'])

    # Run make_pointings.py if pointings file does not exist
    if not os.path.isfile(pntfile):
        print('Create GPS pointings for "%s" event class' % event_class['name'])
        os.system(os.path.expandvars('${CTAGITROOT}/simulation/gps/make_pointings.py'))

    # Substitute IRF names
    inpnt = substitute_irfs(pntfile, event_class)

    # Run csobsdef if observation definition file does not exist
    if not os.path.isfile(outobs):
        print('Create observation definition file for "%s" event class' % event_class['name'])
        obsdef = cscripts.csobsdef()
        obsdef['inpnt']   = inpnt
        obsdef['outobs']  = outobs
        obsdef['rad']     = 5.0
        obsdef['logfile'] = logfile
        obsdef.logFileOpen()
        obsdef.execute()

    # Step out of working directory
    os.chdir(cwd)

    # Return
    return


# =============================== #
# Substitute IRFs for event class #
# =============================== #
def substitute_irfs(pntfile, event_class):
    """
    Substitute IRFs for event class

    Parameters
    ----------
    pntfile : string
        Pointing file
    event_class : dict
        Event class dictionary
    """
    # Read file
    f = open(pntfile, 'r')
    lines = f.readlines()
    f.close()

    # Build new filename
    filename = os.path.splitext(pntfile)[0] + '_%s.dat' % (event_class['name'])

    # Replace IRFs in all lines
    for i, line in enumerate(lines):
        lines[i] = line.replace('_50h', '_%s' % event_class['suffix'], 1)

    # Write file
    f = open(filename, 'w')
    for line in lines:
        f.write(line)
    f.close()

    # Return filename
    return filename


# ===================== #
# Simulate observations #
# ===================== #
def simulate_observations(event_class, processes=1):
    """
    Simulate observations for event class

    Parameters
    ----------
    event_class : dict
        Event class
    processes : int, optional
        Number of parallel processes
    """
    # Dump header
    print('Simulate observations for "%s" event class' % (event_class['name']))

    # Get current working directory
    cwd = os.getcwd()

    # Set filenames
    dirname = 'gps/data/%s'                     % (event_class['name'])
    obsname = '../../obs/gps_obs_%s_obsdef.xml' % (event_class['name'])
    modname = '../../models/models_gps.xml'
    iemname = '../../models/model_iem.xml'

    # Create directory for data
    makedirs(dirname)

    # Step into data directory
    os.chdir(dirname)

    # Load observations
    observations = gammalib.GObservations(obsname)

    # Load models
    models     = gammalib.GModels(modname)
    models_iem = gammalib.GModels(iemname)
    models.extend(models_iem)

    # Create a pool of processes (need to do this after loading the observations)
    if processes > 1:
        pool = Pool(processes=processes)

    # Loop over all observations
    for obs in observations:

        # Define process arguments
        args = (obs, event_class['name'], models)

        # Start processes
        if processes > 1:
            pool.apply_async(simulate_observation, args)
            print('Schedule observation "%s" for "%s" event class' %
                  (obs.id(), event_class['name']))
        else:
            simulate_observation(*args)

        # For testing
        #if obs.id() == '000002':
        #    break

    # Wait for all processes to finish
    if processes > 1:
        pool.close()
        pool.join()

    # Step out of working directory
    os.chdir(cwd)

    # Return
    return


# ==================== #
# Simulate observation #
# ==================== #
def simulate_observation(obs, event_class, models, edisp=True, deadc=0.98,
                         debug=False):
    """
    Simulate events

    Parameters
    ----------
    obs : `~gammalib.GCTAObservation`
        CTA observation
    event_class : string
        Event class
    models : `~gammalib.GModels`
        Models
    edisp : bool, optional
        Enable energy dispersion in simulation
    deadc : float, optional
        Deadtime correction factor
    debug : bool, optional
        Debug simulation
    """
    # Set output filenames
    outevents = 'gps_%s_%s.fits'                   % (event_class, obs.id())
    logfile   = '../../log/gps_ctobssim_%s_%s.log' % (event_class, obs.id())

    # Continue only if output file does not exist
    if not os.path.isfile(outevents+'.gz'):

        # Dump processing
        print('Create "%s" for event class "%s"' % (outevents, event_class))

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

        # Post process event file
        caldb = obs.response().caldb().instrument()
        irf   = obs.response().rspname()
        postprocess_eventfile(outevents, event_class, obs.id(), 'GPS', caldb, irf)

        # Compress FITS file
        os.system('gzip %s' % outevents)

    # ... otherwise notify skipping
    else:
        print('Skip %s (exists already)' % outevents)

    # Return
    return


# ======================= #
# Post process event file #
# ======================= #
def postprocess_eventfile(filename, event_class, id, object, caldb, irf):
    """
    Post process event file

    Parameters
    ----------
    filename : string
        Event file
    event_class : string
        Event class
    id : int
        Observation ID
    object : string
        Observed object
    caldb : string
        Calibration database
    irf : string
        Response string
    """
    # Load FITS file
    fits = gammalib.GFits(filename.rstrip('.gz'))

    # Set OBD_ID keyword
    fits['EVENTS'].card('OBS_ID', id, 'Observation identifier')

    # Set OBJECT keyword
    fits['EVENTS'].card('OBJECT', object, 'Observed object')

    # Set OBSERVER keyword
    fits['EVENTS'].card('OBSERVER', 'CTA Consortium', 'Observer')

    # Set CREATOR keyword
    fits['EVENTS'].card('CREATOR', 'ctobssim (1.6.3)', 'Program which created the file')

    # Set TELLIST keyword
    fits['EVENTS'].card('TELLIST', event_class, 'Telescope IDs')

    # Set CALDB and IRF keywords
    fits['EVENTS'].card('CALDB', caldb, 'Calibration database')
    fits['EVENTS'].card('IRF',   irf,   'Instrument Response Function')

    # Set geographic coordinates
    if 'South' in irf:
        fits['EVENTS'].card('GEOLAT', -24.6272, '[deg] Geographic latitude of array centre')
        fits['EVENTS'].card('GEOLON',  79.4041, '[deg] Geographic longitude of array centre')
        fits['EVENTS'].card('ALTITUDE',   2.15, '[km] Altitude of array centre')
    else:
        fits['EVENTS'].card('GEOLAT', +28.7619, '[deg] Geographic latitude of array centre')
        fits['EVENTS'].card('GEOLON',  17.8900, '[deg] Geographic longitude of array centre')
        fits['EVENTS'].card('ALTITUDE',    2.2, '[km] Altitude of array centre')

    # Get number of events
    nevents = fits['EVENTS'].nrows()

    # Save FITS file
    fits.save(True)

    # Print statistics
    #print('Post processed %s (%d events)' % (filename, nevents))

    # Return
    return


# ============================================ #
# Collect data and write observation container #
# ============================================ #
def collect_data(event_class):
    """
    Collect data and write observation container

    Parameters
    ----------
    event_class : dict
        Event class
    """
    # Dump header
    print('Generate observation definition file for "%s" event class' %
          (event_class['name']))

    # Get current working directory
    cwd = os.getcwd()

    # Set filenames
    inxml      = 'gps_obs_%s_obsdef.xml' % (event_class['name'])
    outxml     = 'gps_obs_%s.xml'        % (event_class['name'])
    reldirname = '../data/%s'            % (event_class['name'])
    absdirname = '$CTADATA/data/%s'      % (event_class['name'])

    # Step into data directory
    os.chdir('gps/obs')

    # Create empty XML file
    xmlout = gammalib.GXml()
    list   = xmlout.append('observation_list title="observation list"')

    # Load XML file
    xml = gammalib.GXml(inxml)

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

        # Set filename
        filename = '%s/gps_%s_%s.fits.gz' % (reldirname, event_class['name'], id)
            
        # Search for a corresponding event file
        files = glob.glob(os.path.expandvars(filename))

        # If we found one file then get the relative filename and append a
        # "EventList" parameter to the XML file
        if len(files) == 1:

            # Get event filename
            eventfile = '%s/%s' % (absdirname, os.path.basename(files[0]))

            # Append XML element
            obs.append('parameter name="EventList" file="%s"' % eventfile)

            # Copy element to output file
            list.append(obs)

        # ... otherwise notify that an event file is missing
        else:
            print('No event file for observation %s' % id)

        # Save XML file
        xmlout.save(outxml)

    # Step out of working directory
    os.chdir(cwd)

    # Return
    return


# ============ #
# Simulate GPS #
# ============ #
def simulate_gps():
    """
    Simulate GPS
    """
    # Create directory tree
    create_directory_tree()

    # Set event types
    #event_classs = [{'name': 'standard',  'suffix': '50h'},
    #                {'name': 'transient', 'suffix': '0.5h'}]
    event_classs = [{'name': 'standard',  'suffix': '50h'}]

    # Loop over event types
    for event_class in event_classs:

        # Create observation definition file
        create_obsdef(event_class)

        # Simulate observations
        simulate_observations(event_class, processes=5)

        # Collect data
        collect_data(event_class)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Dump header
    print('****************')
    print('* Simulate GPS *')
    print('****************')

    # Simulate GPS
    simulate_gps()

    # We are done
    print('... done.')
