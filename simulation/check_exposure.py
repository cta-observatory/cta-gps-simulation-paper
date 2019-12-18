#! /usr/bin/env python
# ==========================================================================
# Generate exposures for Galactic Plane survey
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
import gammalib
import ctools
import cscripts


# ======================= #
# Enter working directory #
# ======================= #
def enter_wd(create=False):
    """
    Enter working directory

    Parameters
    ----------
    create : bool, optional
        Create directory tree

    Returns
    -------
    cwd : string
        Current working directory
    """
    # Get current working directory
    cwd = os.getcwd()

    # Create directory tree
    if create:
        try:
            os.makedirs('gps')
        except:
            pass
        try:
            os.makedirs('gps/exposure')
        except:
            pass

    # Step into working directory
    os.chdir('gps/exposure')

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


# ================================= #
# Create filename for time interval #
# ================================= #
def get_filename(start, stop):
    """
    Create filename for time interval

    Parameters
    ----------
    start : integer
        Start year
    stop : integer
        Stop year

    Returns
    -------
    filename : string
        Filename
    """
    # Set filename
    filename = '_%4.4d-%4.4d' % (start, stop)

    # Return filename
    return filename


# ================================== #
# Create observation definition file #
# ================================== #
def create_obsdef(pntfile, start, stop, rad=5.0):
    """
    Create observation definition file for binary

    Parameters
    ----------
    pntfile : string
        Pointing file
    start : integer
        Start year
    stop : integer
        Stop year
    rad : float, optional
        Simulation radius (deg)
    """
    # Set start and stop times
    tmin = '%4.4d-01-01T00:00:00' % (start)
    tmax = '%4.4d-12-31T23:59:59' % (stop)

    # Set observation definition filename
    outobs  = 'gps_obs' + get_filename(start, stop) + '.xml'
    logfile = 'gps_obs' + get_filename(start, stop) + '.log'

    # Enter working directory
    cwd = enter_wd(True)

    # Run make_pointings.py if pointings file does not exist
    if not os.path.isfile(pntfile):
        print('Create GPS pointings')
        os.system(os.path.expandvars('${CTAGITROOT}/simulation/gps/make_pointings.py'))

    # Run csobsdef if observation definition file does not exist
    if not os.path.isfile('gps_obs.xml'):
        print('Create observation definition file')
        obsdef = cscripts.csobsdef()
        obsdef['inpnt']   = pntfile
        obsdef['outobs']  = 'gps_obs.xml'
        obsdef['rad']     = rad
        obsdef['logfile'] = 'gps_csobsdef.log'
        obsdef.logFileOpen()
        obsdef.execute()

    # Run csobsselect if observation definition file does not exist
    if not os.path.isfile(outobs):
        print('Select pointings from observation definition file for years %d-%d' % (start, stop))
        obsselect = cscripts.csobsselect()
        obsselect['inobs']     = 'gps_obs.xml'
        obsselect['outobs']    = outobs
        obsselect['pntselect'] = 'CIRCLE'
        obsselect['coordsys']  = 'GAL'
        obsselect['glon']      = 0.0
        obsselect['glat']      = 0.0
        obsselect['rad']       = 180.0
        obsselect['tmin']      = tmin
        obsselect['tmax']      = tmax
        obsselect['logfile']   = logfile
        obsselect.logFileOpen()
        obsselect.execute()

    # Exit working directory
    exit_wd(cwd)

    # Return
    return


# ==================== #
# Create exposure cube #
# ==================== #
def create_expcube(start, stop):
    """
    Create exposure cube

    Parameters
    ----------
    start : integer
        Start year
    stop : integer
        Stop year
    """
    # Set observation definition filename
    inobs   = 'gps_obs'     + get_filename(start, stop) + '.xml'
    expcube = 'gps_expcube' + get_filename(start, stop) + '.fits'
    logfile = 'gps_expcube' + get_filename(start, stop) + '.log'

    # Enter working directory
    cwd = enter_wd()

    # Continue only if exposure cube does not exist
    if not os.path.isfile(expcube):

        # Dump header
        print('Create exposure cube for years %d-%d' % (start, stop))

        # Setup task parameters
        ctexpcube = ctools.ctexpcube()
        ctexpcube['inobs']   = inobs
        ctexpcube['incube']  = 'NONE'
        ctexpcube['ebinalg']  = 'LOG'
        ctexpcube['emin']     =   0.1
        ctexpcube['emax']     = 100.0
        ctexpcube['enumbins'] = 3
        ctexpcube['coordsys'] = 'GAL'
        ctexpcube['proj']     = 'CAR'
        ctexpcube['xref']     = 0.0
        ctexpcube['yref']     = 0.0
        ctexpcube['nxpix']    = 3600
        ctexpcube['nypix']    = 200
        ctexpcube['binsz']    = 0.1
        ctexpcube['outcube'] = expcube
        ctexpcube['logfile'] = logfile
        ctexpcube.logFileOpen()

        # Generate exposure cube
        ctexpcube.execute()

    # Exit working directory
    exit_wd(cwd)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Dump header
    print('************************************')
    print('* Check exposure of GPS simulation *')
    print('************************************')

    # Create observation definition file
    create_obsdef('GPS_pointings.dat', 2021, 2021)
    create_obsdef('GPS_pointings.dat', 2022, 2022)
    create_obsdef('GPS_pointings.dat', 2023, 2023)
    create_obsdef('GPS_pointings.dat', 2024, 2024)
    create_obsdef('GPS_pointings.dat', 2025, 2025)
    create_obsdef('GPS_pointings.dat', 2026, 2026)
    create_obsdef('GPS_pointings.dat', 2027, 2027)
    create_obsdef('GPS_pointings.dat', 2028, 2028)
    create_obsdef('GPS_pointings.dat', 2029, 2029)
    create_obsdef('GPS_pointings.dat', 2030, 2030)
    #
    create_obsdef('GPS_pointings.dat', 2021, 2022) # STP
    create_obsdef('GPS_pointings.dat', 2023, 2030) # LTP
    #
    create_obsdef('GPS_pointings.dat', 2021, 2023)
    create_obsdef('GPS_pointings.dat', 2021, 2024)
    create_obsdef('GPS_pointings.dat', 2021, 2025)
    create_obsdef('GPS_pointings.dat', 2021, 2026)
    create_obsdef('GPS_pointings.dat', 2021, 2027)
    create_obsdef('GPS_pointings.dat', 2021, 2028)
    create_obsdef('GPS_pointings.dat', 2021, 2029)
    create_obsdef('GPS_pointings.dat', 2021, 2030)

    # Create exposure cube
    create_expcube(2021, 2021)
    create_expcube(2022, 2022)
    create_expcube(2023, 2023)
    create_expcube(2024, 2024)
    create_expcube(2025, 2025)
    create_expcube(2026, 2026)
    create_expcube(2027, 2027)
    create_expcube(2028, 2028)
    create_expcube(2029, 2029)
    create_expcube(2030, 2030)
    #
    create_expcube(2021, 2022) # STP
    create_expcube(2023, 2030) # LTP
    #
    create_expcube(2021, 2023)
    create_expcube(2021, 2024)
    create_expcube(2021, 2025)
    create_expcube(2021, 2026)
    create_expcube(2021, 2027)
    create_expcube(2021, 2028)
    create_expcube(2021, 2029)
    create_expcube(2021, 2030)
    
    # We are done
    print('... done.')
