#! /usr/bin/env python
# ==========================================================================
# Display phase information from event list
#
# Copyright (C) 2017 Juergen Knoedlseder
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
import gammalib
import cscripts
try:
    import matplotlib.pyplot as plt
    plt.figure()
    plt.close()
except (ImportError, RuntimeError):
    print('This script needs the "matplotlib" module')
    sys.exit()


# =========== #
# Plot phases #
# =========== #
def plot_phases(filename, nphases=40, color='r'):
    """
    Plot phases

    Parameters
    ----------
    filename : str
        Name of lightcurve FITS file
    """
    # Read observation container. If an exception occurs then try loading the
    # file as an event list and build an observation container.
    try:
        obs = gammalib.GObservations(filename)
    except:
        obs = gammalib.GObservations()
        obs.append(gammalib.GCTAObservation(filename))

    # Initialise phases
    phases = []

    # Loop over all observations
    for run in obs:

        # Get events
        events = run.events()

        # Loop over all events
        for event in events:

            # Get phase
            phase = event.phase()

            # Collect phase
            phases.append(phase)
            phases.append(phase+1.0)

    # Plot phase
    plt.hist(phases,bins=nphases,range=(0.0,2.0),facecolor=color)
    plt.xlabel('Phase')
    plt.ylabel('Events')

    # Return
    return


# =========== #
# Show phases #
# =========== #
def show_phases():
    """
    Show phases
    """
    # Set filenames
    file_const = 'gps_obs_selected_phased_const.xml'
    file_phase = 'gps_obs_selected_phased.xml'
    file_model = 'gps_models.xml'

    # Set usage string
    usage = 'show_phases.py [-n nbins]'

    # Set default options
    options = [{'option': '-n', 'value': '10'}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    nphases = int(options[0]['value'])*2

    # Get title
    models = gammalib.GModels(file_model)
    model  = models[0]
    title  = model.name()

    # Create figure
    plt.figure()

    # Plot phases
    plot_phases(file_const, nphases=nphases, color='r')
    plot_phases(file_phase, nphases=nphases, color='g')

    # Plot title
    plt.title(title)

    # Show figure
    plt.show()

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show phases
    show_phases()
