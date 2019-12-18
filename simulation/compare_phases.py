#! /usr/bin/env python
# ==========================================================================
# Compare phase histograms
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


# =========== #
# Plot phases #
# =========== #
def plot_phases(name, constname, phasename, lcname, nphases=40):
    """
    Plot phases

    Parameters
    ----------
    constname : str
        Name of phase curve FITS file for constant source
    phasename : str
        Name of phase curve FITS file for phased source
    """
    # Read observation containers
    obs_const = gammalib.GObservations(constname)
    obs_phase = gammalib.GObservations(phasename)

    # Read phase curve FITS file
    fits = gammalib.GFits(lcname)

    # Extract phase curve
    table    = fits.table(1)
    c_phase  = table['PHASE']
    c_norm   = table['NORM']
    m_phases = [0.0 for i in range(2*table.nrows())]
    m_norm   = [0.0 for i in range(2*table.nrows())]
    for row in range(table.nrows()):
        m_phases[row]               = c_phase.real(row)
        m_phases[row+table.nrows()] = c_phase.real(row)+1.0
        m_norm[row]                 = c_norm.real(row)
        m_norm[row+table.nrows()]   = c_norm.real(row)

    # Initialise phases
    nbins        = 2*nphases
    dphases      = 1.0/float(nphases)
    phases       = [(i+0.5)*dphases for i in range(nbins)]
    phases_const = [0.0 for i in range(nbins)]
    phases_phase = [0.0 for i in range(nbins)]
    source_const = [0.0 for i in range(nbins)]
    source_phase = [0.0 for i in range(nbins)]

    # Loop over all observations
    for run in obs_const:

        # Get events
        events = run.events()

        # Loop over all events
        for event in events:

            # Get phase and MC ID
            phase = event.phase()
            mc_id = event.mc_id()

            # Collect phase
            i1 = int(phase/dphases)
            i2 = int((phase+1.0)/dphases)
            if i1 >=0 and i1 < nbins:
                phases_const[i1] += 1.0
            if i2 >=0 and i2 < nbins:
                phases_const[i2] += 1.0

            # Collect phases for source events
            if mc_id == 1:
                if i1 >=0 and i1 < nbins:
                    source_const[i1] += 1.0
                if i2 >=0 and i2 < nbins:
                    source_const[i2] += 1.0

    # Loop over all observations
    for run in obs_phase:

        # Get events
        events = run.events()

        # Loop over all events
        for event in events:

            # Get phase and MC ID
            phase = event.phase()
            mc_id = event.mc_id()

            # Collect phase
            i1 = int(phase/dphases)
            i2 = int((phase+1.0)/dphases)
            if i1 >=0 and i1 < nbins:
                phases_phase[i1] += 1.0
            if i2 >=0 and i2 < nbins:
                phases_phase[i2] += 1.0

            # Collect phases for source events
            if mc_id == 1:
                if i1 >=0 and i1 < nbins:
                    source_phase[i1] += 1.0
                if i2 >=0 and i2 < nbins:
                    source_phase[i2] += 1.0

    # Derive phase curve
    x      = []
    ratios = []
    errors = []
    for i in range(nbins):
        if source_const[i] > 0.0:
            ratio = source_phase[i]/source_const[i]
            error = math.sqrt(phases_phase[i])/source_const[i]
            x.append(phases[i])
            ratios.append(ratio)
            errors.append(error)

    # Plot phase
    plt.errorbar(x, ratios, yerr=errors, fmt='ro')
    plt.plot(m_phases, m_norm, 'b-')
    plt.xlabel('Phase')
    plt.ylabel('Event ratio')
    plt.title(name)

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
    options = [{'option': '-n', 'value': '10'},
               {'option': '-s', 'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    nphases = int(options[0]['value'])
    name    = options[1]['value']

    # Get phase curve file name
    models   = gammalib.GModels(file_model)
    model    = models[0]
    filename = model.temporal().filename()
    title    = model.name()

    # Create figure
    plt.figure()

    # Plot phases
    plot_phases(name, file_const, file_phase, filename, nphases=nphases)

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
