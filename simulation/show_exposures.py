#! /usr/bin/env python
# ==========================================================================
# Show exposure maps
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
import ctools
import cscripts
import matplotlib.pyplot as plt


# ======== #
# Plot map #
# ======== #
def plot_exposure(inmap, sub, imap, smooth=0.0, labeled=True):
    """
    Plot exposure map

    Parameters
    ----------
    inmap : `~gammalib.GSkyMap()`
        Input sky map
    sub : pyplot
        Frame for map
    imap : int
        Index of map to plot
    smooth : float, optional
        Map smoothing parameter (degrees)
    labeled : bool, optional
        Add map labels
    """
    # Optionally smooth map
    if smooth > 0.0:
        map = inmap.copy()
        map.smooth('GAUSSIAN', smooth)
    else:
        map = inmap
    
    # Create array from skymap
    array = []
    v_max = 0.0
    for iy in range(map.ny()):
        row = []
        for ix in range(map.nx()):
            index = ix+iy*map.nx()
            value = map[index,imap]
            row.append(value)
        array.append(row)

    # Show Aitoff projection
    c    = sub.imshow(array, extent=(180.0,-180.0,-10.0,10.0), cmap='jet')
    if labeled:
        cbar = plt.colorbar(c, orientation='vertical', shrink=0.5)
        cbar.ax.set_ylabel('exp. (cm$^2$ sec)', fontsize=10)
        sub.set_xlabel('Galactic longitude (deg)', fontsize=10)
        sub.set_ylabel('Gal. latitude (deg)', fontsize=10)

    # Return
    return


# ========================= #
# Show exposures for period #
# ========================= #
def show_exposure_period(start, stop):
    """
    Show exposures for period

    Parameters
    ----------
    start : integer
        Start year
    stop : integer
        Steop year
    """
    # Create figure
    fig = plt.figure('Exposure', (12, 7))
    fig.subplots_adjust(left=0.07, right=1.05, top=0.98, bottom=0.05)

    # Build file name and title
    filename = 'gps_expcube_%4.4d-%4.4d.fits' % (start, stop)
    if start == stop:
        title = 'Exposure for year %4.4d' % (start)
    elif start == 2021 and stop == 2022:
        title = 'Exposure for STP (%4.4d-%4.4d)' % (start, stop)
    elif start == 2023 and stop == 2030:
        title = 'Exposure for LTP (%4.4d-%4.4d)' % (start, stop)
    else:
        title = 'Exposure for years %4.4d-%4.4d' % (start, stop)

    # Add title
    fig.suptitle(title, fontsize=16)

    # Load exposure
    exposure = gammalib.GSkyMap(filename)

    # Plot
    ax1 = fig.add_subplot(411)
    plot_exposure(exposure, ax1, 0) # 100 GeV
    ax1.text(-175.0,5.0,'100 GeV', color='w',
             horizontalalignment='right', verticalalignment='center')
    #
    ax2 = fig.add_subplot(412)
    plot_exposure(exposure, ax2, 1) # 1 TeV
    ax2.text(-175.0,5.0,'1 TeV', color='w',
             horizontalalignment='right', verticalalignment='center')
    #
    ax3 = fig.add_subplot(413)
    plot_exposure(exposure, ax3, 2) # 10 TeV
    ax3.text(-175.0,5.0,'10 TeV', color='w',
             horizontalalignment='right', verticalalignment='center')
    #
    ax4 = fig.add_subplot(414)
    plot_exposure(exposure, ax4, 3) # 100 TeV
    ax4.text(-175.0,5.0,'100 TeV', color='w',
             horizontalalignment='right', verticalalignment='center')

    # Show figure
    plt.show()

    # Return
    return


# ====================== #
# Show exposure per year #
# ====================== #
def show_yearly_exposures():
    """
    Show exposure per year
    """
    # Create figure
    fig = plt.figure('Exposure', (12, 7))
    fig.subplots_adjust(left=0.05, right=1.0, top=0.94, bottom=0.05)

    # Set map index and energy
    imap  = 1
    title = 'Yearly exposures at 1 TeV'

    # Add title
    fig.suptitle(title, fontsize=16)

    # Loop over years
    for index, year in enumerate(range(2021,2031)):

        # Build file name and title
        filename = 'gps_expcube_%4.4d-%4.4d.fits' % (year, year)

        # Load exposure
        exposure = gammalib.GSkyMap(filename)

        # Plot
        ax = fig.add_subplot(10,1,index+1)
        plot_exposure(exposure, ax, imap, labeled=False)
        ax.text(-175.0,5.0,'%4.4d' % year, color='w',
                horizontalalignment='right', verticalalignment='center')

    # Show figure
    plt.show()

    # Return
    return


# =================================== #
# Show cumulative exposure after year #
# =================================== #
def show_cumulative_exposure():
    """
    Show cumulative exposure after year
    """
    # Create figure
    fig = plt.figure('Cumulative exposure', (12, 7))
    fig.subplots_adjust(left=0.05, right=1.0, top=0.94, bottom=0.05)

    # Set map index and energy
    imap  = 1
    title = 'Cumulative exposure at 1 TeV'

    # Add title
    fig.suptitle(title, fontsize=16)

    # Loop over years
    for index, year in enumerate(range(2021,2031)):

        # Build file name and title
        filename = 'gps_expcube_%4.4d-%4.4d.fits' % (2021, year)

        # Load exposure
        exposure = gammalib.GSkyMap(filename)

        # Set text
        if year == 2021:
            text = '%4.4d' % year
        else:
            text = '2021-%4.4d' % year

        # Plot
        ax = fig.add_subplot(10,1,index+1)
        plot_exposure(exposure, ax, imap, labeled=False)
        ax.text(-175.0, 5.0, text, color='w',
                horizontalalignment='right', verticalalignment='center')

    # Show figure
    plt.show()

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show exposure for period
    show_exposure_period(2021, 2022)
    show_exposure_period(2023, 2030)

    # Show exposure per year
    show_yearly_exposures()

    # Show cumulative exposure after year
    show_cumulative_exposure()

