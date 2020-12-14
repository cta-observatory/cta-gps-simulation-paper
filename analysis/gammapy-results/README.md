# CTA-GPS catalogue production with gammapy

## Content

- draft of document containing the analysis description and results, "CTA-GPS_gammapy.pdf"
- main-catalogue, "SeedList_fitted_step1.fits", also available in .yaml format
- catalogue subset containing further informations for Pevatron studies, "SeedList_fitted_step1_selection_criterion_pevatron.fits"

## Acces to reduced datasets:
The global datasets and models produced by the gammapy analysis can be found [here](https://1drv.ms/f/s!AiTtm00zHBzSiJZBumhtQEqOETxRkg)
(due to the file size limitation, it cannot be share on GitHub).

How to use :
- Download and unzip the folder "CTA-GPS_read_v0-18-2.zip"
- Install gammapy version v0.18.2 using the environment provided in the [documentation](https://docs.gammapy.org/0.18.2/install/index.html)
- Read the datasets using the script CTA-GPS_read_v0-18-2.py

This script includes several functions that "patch" gammapy to ensure the compatibility of the datasets and models files that were generated with a different version and custom code.
I also added a function that extract a sub-region centred on a given model from the global dataset (as an example I extracted a sub-datasets containing Vela-Junior and the GC region).
So from this you can re-analyse any sub-region surrounding your sources of interest. However the spatial and spectral binning will be the same than the ones of the global GPS analysis.
If you want to change the analysis binning you have to build your own datasets from the observations files.
 
## Other informations

### Contact
Questions or suggestions are welcome, you can send them to me on slack in the CTA or gammapy channels, or by mail at quentin.remy@mpi-hd.mpg.de

### About the main catalogue production

See "CTA-GPS_gammapy.pdf".
Further details on the analysis, description of the catalogues column content, and more results figures will be added to this document progressively.

### About the selection for the Pevatron study

I made a minimal selection so the list of sources is shorter than main catalogue.
First I removed the sources with TSnull_postTeV<25 (not significant in the 1-200 TeV range).
Then I estimated the excess counts for each sources in the 10-200 / 50-200 TeV / 100-200 TeV ranges
(no simulated data above 200 TeV), and selected only the ones with at least one count above 50 TeV.
Finally for the remaining sources an exponential cut-off model have been fitted instead of the baseline log-parabola model.
The table contains:
- sources names
- longitudes
- latitudes
(I did not copy the others spatial parameters that appear in the main catalogue)
- the exponential cut-off model models parameters (amplitude, index, Ecut, Eref)
- the number of excess events in energies bands (Npred_...) and the integrated photon flux associated (F_...)
- Test statistic between the  log-parabola  and exponential cut-off models (TS_lp-ecpl > 0 if ecpl is prefered)
 
### Units
For now units are included in the .yaml files but missing in the fits file header.
- amplitude: '1 / (cm2 s TeV)'
- Integrated flux : '1 / (cm2 s)'
- Ecut and Eref : 'TeV'