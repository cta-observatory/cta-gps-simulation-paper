# Synthetic binaries

Based on Dubus et al. 2017 (DOI:10.1051/0004-6361/201731084)

Description by Guillaume

(1) Binary parameters

I?ve assumed the Faucher parametrisation of spirals with a radial distribution corresponding to OB stars, as implemented in the python codes pointed out by Andrea. This is slightly different from I assumed in the paper in that it the distribution is more peaked within -60<l<60 degrees.
The other choices of binary parameters are as in my paper. Note that galactic latitude = 0 degrees for all of them (you can pick a random number between +-2 degrees).

(2) Lightcurves

My orbital lightcurves are phase and flux in ph/cm2/s > 1 TeV.

I also give in mCrab with 1 Crab = 1.823e-11 ph/cm2/s > 1 TeV 

I do not have spectral info and I think this would be overkill for the present purposes (quite complicated to do, not necessarily more realistic = not worth the effort in my opinion). I would assume a powerlaw spectrum with a photon index of 2.5 regardless of orbital phase.

(3) Population

I can generate a population of 200 binaries (a reasonable number of the overall population expected in our Galaxy) and to send you the resulting 200 files. Several will have peak flux <0.1 mCrab but you can just remove them as necessary (or I can do that for you).