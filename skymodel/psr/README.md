This directory contains xml models of 38 pulsars prepared for the GPS simulations.
Here are some information on the produced models:

1. Pulsars are sub-divided in two classes: 
   - Vela-like pulsars (J0007+7303, J0633+1746, J0835-4510, J1057-5226, J1413-6205, J1418-6058, J1709-4429, J1732-3131, J1813-1246, J1846+0919, J1952+3252, J2021+3651, J2021+4026) 
   - Crab-like pulsars (all others)

2. The spectra of Crab-like pulsars are taken from 3FHL and correspond to the power-law function.
For the Crab pulsar we adopted the power-law spectrum of the peak emission (P1+P2) calculated for the joint Fermi-LAT and MAGIC spectral data as reported in Ansoldi+2016.

3. The spectra of the Vela-like pulsars contains HE and VHE components.
HE components of the Vela-like pulsars corresponds to the
   - LogParaola function from 3FHL for J0007+7303.
   - Power-law function with exponential cutoff and free beta-parameter (PLEC) from Ahen+2016 for J0633+1746 (Geminga). Phase-averaged spectrum.
   - PLEC from HESS Collaboration+2018 paper for J0835-4510 (Vela). Phase-averaged spectrum.
   - Power-law function with exponential cutoff and the beta-parameter fixed to 1 (PLEC1) from 2PC for J1057-5226.
   - PLEC1 from 2PC for J1413-6205.
   - PLEC1 from 2PC for J1418-6058.
   - PLEC from 2PC for J1709-4429.
   - PLEC1 from 2PC for J1732-3131.
   - PLEC1 from 2PC for J1813-1246.
   - PLEC1 from 2PC for J1846+0919.
   - PLEC from 2PC for J1952+3252
   - PLEC from 2PC for J2021+3651
   - PLEC from 2PC for J2021+4026

4. 64-bin pulse profiles are provided by David A. Smith and collaborators. They are are normalised to 1 (maximum=1), and converted to the fits format. Noisy profiles were re-binned using 32 bins for for J1119-6127, and 16 bins for J1459-6053 and J1813-1246

5. Timing properties (MJD, F0, F1, F2) are taken from Kerr+2015 (https://confluence.slac.stanford.edu/display/GLAMCOG/LAT+Gamma-ray+Pulsar+Timing+Models). These timing solutions most probably refer to the same reference date as those used by David A. Smith because the resulting pulse profiles are not shifted in phase.
5 pulsars are missing in Kerr+2015: J1119-6127, J1514-4946, J1648-4611, J1838-0537 and J2215+5135. For them I used 2PC timing parameters. However, there is an off-set between 2PC and David's profiles of these pulsars.

6. Coordinates are taken from Kerr+2015/2PC ephemerides.
