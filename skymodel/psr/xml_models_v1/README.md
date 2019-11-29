This directory contains xml models of pulsars for the GPS simulations.
This is a preliminary version, which contains pulsars with a Crab-like spectra.
The Vela-like components are still to be provided by Arache.

Here are some additional comments:
1. For the moment all pulsars has a power-law spectrum from 3FHL.
2. Even for those pulsars which are better fitted with LogParabola, I used the power-law spectral indexes, reported in 3FHL catalog as well.
3. Coordinates are taken from 3FHL.
4. Since Crab pulsar is not included in 3FHL, its spectral parameters are adopted from joint Fermi-LAT and MAGIC fit (Ansoldi+2016). Coordinates of the Crab pulsar are taken as those reported for the Crab nebula in 3FHL.
5. Pulse profiles provided by David are normalised to 1 (maximum=1), and converted to the fits format.
6. Timing properties (MJD, F0, F1, F2) are taken from Kerr+2015 (https://confluence.slac.stanford.edu/display/GLAMCOG/LAT+Gamma-ray+Pulsar+Timing+Models). 
These timing solutions most probably refer to the same reference date as those used by David, because the resulting pulse profiles are not shifted in phase.
7. 5 pulsars are missing in Kerr+2015: J1119-6127, J1514-4946, J1648-4611, J1838-0537 and J2215+5135. For them I used 2PC timing parameters. However, there is an offset between 2PC and David's profiles of these pulsars.
