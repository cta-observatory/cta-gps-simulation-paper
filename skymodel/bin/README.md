Binaries

___LSI +61 303:___

data for the maximal and quiescent flux are taken from
https://pos.sissa.it/301/712/pdf

No significant spectral variability is observed. On the top of the orbital variability observed by Veritas the peak flux is supposed to have a shift in the super orbital time scale, phasecurve_LSI.fits. 

___PSR B1259-63___

HESS data are taken from https://arxiv.org/pdf/1708.00895.pdf

No significant spectral variability is observed. Orbital modulation is in the file phasecurve_PSRB1259.fits.



____1FGL J1018.6?5856  === HESS J1018-589____


Data from http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1503.02711 

No spectral variability is known. Orbital modulation is in the file phasecurve_HESSJ1018_DC.fits (there was no need to change the file used for DC1 ).


___LMC P3___

HESS observations are summarised in 
http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1801.06322 

on-peak and off-peak phases seems to have different slopes, thus there are two inputs for this source in xml file (?LMC P3? off and ?LMC P3 on?). Orbital modulations for the on and off phase are in phasecurve_LMCP3_on.fits and phasecurve_LMCP3_off.fits correspondingly. For the off-phase fluxes were chosen to be consistent with the observed orbital-average flux.

___N1___

New fake binary in the crowded field (Galactic Centre). The selected flux is of the order of the flux of recently discovered in the GC region HESS J1746?285,
http://adsabs.harvard.edu/abs/2017arXiv170604535H, as source with a higher fluxes would be probably already known.
 

___ LS 5039___

LS 5039 was one of the first gamma-ray binaries to be discovered [1]. The source displays spectral variability as a unction of the orbital phase [2], related to geometry given by the position of the compact object, the companion star providing the target photon field, and the line of sight. Two broad intervals, SUPC and INFC, have been considered in the respective LC .fits files, which are respectivelly best-fit with a PL and an ECPL. Also, two separate entries in the xml file have been implemented.

[1]Aharonian et al. 2005: http://adsabs.harvard.edu/abs/2005Sci...309..746A
[2]Aharonian et al. 2006: http://adsabs.harvard.edu/abs/2006A%26A...460..743A


___PSR J2032+4127___

This is the last gamma-ray binari discovered so far. The system lies on top of TeV 2032+4130, an extended VHE source. PSR J2032 displays a period of approximately 50 years, with the most recent periastron having occurred on November 2017. Therefore, CTA may not be in a timely situation to observe it during periastron. The source is included in the binary sky model nevertheless for completeness. 

VERITAS and MAGIC reported on the detection of the source [1], retrieving also two distinct spectral behaviours when approaching the periastron phase and right at periastron. Since the system has a very large orbital period, however, we considered only the emission in a single orbital phase range (with a 0.1 width) with an averaged flux and assuming an ECPL model, as suggested in [1].

[1]:VERITAS & MAGIC Collaborations 2018: http://adsabs.harvard.edu/abs/2018arXiv181005271T


___HESS J1832-093___

HESS J1832 is a gamma-ray binary candidate discovered with HESS [1]. Despite that there is so far no indication of variability/periodicity at VHEs, the source was followed up at X-rays, yielding a statistically significant variability at this energy range [2] (see also [3]). Together with its VHE point-like nature, this suggested the possibility of the source to be a new gamma-ray binary system. 

[1] HESS Collaboration 2014: http://adsabs.harvard.edu/abs/2014arXiv1411.0572H
[2] Eger et al. 2016: http://adsabs.harvard.edu/abs/2016MNRAS.457.1753E
[3] Mori et al. 20017: http://adsabs.harvard.edu/abs/2017ApJ...848...80M


___HESS J0632+057___

HESS J0632+057 was discovered with HESS [1], and displays correlated variability in X-rays and TeVs. The source displays several emission states along its 315days orbital period, consisting of a flaring phase, a dip, a second maximum, and a baseline. These states do not display however any significant spectral variability. Therefore, a single phasecurve_HESSJ0632.fits file has been produced. 

[1] Aharonian et al. 2007: http://adsabs.harvard.edu/abs/2007A%26A...469L...1A
[2] HESS and VERITAS Collaborations 2014: http://adsabs.harvard.edu/abs/2014ApJ...780..168A

___N2___

Fake binary system located in a non-crowded region and with a moderately low flux level. This test source is considered to probe the capabilities of CTA for faint sources. A system similar to LMC P3 has been considered, but 10 times fainter, with a ~10d orbital period and gamma-rays being mainly produced in only two orbital phase bins. Contrary to the current behaviour observed in LMC P3, though, we considered also a "baseline" emission, of about 1/3 of the flux proposed for the maximum. 
