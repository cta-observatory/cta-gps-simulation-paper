# Sky model

The final model format follows the format required to simulate using ctools

http://cta.irap.omp.eu/ctools/users/user_manual/models.html

For example models see the CTA DC-1 products. You can find instructions to get them at

https://forge.in2p3.fr/projects/data-challenge-1-dc-1/wiki/Getting_data

## Note on the PWN/SNR sky model
For the PWN and SNR files. The synthetic galaxy simulation files that were used are the [`PWNe_final_population.txt`](https://github.com/cta-observatory/cta-gps-simulation-paper/tree/master/skymodel/pwn) and [`results_0.txt`](https://github.com/cta-observatory/cta-gps-simulation-paper/tree/master/skymodel/snr/ALL_FILES_0) files.  
The position of each PWN was initiated from the position of a core-collapse SN. ID number of the PWN and SNR are consistent in the above txt files. (X,Y) positions match for the SNRs and the PWNe.  
**WARNING:** Only in the XML files, the ID number of SNR and PWN are offset by 1 (e.g. snr_86 = pwn85).  

Note also that the transformation (X,Y) ==> (GLON, GLAT) is computed using a [gammapy](https://docs.gammapy.org/0.18.2/api/gammapy.utils.coordinates.galactic.html#gammapy.utils.coordinates.galactic) routine. The routine assumes Sun at Y=-8.5 kpc but at function call, axes were swapped to assume the Sun at the position X,Y = (8.5,0) kpc ([astropy](https://docs.astropy.org/en/stable/api/astropy.coordinates.Galactocentric.html) assumes X=-8.1 kpc). 
