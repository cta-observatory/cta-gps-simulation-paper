from evo_gelfand_model import *

age = 3000. #age of system
dens = 0.01 #density of ISM
tc = 5305. #Pulsar characteristic age
nn = 3. #braking index
t0 = ((2*tc)/(nn-1))-age #Pulsar initial spin-down timescale
e0 = 1.0e51 #Energy of the SN
mej2 = ((9.*u.M_sun).to('g')).value #SN ejected mass
eta = 0.0313 #Magnetization of the PWN
eps = 0.25 #Confinment factor
l0 = 2.29e38 #Pulsar initial spin-down
ebreak = 0.045*gp.TeV_to_erg #energy break of the injection spectrum (Broken PWL)
alpha1 = 1.10 #first index of the injection spectrum (Broken PWL)
alpha2 = 2.52 #second index of the injection spectrum (Broken PWL)
Tfir = 30 #Temperature of FIR field in K
Ufir = 3.8*gp.eV_to_erg #density of FIR field in erg/cm^3
Tnir = 3000 #Temperature of NIR field in K
Unir = 25.*gp.eV_to_erg #density of NIR field in erg/cm^3
dist=13300. #distance of the system

launch_evo(age, l0, e0, mej2, nn, eta, t0, eps, dens, Tfir, Ufir, Tnir, Unir, 250, 200.,  ebreak, alpha1, alpha2, dist)
