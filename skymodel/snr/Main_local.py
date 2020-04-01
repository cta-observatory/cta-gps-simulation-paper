#---------
# Main SNR TEHANU new version 2018.
#--------

import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as mticker
import matplotlib.ticker as ticker
from itertools import zip_longest
import os
import csv
from random import choices
import random
from matplotlib import rc
import sys
import os
from SNR import *  # SNR.py is where all functions are


#--------------------------------------------------------------------------#
#           takes one file (.txt), one number (int), and returns one file (.txt)
#           Takes parameters.txt  - > creates output.txt
#--------------------------------------------------------------------------#

#fixing seed of random generator:
# change
number_seed=1
random.seed(number_seed)





#input
file_parameters=sys.argv[1]
NUMBER=np.int(sys.argv[2])  # number of MC realizations



index_of_parameter=get_first_nbr_from_str(file_parameters)
data=np.loadtxt(file_parameters,usecols=1)
alpha=data[0]
Kep=data[1]
D=data[2]
eta=data[3]
KNEE=data[4]

# output:
FOLDER_FOR_RESULTS='ALL_FILES_'+str(index_of_parameter)
if not os.path.exists(FOLDER_FOR_RESULTS):
    os.makedirs(FOLDER_FOR_RESULTS)


#test for one SNR:


## Parameters:
#a=2.1   # slope accelerated particles, in energy
#Kep=10**-2
#eta=0.1  #CR efficiency at the shock
#n0=1. # cm**-3
#B0=4.3 # muG
## Photon field directly in SNR.py
#
#ENERGY_GRID=np.logspace(-1.5,2.3,40)* u.TeV
#TIME_GRID=np.logspace(-1,np.log10(1),40) # kyear
#
#
#def save_one_SNR_to_file(SNR,file):
#    with open(file, 'w') as text_file:
#        writer = csv.writer(text_file, delimiter='\t')
#        writer.writerow(["Time[kyr]","Energy[TeV]", \
#                 "Diff_spectrum[TeV/s]"])
#
#
#        for i in range (0,len(SNR.TIME)):
#            ENERGY=np.array(SNR.ENERGY)
#            for j in range (0,len(SNR.ENERGY)):
#                writer.writerow((SNR.TIME[i], ENERGY[j] , SNR.LGAMMA_DIFF_T[i][j]))
#
#
#
#print( 'Type Ia .... ')
#SNR1=SNR()
#SNR1.type=1
#SNR1.n0=n0
#SNR1.B0=B0
#SNR1.set_density()
#SNR1.TIME=TIME_GRID
#SNR1.ENERGY=ENERGY_GRID
#SNR1.typical_associated_parameters()
#SNR1.alpha=a
#print(' SNR1.alpha = ', SNR1.alpha, 'SNR1.n0 =' , SNR1.n0 , 'SNR1.rho0 =' , SNR1.rho0 )
#SNR1.Kep=Kep
#SNR1.dist=1.
#SNR1.calculate_diff_spectrum_TIME()
#save_one_SNR_to_file(SNR1,FOLDER_FOR_RESULTS+'/'+'TypeI.txt')
#
#
#
#
#
#
#exit()
#


#calculations + writing to files:
for i in range (0,NUMBER):
    print( 'Sim = ', i)
    LIST_SNR=one_realization_only_pevatrons(alpha,Kep, D, eta, KNEE)
    file_output=FOLDER_FOR_RESULTS+'/'+'results_'+str(i)+'.txt'
    print( 'LIST_SNR[0].TIME =', len(LIST_SNR[0].TIME)-1)
    save_one_LIST_to_file_cyril(LIST_SNR,file_output)


exit()
