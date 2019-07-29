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
from matplotlib import rc

import sys
import os
from SNR import *  # SNR.py is where all functions are


#--------------------------------------------------------------------------#
#           takes one file (.txt), one number (int), and returns one file (.txt)
#           Takes parameters.txt  - > creates output.txt
#--------------------------------------------------------------------------#

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

#calculations + writing to files:
for i in range (0,NUMBER):
    print( 'Sim = ', i)
    LIST_SNR=one_realization_only_pevatrons(alpha,Kep, D, eta, KNEE)
    file_output=FOLDER_FOR_RESULTS+'/'+'results_'+str(i)+'.txt'
    print( 'LIST_SNR[0].TIME =', len(LIST_SNR[0].TIME)-1)
    save_one_LIST_to_file(LIST_SNR,file_output)


exit()
