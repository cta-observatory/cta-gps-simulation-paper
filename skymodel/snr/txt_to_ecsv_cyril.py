"""
Reformat data for simulated SNRs by Pierre Cristofari from TXT to ECSV
"""
import numpy as np
import astropy.units as u
from astropy.table import Table, Column
from gammapy.utils.coordinates import galactic as compute_galactic_coordinates
import numpy as np
import os



def test_not_empty(version):
    filename = 'ALL_FILES_0/results_{}.txt'.format(version)
    data2 = np.genfromtxt(filename, defaultfmt='%0.6f',skip_header=1)

    if(len(data2)==0):
        return 0
    else:
        return 1

def read_txt_files(version):
    filename = 'ALL_FILES_0/results_{}.txt'.format(version)
    print('Reading {}'.format(filename))
    data = np.genfromtxt(filename, defaultfmt='%0.6f', names=True)

    # Trick to read the Energy grid
    data2 = np.genfromtxt(filename, defaultfmt='%0.6f',skip_header=1)
    

#    if (len(data2)==0):
#        print( 'HOIHIFH')
#        t=[]
#        return t
    #    if (len(data2[:,0]) = 0 ):
#        print( ' EMPTTY ')
#        exit()

    EGRID = data2[:,11]     #Ag
    p=0
#    while(data2[p][0]==0 ):
#        p=p+1
    p=40
    EGRID_SMALL=[]
    for i in range (0,p):
        EGRID_SMALL.append(EGRID[i])


    number_of_SNRs = int(data2[len(data2[:,0])-1][0]+1) # number of SNRs

    print( 'number of SNRs =', number_of_SNRs)
    size = len(data.dtype.names)

    # Reformat spectral data into arrays
    test = [None] * number_of_SNRs
    for j in range(0, number_of_SNRs):
        test[j] = [data[j*len(EGRID_SMALL)][0]]
        for i in range(1, size-4):   # -2 to not take the E and diff spectrum column    #Ag
            test[j] = np.append(test[j], data[j*len(EGRID_SMALL)][i])
        
        
        for e in range (j*len(EGRID_SMALL),j*len(EGRID_SMALL)+len(EGRID_SMALL)):
            test[j]=np.append(test[j], data[e][size-1])
    #   print( 'test[j]=',  test[j])
        #  print ('len = ,', len(test[j]))
#   print ('-------')

    NAMES = data.dtype.names[1]
    temp = test[0]
    for j in range(0, number_of_SNRs):
        temp = np.vstack([temp, test[j]])
    
        NAMES = data.dtype.names[0]
        for x in range(1, size-4):       #Ag
            NAMES = np.append(NAMES, data.dtype.names[x])
        for x in range(0, len(EGRID_SMALL)):
            NAMES = np.append(NAMES, EGRID_SMALL[x])

    print( ' NAMES ==', NAMES)
#    exit()

    EGRID_SMALL=np.array(EGRID_SMALL)
    t = Table(temp, names=NAMES)

    # Parameters used for this simulation
    t.meta['simulation_parameters'] = [
        dict(name='alpha', value=4.1, description='Spectral index'),
        dict(name='Kep', value=1e-2, description='Electron to proton ratio'),
        dict(name='xi', value=0.1, description='CR efficiency'),
    ]

    t.meta['energy_array'] = EGRID_SMALL.tolist()

    t[NAMES[0]].unit = ''
    t[NAMES[0]].format = '%i'
    t[NAMES[0]].description = 'Simulation number'

    t[NAMES[1]].unit = 'kpc'
    t[NAMES[1]].format = '%0.4f'
    t[NAMES[1]].description = 'X Position'

    t[NAMES[2]].unit = 'kpc'
    t[NAMES[2]].format = '%0.4f'
    t[NAMES[2]].description = 'Y Position'

    t[NAMES[3]].unit = 'kpc'
    t[NAMES[3]].format = '%0.4f'
    t[NAMES[3]].description = 'Z Position'

    t[NAMES[4]].unit = u.Unit('cm-3')
    t[NAMES[4]].format = '%0.4f'
    t[NAMES[4]].description = 'ISM density'

    t[NAMES[5]].unit = u.Unit('cm-3')                  #Ag adding this field
    t[NAMES[5]].format = '%0.4f'
    t[NAMES[5]].description = 'Shock density'

    t[NAMES[6]].unit = ''
    t[NAMES[6]].format = '%i'
    t[NAMES[6]].description = 'type of progenitor'

    t[NAMES[7]].unit = u.Unit('solMass')                  #Ag adding this field
    t[NAMES[7]].format = '%0.4f'
    t[NAMES[7]].description = 'Mass of Ejecta'

    t[NAMES[8]].unit = 'kyear'
    t[NAMES[8]].format = '%0.4f'
    t[NAMES[8]].description = 'age of the SNR'

    t[NAMES[9]].unit = 'arcmin'
    t[NAMES[9]].format = '%0.4f'
    t[NAMES[9]].description = 'apparent size of the SNR'

    t[NAMES[10]].unit = 'pc'
    t[NAMES[10]].format = '%0.4f'
    t[NAMES[10]].description = 'SNR Radius'


    for x in range(11, 11+p):       #Ag
        t[NAMES[x]].unit = 'cm-2 s-1 TeV'
        t[NAMES[x]].description = 'Differential spectrum'

    return t


def add_extra_info(table):
    # Change Pierre's XYZ to the one used in Gammapy at the moment
    # This was checked to be correct in https://github.com/gammasky/cta-dc/issues/17
    table['galactocentric_x'] = Column(table['POS_Y'].data, unit='kpc', description='Galactocentric X', format='%0.5f')
    table['galactocentric_y'] = Column(-table['POS_X'].data, unit='kpc', description='Galactocentric Y', format='%0.5f')
    table['galactocentric_z'] = Column(table['POS_Z'].data, unit='kpc', description='Galactocentric Y', format='%0.5f')
    table.remove_columns(['POS_X', 'POS_Y', 'POS_Z'])
#    table['flux_1_10'] = Column(table['flux_1_10'].data, unit='cm-2 s-1', description='Integral flux from 1 to 10 TeV',
#                                format='%0.4e')

#  table.rename_column('Radius', 'size_physical')
    table.rename_column('size', 'sigma')

    r = np.sqrt(table['galactocentric_x'] ** 2 + table['galactocentric_y'] ** 2)
    table['galactocentric_r'] = Column(r, unit='kpc', description='Galactocentric radius in the xy plan')

    distance, glon, glat = compute_galactic_coordinates(
        x=table['galactocentric_x'].quantity,
        y=table['galactocentric_y'].quantity,
        z=table['galactocentric_z'].quantity,
    )
    table['distance'] = Column(distance, unit='kpc', description='Distance from Earth')
    table['distance'].format = '%.5f'

    table['glon'] = Column(glon, unit='deg', description='Galactic longitude')
    table['glon'].format = '%.5f'
    table['glat'] = Column(glat, unit='deg', description='Galactic latitude')
    table['glat'].format = '%.5f'
    zarray = np.zeros(len(table['distance']))
    table['skip'] = Column(zarray, description='Skip boolean, 1 skip 0 keep')
    table['skip'].format = '%.d'

    #print(table['flux_1_10'])
    #print(table.info())
    return table


#
#filename = 'ALL_FILES_0/results_0.txt'
##print('Reading {}'.format(filename))
#data = np.genfromtxt(filename, defaultfmt='%0.6f', names=True)
#
    # Trick to read the Energy grid
#data2 = np.genfromtxt(filename, defaultfmt='%0.6f',skip_header=1)

#EGRID = data2[:,9]
#print ('EGRID', EGRID)
#p=0
#while(data2[p][0]==0):
#    p=p+1
#
#EGRID_SMALL=[]
#print('last p =' , p )
#for i in range (0,p):
#    print ('EGRID ', EGRID[i])
#    EGRID_SMALL.append(EGRID[i])
#
#print (' EGRID_SMALL' , EGRID_SMALL)
#
#
#print ( ' len(data2[:,0]) = ' , len(data2[:,0]))
#print( 'data2[29,0]', int(data2[29,0]))
#number_of_SNRs = int(data2[len(data2[:,0])-1][0]+1) # number of SNRs
#print( 'number_of_SNRs  = ', number_of_SNRs)
#size = len(data.dtype.names)
#print ('size =', size)
#size_and_EGRID

    # Reformat spectral data into arrays
#test = [None] * number_of_SNRs
#for j in range(0, number_of_SNRs):
#    test[j] = [data[j*len(EGRID_SMALL)][0]]
#    for i in range(1, size-2):   # -2 to not take the E and diff spectrum column
#        test[j] = np.append(test[j], data[j+len(EGRID_SMALL)][i])
#    for e in range (j*len(EGRID_SMALL),j*len(EGRID_SMALL)+len(EGRID_SMALL)):
#        test[j]=np.append(test[j], data[e][size-1])
#    print( 'test[j]=',  test[j])
#    print ('len = ,', len(test[j]))
#    print ('-------')
#
#
#NAMES = data.dtype.names[1]
#print ('NAMES = ',NAMES)
#temp = test[0]
#for j in range(1, number_of_SNRs):
#    temp = np.vstack([temp, test[j]])
#
#    NAMES = data.dtype.names[0]
#    for x in range(1, size-2):
#        NAMES = np.append(NAMES, data.dtype.names[x])
#    for x in range(0, len(EGRID_SMALL)):
#        NAMES = np.append(NAMES, EGRID_SMALL[x])
#





FOLDER_FOR_RESULTS='OUTPUT_FILES_1'
if not os.path.exists(FOLDER_FOR_RESULTS):
    os.makedirs(FOLDER_FOR_RESULTS)

if __name__ == '__main__':
    for version in range (0,1): # to 100 for cyril

        not_empty=test_not_empty(version)
        if ( not_empty==1):
            table = read_txt_files(version=version)
            table = add_extra_info(table)
            
            filename = FOLDER_FOR_RESULTS+'/ctadc_skymodel_gps_sources_pevatron_{}.ecsv'.format(version)
            print('Writing {}'.format(filename))
            table.write(filename, format='ascii.ecsv', overwrite=True)
