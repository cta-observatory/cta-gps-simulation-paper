#!/usr/bin/env python3
"""
Created on Fri Apr 24 11:48:47 2020

@author: Andrea  ( andrea.giuliani@inaf.it )
"""

import sys, os
from numpy import *
from matplotlib.pyplot import *
from astropy.table import Table 
#from astropy import units as u
#from astropy.io import fits
#from astropy.coordinates import SkyCoord
import aglib as ag


def createXml_disk(srcname='source0',specfile='specfile.txt',lon=0.0,lat=0.0,radius=1.0):
  code=[
  '   <source name="source0" type="ExtendedSource""  tscalc="1">\n',
  '       <spectrum file="model0.txt" type="FileFunction">\n',
  '           <parameter free="0" max="1000.0" min="0.0" name="Normalization" scale="1.0" value="1.0" >\n',
  '       </spectrum>\n',
  '      <spatialModel type="RadialDisk">\n',
  ' 		 <parameter name="GLON"    scale="1.0" value="888" min="-360" max="360" free="0" >\n',
  ' 		 <parameter name="GLAT"   scale="1.0" value="777" min="-90"  max="90"  free="0" >\n',
  '  	 	 <parameter name="Radius" scale="1.0" value="111"    min="0.005" max="10"  free="1" >\n',
  '	</spatialModel>\n',  
   '   </source>\n' ]


  code[0] = code[0].replace('source0',srcname)  
  code[1] = code[1].replace('model0.txt',specfile)
  code[5] = code[5].replace('888',str(lon))
  code[6] = code[6].replace('777',str(lat))
  code[7] = code[7].replace('111',str(radius))

  return code

###
  
pwne = Table.read('PWNe_final_population.txt', format='ascii') 

xml=[
'<?xml version="1.0" standalone="no"?>\n',
'<source_library title="interacting supernova remnants">\n'
]

for pwn in pwne:
  
  srcname='pwn'+str(pwn['N'])
  print('  --> Adding ',srcname)   
  xml_source = createXml_disk(srcname=srcname, specfile=pwn['filename'],
                lon=pwn['GLON'], lat=pwn['GLAT'], radius=pwn['R_pwn_deg'])
  xml += xml_source 

xml += ['</source_library>\n']

outfile='pwn.xml'
mod = open(outfile, 'w')
mod.writelines(xml)
mod.close()



