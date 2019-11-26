import sys
import gammalib 

#xmlfile='out/isnr.xml'
xmlfile=sys.argv[1]

try:
  q=gammalib.GModels(xmlfile)
  mss='  '+xmlfile+' loaded!'
  print(mss)
  
except:
  mss='  '+'problem in loading '+ xmlfile
  print(mss)
  
