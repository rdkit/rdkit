#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#

def parseXPM(fileName):
  inD = open(fileName,'r').read()
  startIdx = inD.find('{')+1
  endIdx = inD.rfind('}')
  return inD[startIdx:endIdx]

header="""#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#

# icons for use in the visual programming GUI
"""  

def go(fileNames,outFName='VPIcons.py'):
  outF = open(outFName,'w+')
  outF.write(header)
  for fileName in fileNames:
    base = fileName.split('.xpm')[0]
    xpmData = parseXPM(fileName)
    resName = '%sIcon_xpmData'%base
    outF.write('%s=[%s]\n'%(resName,xpmData))
  outF.close()

if __name__=='__main__':
  import sys
  go(sys.argv[1:])
    
    
  
  
  
  
