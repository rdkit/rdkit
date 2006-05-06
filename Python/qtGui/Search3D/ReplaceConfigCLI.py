import cPickle
import RDLogger 
import sys,getopt,imp
from qtGui.Search3D import Persist
logger = RDLogger.logger()

if len(sys.argv)<4:
  logger.error('Usage: ReplaceConfigCLI <old p3d filename> <config filename> <new p3d filename>')
  sys.exit(-1)

pklFilename=sys.argv[1]
try:
  inD = file(pklFilename,'rb').read()
except:
  logger.error('Could not open file %s for reading.'%pklFilename)
  sys.exit(-1)

configFilename=sys.argv[2]
try:
  configMod = imp.load_source('configMod',configFilename,file(configFilename,'r'))
except:
  logger.error('Could not import parameter file: %s'%configFilename,
               exc_info=True)

obj=cPickle.loads(inD)
obj['configModule'] = Persist.ModuleVarsToDict(configMod)

outFilename=sys.argv[3]
try:
  file(outFilename,'wb+').write(cPickle.dumps(obj,True))
except:
  logger.error('Could not open file %s for writing.'%outFilename)
  sys.exit(-1)

