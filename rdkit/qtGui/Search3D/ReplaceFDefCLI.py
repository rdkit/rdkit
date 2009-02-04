import Chem
import cPickle
import RDLogger 
import sys,getopt
logger = RDLogger.logger()

if len(sys.argv)<4:
  logger.error('Usage: ReplaceFDefCLI <old p3d filename> <fdef filename> <new p3d filename>')
  sys.exit(-1)

pklFilename=sys.argv[1]
try:
  inD = file(pklFilename,'rb').read()
except:
  logger.error('Could not open file %s for reading.'%pklFilename)
  sys.exit(-1)

fdefFilename=sys.argv[2]
try:
  fdefData = file(fdefFilename,'r').read()
except:
  logger.error('Could not read from fdef file: %s'%fdefFilename)

obj=cPickle.loads(inD)
obj['fdefData'] = fdefData

outFilename=sys.argv[3]
try:
  file(outFilename,'wb+').write(cPickle.dumps(obj,True))
except:
  logger.error('Could not open file %s for writing.'%outFilename)
  sys.exit(-1)

