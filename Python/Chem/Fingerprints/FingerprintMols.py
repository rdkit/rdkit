# $Id$
#
#  Copyright (c) 2003-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" utility functionality for fingerprinting sets of molecules
 includes a command line app for working with fingerprints
 and databases


Sample Usage:

  python FingerprintMols.py  -d data.gdb \
        -t 'raw_dop_data' --smilesName="Structure" --idName="Mol_ID"  \
        --outTable="daylight_sig"


"""
import Chem
from Chem import MACCSkeys
from Dbase.DbConnection import DbConnect
from Dbase import DbInfo,DbUtils,DbModule
from ML.Data import DataUtils
from ML.Cluster import Murtagh
import DataStructs
import sys
import cPickle

_cvsVersion="$Id$"
idx1 = _cvsVersion.find(':')+1
idx2 = _cvsVersion.rfind('$')
__VERSION_STRING="%s"%(_cvsVersion[idx1:idx2])


def error(msg):
  sys.stderr.write(msg)
def message(msg):
  sys.stderr.write(msg)


def FingerprintMol(mol,
                   fingerprinter=Chem.DaylightFingerprint,
                   **fpArgs):
  if fingerprinter != Chem.DaylightFingerprint:
    fp = fingerprinter(mol,**fpArgs)
  else:
    fp = fingerprinter(mol,fpArgs['minPath'],fpArgs['maxPath'],
                       fpArgs['fpSize'],fpArgs['bitsPerHash'],
                       fpArgs['useHs'])
  nOn = fp.GetNumOnBits()
  nTot = fp.GetNumBits()
  while( float(nOn)/nTot < fpArgs['tgtDensity'] ):
    if nTot / 2 > fpArgs['minSize']:
      fp = DataStructs.FoldFingerprint(fp,2)
      nOn = fp.GetNumOnBits()
      nTot = fp.GetNumBits()
    else:
      break
  return fp

def FingerprintsFromSmiles(dataSource,idCol,smiCol,
                           fingerprinter=Chem.DaylightFingerprint,
                           reportFreq=10,maxMols=-1,
                           **fpArgs):
  """ fpArgs are passed as keyword arguments to the fingerprinter

  Returns a list of 2-tuples: (id,fp)
  
  """
  res = []
  nDone = 0
  for entry in dataSource:
    id,smi = str(entry[idCol]),str(entry[smiCol])
    try:
      mol = Chem.MolFromSmiles(smi)
    except:
      mol = None
    if mol:
      fp = FingerprintMol(mol,fingerprinter,**fpArgs)
      res.append((id,fp))
      nDone += 1
      if reportFreq>0 and not nDone % reportFreq:
        message('Done %d molecules\n'%(nDone))
      if maxMols > 0 and nDone >= maxMols:
        break
    else:
      error('Problems parsing SMILES: %s\n'%smi)
  return res

def FingerprintsFromDetails(details):
  data = None
  if details.dbName and details.tableName:
    try:
      conn = DbConnect(details.dbName,details.tableName)
    except:
      import traceback
      error('Problems establishing connection to database: %s|%s\n'%(details.dbName,
                                                                     details.tableName))
      traceback.print_exc()
    if not details.idName:
      details.idName=DbInfo.GetColumnNames(details.dbName,details.tableName)[0]
    dataSet = DataUtils.DBToData(details.dbName,details.tableName,
                                 what='%s,%s'%(details.idName,details.smilesName))
    idCol = 0
    smiCol = 1
  elif details.inFileName:
    conn = None
    if not details.idName:
      details.idName='ID'
    try:
      dataSet = DataUtils.TextFileToData(details.inFileName,
                                         onlyCols=[details.idName,details.smilesName])
    except IOError:
      import traceback
      error('Problems reading from file %s\n'%(details.inFileName))
      traceback.print_exc()
      
    idCol = 0
    smiCol = 1
  else:
     dataSet = None
      
  fps = None
  if dataSet:
    data = dataSet.GetNamedData()
    fps = apply(FingerprintsFromSmiles,(data,idCol,smiCol),
                details.__dict__)
  if fps:
    if details.outFileName:
      outF = open(details.outFileName,'wb+')
      for i in range(len(fps)):
        cPickle.dump(fps[i],outF)
      outF.close()
    dbName = details.outDbName or details.dbName
    if details.outTableName and dbName:
      conn = DbConnect(dbName)
      #
      #  We don't have a db open already, so we'll need to figure out
      #    the types of our columns...
      #
      colTypes = DbUtils.TypeFinder(data,len(data),len(data[0]))
      typeStrs = DbUtils.GetTypeStrings([details.idName,details.smilesName],colTypes,
                                        keyCol=details.idName)
      cols = '%s, %s %s'%(typeStrs[0],details.fpColName,DbModule.binaryTypeName)

      # FIX: we should really check to see if the table
      #  is already there and, if so, add the appropriate
      #  column.

      #
      # create the new table
      #
      if details.replaceTable or \
         details.outTableName.upper() not in [x.upper() for x in conn.GetTableNames()]:
        conn.AddTable(details.outTableName,cols)
      
      #
      # And add the data
      #
      for id,fp in fps:
        tpl = id,DbModule.binaryHolder(fp.ToBinary())
        conn.InsertData(details.outTableName,tpl)
      conn.Commit()  
  return fps
# ------------------------------------------------
#
#  Command line parsing stuff
#
# ------------------------------------------------

class FingerprinterDetails(object):
  """ class for storing the details of a fingerprinting run,
     generates sensible defaults on construction

  """
  def __init__(self):
    self._fingerprinterInit()
    self._screenerInit()
    self._clusterInit()

  def _fingerprinterInit(self):
    self.fingerprinter = Chem.DaylightFingerprint
    self.fpColName="AutoFragmentFP"
    self.idName=''
    self.dbName=''
    self.outDbName=''
    self.tableName=''
    self.minSize=64
    self.fpSize=2048
    self.tgtDensity=0.3
    self.minPath=1
    self.maxPath=7
    self.discrimHash=0
    self.useHs=0
    self.useValence=0
    self.bitsPerHash=4
    self.smilesName='SMILES'
    self.maxMols=-1
    self.outFileName=''
    self.outTableName=''
    self.inFileName=''
    self.replaceTable=True

  def _screenerInit(self):
    self.metric = DataStructs.TanimotoSimilarity
    self.doScreen=''
    self.topN=10
    self.screenThresh=0.75
    self.doThreshold=0
    self.smilesTableName=''
    self.probeSmiles=''
    self.probeMol=None
    self.noPickle=0
    
  def _clusterInit(self):
    self.clusterAlgo = Murtagh.WARDS
    self.actTableName = ''
    self.actName = ''
    
  def GetMetricName(self):
    if self.metric == DataStructs.TanimotoSimilarity:
      return 'Tanimoto'
    elif self.metric == DataStructs.DiceSimilarity:
      return 'Dice'
    elif self.metric == DataStructs.CosineSimilarity:
      return 'Cosine'
    elif self.metric:
      return self.metric
    else:
      return 'Unknown'
  def SetMetricFromName(self,name):
    name = name.upper()
    if name=="TANIMOTO":
      self.metric = DataStructs.TanimotoSimilarity
    elif name=="DICE":
      self.metric = DataStructs.DiceSimilarity
    elif name=="COSINE":
      self.metric = DataStructs.CosineSimilarity
    
def Usage():
  """  prints a usage string and exits

  """
  print _usageDoc
  sys.exit(-1)

_usageDoc="""
Usage: FingerprintMols.py [args] <fName>

  If <fName> is provided and no tableName is specified (see below),
  data will be read from the text file <fName>.  Text files delimited
  with either commas (extension .csv) or tabs (extension .txt) are
  supported. 

  Command line arguments are:
    - -d _dbName_: set the name of the database from which
      to pull input molecule information.  If output is 
      going to a database, this will also be used for that
      unless the --outDbName option is used.

    - -t _tableName_: set the name of the database table 
      from which to pull input molecule information

    - --smilesName=val: sets the name of the SMILES column
      in the input database.  Default is *SMILES*.

    - --idName=val: sets the name of the id column in the input
      database.  Defaults to be the name of the first db column
      (or *ID* for text files).
      
    - -o _outFileName_:  name of the output file (output will
      be a pickle file with one label,fingerprint entry for each
      molecule).  
      
    - --outTable=val: name of the output db table used to store
      fingerprints.  If this table already exists, it will be
      replaced.

    - --outDbName: name of output database, if it's being used.
      Defaults to be the same as the input db.

    - --fpColName=val: name to use for the column which stores
      fingerprints (in pickled format) in the output db table.
      Default is *AutoFragmentFP*
      
    - --maxSize=val:  base size of the fingerprints to be generated
      Default is *2048*

    - --minSize=val: minimum size of the fingerprints to be generated
      (limits the amount of folding that happens).  Default is *64*

    - --density=val: target bit density in the fingerprint.  The
      fingerprint will be folded until this density is
      reached. Default is *0.3*

    - --minPath=val:  minimum path length to be included in
      fragment-based fingerprints. Default is *1*.

    - --maxPath=val:  maximum path length to be included in
      fragment-based fingerprints. Default is *7*.
      
    - --nBitsPerHash: number of bits to be set in the output
      fingerprint for each fragment. Default is *4*.

    - --discrim: use of path-based discriminators to hash bits.
      Default is *false*.

    - -V: include valence information in the fingerprints
      Default is *false*.
      
    - -H: include Hs in the fingerprint 
      Default is *false*.

    - --maxMols=val: sets the maximum number of molecules to be
      fingerprinted.

    - --useMACCS: use the public MACCS keys to do the fingerprinting
      (instead of a daylight-type fingerprint)

"""

def ParseArgs(details=None):
  """ parses the command line arguments and returns a
   _FingerprinterDetails_ instance with the results.

   **Note**:

     - If you make modifications here, please update the global
       _usageDoc string so the Usage message is up to date.

     - This routine is used by both the fingerprinter, the clusterer and the
       screener; not all arguments make sense for all applications.

  """
  import sys,getopt
  try:
    args = sys.argv[1:]
  except:
    Usage()
  try:
    args,extras = getopt.getopt(args,'HVs:d:t:o:h',
                                [
                                 'minSize=','maxSize=',
                                 'density=',
                                 'minPath=','maxPath=',
                                 'bitsPerHash=',
                                 'smilesName=',
                                 'idName=',
                                 'discrim',
                                 'outTable=',
                                 'outDbName=',
                                 'fpColName=',
                                 'maxMols=',
                                 'useMACCS',
                                 'keepTable',
                                 # SCREENING:
                                 'smilesTable=',
                                 'doScreen=',
                                 'topN=',
                                 'thresh=',
                                 'smiles=',
                                 'dice',
                                 'cosine',
                                 # CLUSTERING:
                                 'actTable=',
                                 'actName=',
                                 'SLINK',
                                 'CLINK',
                                 'UPGMA',

                                 ])
  except:
    import traceback
    traceback.print_exc()
    Usage()
                             
  if details is None:
    details = FingerprinterDetails()
  if len(extras):
    details.inFileName=extras[0]
    
  for arg,val in args:
    if arg=='-H':
      details.useHs=1
    elif arg=='-V':
      details.useValence=1
    elif arg=='-d':
      details.dbName = val
    elif arg=='-t':
      details.tableName = val
    elif arg=='-o':
      details.outFileName = val
    elif arg=='--minSize':
      details.minSize= int(val)
    elif arg=='--maxSize':
      details.fpSize= int(val)
    elif arg=='--density':
      details.tgtDensity = float(val)
    elif arg=='--outTable':
      details.outTableName = val
    elif arg=='--outDbName':
      details.outDbName = val
    elif arg=='--fpColName':
      details.fpColName = val
    elif arg=='--minPath':
      details.minPath= int(val)
    elif arg=='--maxPath':
      details.maxPath= int(val)
    elif arg=='--nBitsPerHash':
      details.bitsPerHash= int(val)
    elif arg=='--discrim':
      details.discrimHash=1
    elif arg=='--smilesName':
      details.smilesName = val
    elif arg=='--idName':
      details.idName = val
    elif arg=='--maxMols':
      details.maxMols = int(val)
    elif arg=='--useMACCS':
      details.fingerprinter = MACCSkeys.GenMACCSKeys
    elif arg=='--keepTable':
      details.replaceTable=False

    # SCREENER:
    elif arg=='--smilesTable':
      details.smilesTableName=val;
    elif arg=='--topN':
      details.doThreshold=0
      details.topN=int(val)
    elif arg=='--thresh':
      details.doThreshold=1
      details.screenThresh=float(val)
    elif arg=='--smiles':
      details.probeSmiles=val;
    elif arg=='--dice':
      details.metric = DataStructs.DiceSimilarity
    elif arg=='--cosine':
      details.metric = DataStructs.CosineSimilarity

    # CLUSTERS:
    elif arg=='--SLINK':
      details.clusterAlgo = Murtagh.SLINK
    elif arg=='--CLINK':
      details.clusterAlgo = Murtagh.CLINK
    elif arg=='--UPGMA':
      details.clusterAlgo = Murtagh.UPGMA
    elif arg=='--actTable':
      details.actTableName = val
    elif arg=='--actName':
      details.actName = val
    elif arg=='-h':
      Usage()
  return details

if __name__ == '__main__':
  message("This is FingerprintMols version %s\n\n"%(__VERSION_STRING))
  details = ParseArgs()
  FingerprintsFromDetails(details)

    



