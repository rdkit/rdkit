#
#  Copyright (c) 2003-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" utility functionality for fingerprinting sets of molecules
 includes a command line app for working with fingerprints
 and databases


Sample Usage:

  python FingerprintMols.py  -d data.gdb \
        -t 'raw_dop_data' --smilesName="Structure" --idName="Mol_ID"  \
        --outTable="daylight_sig"


"""

import getopt
import pickle
import sys

from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys
from rdkit.ML.Cluster import Murtagh


def error(msg):
  sys.stderr.write(msg)


def message(msg):
  sys.stderr.write(msg)


def GetRDKFingerprint(mol):
  """ uses default parameters """
  details = FingerprinterDetails()
  return FingerprintMol(mol, **details.__dict__)


def FoldFingerprintToTargetDensity(fp, **fpArgs):
  nOn = fp.GetNumOnBits()
  nTot = fp.GetNumBits()
  while float(nOn) / nTot < fpArgs['tgtDensity']:
    if nTot / 2 > fpArgs['minSize']:
      fp = DataStructs.FoldFingerprint(fp, 2)
      nOn = fp.GetNumOnBits()
      nTot = fp.GetNumBits()
    else:
      break
  return fp


def FingerprintMol(mol, fingerprinter=Chem.RDKFingerprint, **fpArgs):
  if not fpArgs:
    fpArgs = FingerprinterDetails().__dict__

  if fingerprinter != Chem.RDKFingerprint:
    fp = fingerprinter(mol, **fpArgs)
    return FoldFingerprintToTargetDensity(fp, **fpArgs)

  return fingerprinter(mol, fpArgs['minPath'], fpArgs['maxPath'], fpArgs['fpSize'],
                       fpArgs['bitsPerHash'], fpArgs['useHs'], fpArgs['tgtDensity'],
                       fpArgs['minSize'])


def FingerprintsFromSmiles(dataSource, idCol, smiCol, fingerprinter=Chem.RDKFingerprint,
                           reportFreq=10, maxMols=-1, **fpArgs):
  """ fpArgs are passed as keyword arguments to the fingerprinter

  Returns a list of 2-tuples: (ID,fp)

  """
  res = []
  nDone = 0
  for entry in dataSource:
    ID, smi = str(entry[idCol]), str(entry[smiCol])
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
      fp = FingerprintMol(mol, fingerprinter, **fpArgs)
      res.append((ID, fp))
      nDone += 1
      if reportFreq > 0 and not nDone % reportFreq:
        message(f'Done {nDone} molecules\n')
      if maxMols > 0 and nDone >= maxMols:
        break
    else:
      error(f'Problems parsing SMILES: {smi}\n')
  return res


def FingerprintsFromMols(mols, fingerprinter=Chem.RDKFingerprint, reportFreq=10, maxMols=-1,
                         **fpArgs):
  """ fpArgs are passed as keyword arguments to the fingerprinter

  Returns a list of 2-tuples: (ID,fp)

  """
  res = []
  nDone = 0
  for ID, mol in mols:
    if mol:
      fp = FingerprintMol(mol, fingerprinter, **fpArgs)
      res.append((ID, fp))
      nDone += 1
      if reportFreq > 0 and not nDone % reportFreq:
        message(f'Done {nDone} molecules\n')
      if maxMols > 0 and nDone >= maxMols:
        break
    else:
      error(f'Problems parsing SMILES: {smi}\n')
  return res


def FingerprintsFromPickles(dataSource, idCol, pklCol, fingerprinter=Chem.RDKFingerprint,
                            reportFreq=10, maxMols=-1, **fpArgs):
  """ fpArgs are passed as keyword arguments to the fingerprinter

  Returns a list of 2-tuples: (ID,fp)

  """
  res = []
  nDone = 0
  for entry in dataSource:
    ID, pkl = str(entry[idCol]), str(entry[pklCol])
    mol = Chem.Mol(pkl)
    if mol is not None:
      fp = FingerprintMol(mol, fingerprinter, **fpArgs)
      res.append((ID, fp))
      nDone += 1
      if reportFreq > 0 and not nDone % reportFreq:
        message(f'Done {nDone} molecules\n')
      if maxMols > 0 and nDone >= maxMols:
        break
    else:
      error(f'Problems parsing pickle for ID: {ID}\n')
  return res


def FingerprintsFromDetails(details, reportFreq=10):
  data = None
  if details.dbName and details.tableName:
    from rdkit.Dbase import DbInfo
    from rdkit.Dbase.DbConnection import DbConnect
    from rdkit.ML.Data import DataUtils
    try:
      conn = DbConnect(details.dbName, details.tableName)
    except Exception:
      import traceback
      error(f'Problems establishing connection to database: {details.dbName}|{details.tableName}\n')
      traceback.print_exc()
    if not details.idName:
      details.idName = DbInfo.GetColumnNames(details.dbName, details.tableName)[0]
    dataSet = DataUtils.DBToData(details.dbName, details.tableName,
                                 what=f'{details.idName},{details.smilesName}')
    idCol = 0
    smiCol = 1
  elif details.inFileName and details.useSmiles:
    from rdkit.ML.Data import DataUtils
    conn = None
    if not details.idName:
      details.idName = 'ID'
    try:
      dataSet = DataUtils.TextFileToData(details.inFileName,
                                         onlyCols=[details.idName, details.smilesName])
    except IOError:
      import traceback
      error(f'Problems reading from file {details.inFileName}\n')
      traceback.print_exc()

    idCol = 0
    smiCol = 1
  elif details.inFileName and details.useSD:
    conn = None
    if not details.idName:
      details.idName = 'ID'
    dataSet = []
    try:
      s = Chem.SDMolSupplier(details.inFileName)
    except Exception:
      import traceback
      error(f'Problems reading from file {details.inFileName}\n')
      traceback.print_exc()
    else:
      while 1:
        try:
          m = s.next()
        except StopIteration:
          break
        if m:
          dataSet.append(m)
          if reportFreq > 0 and not len(dataSet) % reportFreq:
            message(f'Read {len(dataSet)} molecules\n')
            if 0 < details.maxMols <= len(dataSet):
              break

    for i, mol in enumerate(dataSet):
      if mol.HasProp(details.idName):
        nm = mol.GetProp(details.idName)
      else:
        nm = mol.GetProp('_Name')
      dataSet[i] = (nm, mol)
  else:
    dataSet = None

  fps = None
  if dataSet and not details.useSD:
    data = dataSet.GetNamedData()
    if not details.molPklName:
      fps = FingerprintsFromSmiles(data, idCol, smiCol, **details.__dict__)
    else:
      fps = FingerprintsFromPickles(data, idCol, smiCol, **details.__dict__)
  elif dataSet and details.useSD:
    fps = FingerprintsFromMols(dataSet, **details.__dict__)

  if fps:
    if details.outFileName:
      outF = open(details.outFileName, 'wb+')
      for i in range(len(fps)):
        pickle.dump(fps[i], outF)
      outF.close()
    dbName = details.outDbName or details.dbName
    if details.outTableName and dbName:
      from rdkit.Dbase import DbModule, DbUtils
      from rdkit.Dbase.DbConnection import DbConnect
      conn = DbConnect(dbName)
      #
      #  We don't have a db open already, so we'll need to figure out
      #    the types of our columns...
      #
      colTypes = DbUtils.TypeFinder(data, len(data), len(data[0]))
      typeStrs = DbUtils.GetTypeStrings([details.idName, details.smilesName], colTypes,
                                        keyCol=details.idName)
      cols = f'{typeStrs[0]}, {details.fpColName} {DbModule.binaryTypeName}'

      # FIX: we should really check to see if the table
      #  is already there and, if so, add the appropriate
      #  column.

      #
      # create the new table
      #
      if details.replaceTable or \
         details.outTableName.upper() not in [x.upper() for x in conn.GetTableNames()]:
        conn.AddTable(details.outTableName, cols)

      #
      # And add the data
      #
      for ID, fp in fps:
        tpl = ID, DbModule.binaryHolder(fp.ToBinary())
        conn.InsertData(details.outTableName, tpl)
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
    self.fingerprinter = Chem.RDKFingerprint
    self.fpColName = "AutoFragmentFP"
    self.idName = ''
    self.dbName = ''
    self.outDbName = ''
    self.tableName = ''
    self.minSize = 64
    self.fpSize = 2048
    self.tgtDensity = 0.3
    self.minPath = 1
    self.maxPath = 7
    self.discrimHash = 0
    self.useHs = 0
    self.useValence = 0
    self.bitsPerHash = 2
    self.smilesName = 'SMILES'
    self.maxMols = -1
    self.outFileName = ''
    self.outTableName = ''
    self.inFileName = ''
    self.replaceTable = True
    self.molPklName = ''
    self.useSmiles = True
    self.useSD = False

  def _screenerInit(self):
    self.metric = DataStructs.TanimotoSimilarity
    self.doScreen = ''
    self.topN = 10
    self.screenThresh = 0.75
    self.doThreshold = 0
    self.smilesTableName = ''
    self.probeSmiles = ''
    self.probeMol = None
    self.noPickle = 0

  def _clusterInit(self):
    self.clusterAlgo = Murtagh.WARDS
    self.actTableName = ''
    self.actName = ''

  def GetMetricName(self):
    # DataStructs.TverskySimilarity: 'Tversky'
    metricDict = {
      DataStructs.DiceSimilarity: 'Dice',
      DataStructs.TanimotoSimilarity: 'Tanimoto',
      DataStructs.CosineSimilarity: 'Cosine',
    }
    metric = metricDict.get(self.metric, self.metric)
    if metric:
      return metric
    return 'Unknown'

  def SetMetricFromName(self, name):
    # 'TVERSKY': DataStructs.TverskySimilarity,
    metricDict = {
      'DICE': DataStructs.DiceSimilarity,
      'TANIMOTO': DataStructs.TanimotoSimilarity,
      'COSINE': DataStructs.CosineSimilarity,
    }
    self.metric = metricDict.get(name.upper(), self.metric)


def Usage():
  """  prints a usage string and exits

  """
  print(_usageDoc)
  sys.exit(-1)


_usageDoc = """
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

    - --useSD:  Assume that the input file is an SD file, not a SMILES
       table.

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
      fingerprint for each fragment. Default is *2*.

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
  args = sys.argv[1:]
  try:
    args, extras = getopt.getopt(
      args,
      'HVs:d:t:o:h',
      [
        'minSize=',
        'maxSize=',
        'density=',
        'minPath=',
        'maxPath=',
        'bitsPerHash=',
        'smilesName=',
        'molPkl=',
        'useSD',
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
  except Exception:
    import traceback
    traceback.print_exc()
    Usage()

  if details is None:
    details = FingerprinterDetails()
  if len(extras):
    details.inFileName = extras[0]

  for arg, val in args:
    if arg == '-H':
      details.useHs = 1
    elif arg == '-V':
      details.useValence = 1
    elif arg == '-d':
      details.dbName = val
    elif arg == '-t':
      details.tableName = val
    elif arg == '-o':
      details.outFileName = val
    elif arg == '--minSize':
      details.minSize = int(val)
    elif arg == '--maxSize':
      details.fpSize = int(val)
    elif arg == '--density':
      details.tgtDensity = float(val)
    elif arg == '--outTable':
      details.outTableName = val
    elif arg == '--outDbName':
      details.outDbName = val
    elif arg == '--fpColName':
      details.fpColName = val
    elif arg == '--minPath':
      details.minPath = int(val)
    elif arg == '--maxPath':
      details.maxPath = int(val)
    elif arg == '--nBitsPerHash':
      details.bitsPerHash = int(val)
    elif arg == '--discrim':
      details.discrimHash = 1
    elif arg == '--smilesName':
      details.smilesName = val
    elif arg == '--molPkl':
      details.molPklName = val
    elif arg == '--useSD':
      details.useSmiles = False
      details.useSD = True
    elif arg == '--idName':
      details.idName = val
    elif arg == '--maxMols':
      details.maxMols = int(val)
    elif arg == '--useMACCS':
      details.fingerprinter = MACCSkeys.GenMACCSKeys
    elif arg == '--keepTable':
      details.replaceTable = False

    # SCREENER:
    elif arg == '--smilesTable':
      details.smilesTableName = val
    elif arg == '--topN':
      details.doThreshold = 0
      details.topN = int(val)
    elif arg == '--thresh':
      details.doThreshold = 1
      details.screenThresh = float(val)
    elif arg == '--smiles':
      details.probeSmiles = val
    elif arg == '--dice':
      details.metric = DataStructs.DiceSimilarity
    elif arg == '--cosine':
      details.metric = DataStructs.CosineSimilarity

    # CLUSTERS:
    elif arg == '--SLINK':
      details.clusterAlgo = Murtagh.SLINK
    elif arg == '--CLINK':
      details.clusterAlgo = Murtagh.CLINK
    elif arg == '--UPGMA':
      details.clusterAlgo = Murtagh.UPGMA
    elif arg == '--actTable':
      details.actTableName = val
    elif arg == '--actName':
      details.actName = val
    elif arg == '-h':
      Usage()
  return details


if __name__ == '__main__':
  message("This is FingerprintMols\n\n")
  details = ParseArgs()
  FingerprintsFromDetails(details)
