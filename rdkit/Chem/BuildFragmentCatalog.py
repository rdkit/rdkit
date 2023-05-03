#
#  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""  command line utility for working with FragmentCatalogs (CASE-type analysis)

**Usage**

  BuildFragmentCatalog [optional args] <filename>

 filename, the name of a delimited text file containing InData, is required
 for some modes of operation (see below)

**Command Line Arguments**

 - -n *maxNumMols*:  specify the maximum number of molecules to be processed

 - -b: build the catalog and OnBitLists
    *requires InData*

 - -s: score compounds
    *requires InData and a Catalog, can use OnBitLists*

 - -g: calculate info gains
    *requires Scores*

 - -d: show details about high-ranking fragments
    *requires a Catalog and Gains*

 - --catalog=*filename*: filename with the pickled catalog.
    If -b is provided, this file will be overwritten.

 - --onbits=*filename*: filename to hold the pickled OnBitLists.
   If -b is provided, this file will be overwritten

 - --scores=*filename*: filename to hold the text score data.
   If -s is provided, this file will be overwritten

 - --gains=*filename*: filename to hold the text gains data.
   If -g is provided, this file will be overwritten

 - --details=*filename*: filename to hold the text details data.
   If -d is provided, this file will be overwritten.

 - --minPath=2: specify the minimum length for a path

 - --maxPath=6: specify the maximum length for a path

 - --smiCol=1: specify which column in the input data file contains
     SMILES

 - --actCol=-1: specify which column in the input data file contains
     activities

 - --nActs=2: specify the number of possible activity values

 - --nBits=-1: specify the maximum number of bits to show details for

"""

import os
import pickle
import sys

import numpy

from rdkit import RDConfig
from rdkit.Chem import FragmentCatalog
from rdkit.Dbase.DbConnection import DbConnect
from rdkit.ML import InfoTheory


def message(msg, dest=sys.stdout):
  dest.write(msg)


def BuildCatalog(suppl, maxPts=-1, groupFileName=None, minPath=2, maxPath=6, reportFreq=10):
  """ builds a fragment catalog from a set of molecules in a delimited text block

      **Arguments**

        - suppl: a mol supplier

        - maxPts: (optional) if provided, this will set an upper bound on the
          number of points to be considered

        - groupFileName: (optional) name of the file containing functional group
          information

        - minPath, maxPath: (optional) names of the minimum and maximum path lengths
          to be considered

        - reportFreq: (optional) how often to display status information

      **Returns**

        a FragmentCatalog

    """
  if groupFileName is None:
    groupFileName = os.path.join(RDConfig.RDDataDir, "FunctionalGroups.txt")

  fpParams = FragmentCatalog.FragCatParams(minPath, maxPath, groupFileName)
  catalog = FragmentCatalog.FragCatalog(fpParams)
  fgen = FragmentCatalog.FragCatGenerator()
  if maxPts > 0:
    nPts = maxPts
  else:
    if hasattr(suppl, '__len__'):
      nPts = len(suppl)
    else:
      nPts = -1
  for i, mol in enumerate(suppl):
    if i == nPts:
      break
    if i and not i % reportFreq:
      if nPts > -1:
        message('Done %d of %d, %d paths\n' % (i, nPts, catalog.GetFPLength()))
      else:
        message('Done %d, %d paths\n' % (i, catalog.GetFPLength()))
    fgen.AddFragsFromMol(mol, catalog)
  return catalog


def ScoreMolecules(suppl, catalog, maxPts=-1, actName='', acts=None, nActs=2, reportFreq=10):
  """ scores the compounds in a supplier using a catalog

      **Arguments**

        - suppl: a mol supplier

        - catalog: the FragmentCatalog

        - maxPts: (optional) the maximum number of molecules to be
          considered

        - actName: (optional) the name of the molecule's activity property.
          If this is not provided, the molecule's last property will be used.

        - acts: (optional) a sequence of activity values (integers).
          If not provided, the activities will be read from the molecules.

        - nActs: (optional) number of possible activity values

        - reportFreq: (optional) how often to display status information

      **Returns**

        a 2-tuple:

          1) the results table (a 3D array of ints nBits x 2 x nActs)

          2) a list containing the on bit lists for each molecule

    """
  nBits = catalog.GetFPLength()
  resTbl = numpy.zeros((nBits, 2, nActs), numpy.int32)
  obls = []

  if not actName and not acts:
    actName = suppl[0].GetPropNames()[-1]

  fpgen = FragmentCatalog.FragFPGenerator()
  suppl.reset()
  i = 1
  for mol in suppl:
    if i and not i % reportFreq:
      message('Done %d.\n' % (i))
    if mol:
      if not acts:
        act = int(mol.GetProp(actName))
      else:
        act = acts[i - 1]
      fp = fpgen.GetFPForMol(mol, catalog)
      obls.append([x for x in fp.GetOnBits()])
      for j in range(nBits):
        resTbl[j, 0, act] += 1
      for id_ in obls[i - 1]:
        resTbl[id_ - 1, 0, act] -= 1
        resTbl[id_ - 1, 1, act] += 1
    else:
      obls.append([])
    i += 1
  return resTbl, obls


def ScoreFromLists(bitLists, suppl, catalog, maxPts=-1, actName='', acts=None, nActs=2,
                   reportFreq=10):
  """  similar to _ScoreMolecules()_, but uses pre-calculated bit lists
      for the molecules (this speeds things up a lot)


      **Arguments**

        - bitLists: sequence of on bit sequences for the input molecules

        - suppl: the input supplier (we read activities from here)

        - catalog: the FragmentCatalog

        - maxPts: (optional) the maximum number of molecules to be
          considered

        - actName: (optional) the name of the molecule's activity property.
          If this is not provided, the molecule's last property will be used.

        - nActs: (optional) number of possible activity values

        - reportFreq: (optional) how often to display status information

      **Returns**

         the results table (a 3D array of ints nBits x 2 x nActs)

    """
  nBits = catalog.GetFPLength()
  if maxPts > 0:
    nPts = maxPts
  else:
    nPts = len(bitLists)
  resTbl = numpy.zeros((nBits, 2, nActs), numpy.int32)
  if not actName and not acts:
    actName = suppl[0].GetPropNames()[-1]
  suppl.reset()
  for i in range(1, nPts + 1):
    mol = next(suppl)
    if not acts:
      act = int(mol.GetProp(actName))
    else:
      act = acts[i - 1]
    if i and not i % reportFreq:
      message('Done %d of %d\n' % (i, nPts))
    ids = set()
    for id_ in bitLists[i - 1]:
      ids.add(id_ - 1)
    for j in range(nBits):
      resTbl[j, 0, act] += 1
    for id_ in ids:
      resTbl[id_, 0, act] -= 1
      resTbl[id_, 1, act] += 1
  return resTbl


def CalcGains(suppl, catalog, topN=-1, actName='', acts=None, nActs=2, reportFreq=10, biasList=None,
              collectFps=0):
  """ calculates info gains by constructing fingerprints
      *DOC*

      Returns a 2-tuple:
         1) gains matrix
         2) list of fingerprints

    """
  nBits = catalog.GetFPLength()
  if topN < 0:
    topN = nBits
  if not actName and not acts:
    actName = suppl[0].GetPropNames()[-1]

  if hasattr(suppl, '__len__'):
    nMols = len(suppl)
  else:
    nMols = -1
  fpgen = FragmentCatalog.FragFPGenerator()
  # ranker = InfoTheory.InfoBitRanker(nBits,nActs,InfoTheory.InfoType.ENTROPY)
  if biasList:
    ranker = InfoTheory.InfoBitRanker(nBits, nActs, InfoTheory.InfoType.BIASENTROPY)
    ranker.SetBiasList(biasList)
  else:
    ranker = InfoTheory.InfoBitRanker(nBits, nActs, InfoTheory.InfoType.ENTROPY)
  i = 0
  fps = []
  for mol in suppl:
    if not acts:
      try:
        act = int(mol.GetProp(actName))
      except KeyError:
        message('ERROR: Molecule has no property: %s\n' % (actName))
        message('\tAvailable properties are: %s\n' % (str(mol.GetPropNames())))
        raise KeyError(actName)
    else:
      act = acts[i]
    if i and not i % reportFreq:
      if nMols > 0:
        message('Done %d of %d.\n' % (i, nMols))
      else:
        message('Done %d.\n' % (i))
    fp = fpgen.GetFPForMol(mol, catalog)
    ranker.AccumulateVotes(fp, act)
    i += 1
    if collectFps:
      fps.append(fp)
  gains = ranker.GetTopN(topN)
  return gains, fps


def CalcGainsFromFps(suppl, fps, topN=-1, actName='', acts=None, nActs=2, reportFreq=10,
                     biasList=None):
  """ calculates info gains from a set of fingerprints

      *DOC*

    """
  nBits = len(fps[0])
  if topN < 0:
    topN = nBits
  if not actName and not acts:
    actName = suppl[0].GetPropNames()[-1]

  if hasattr(suppl, '__len__'):
    nMols = len(suppl)
  else:
    nMols = -1
  if biasList:
    ranker = InfoTheory.InfoBitRanker(nBits, nActs, InfoTheory.InfoType.BIASENTROPY)
    ranker.SetBiasList(biasList)
  else:
    ranker = InfoTheory.InfoBitRanker(nBits, nActs, InfoTheory.InfoType.ENTROPY)
  for i, mol in enumerate(suppl):
    if not acts:
      try:
        act = int(mol.GetProp(actName))
      except KeyError:
        message('ERROR: Molecule has no property: %s\n' % (actName))
        message('\tAvailable properties are: %s\n' % (str(mol.GetPropNames())))
        raise KeyError(actName)
    else:
      act = acts[i]
    if i and not i % reportFreq:
      if nMols > 0:
        message('Done %d of %d.\n' % (i, nMols))
      else:
        message('Done %d.\n' % (i))
    fp = fps[i]
    ranker.AccumulateVotes(fp, act)
  gains = ranker.GetTopN(topN)
  return gains


def OutputGainsData(outF, gains, cat, nActs=2):
  actHeaders = ['Act-%d' % (x) for x in range(nActs)]
  if cat:
    outF.write('id,Description,Gain,%s\n' % (','.join(actHeaders)))
  else:
    outF.write('id,Gain,%s\n' % (','.join(actHeaders)))
  for entry in gains:
    id_ = int(entry[0])
    outL = [str(id_)]
    if cat:
      descr = cat.GetBitDescription(id_)
      outL.append(descr)
    outL.append('%.6f' % entry[1])
    outL += ['%d' % x for x in entry[2:]]
    outF.write(','.join(outL))
    outF.write('\n')


def ProcessGainsData(inF, delim=',', idCol=0, gainCol=1):
  """ reads a list of ids and info gains out of an input file

    """
  res = []
  _ = inF.readline()
  for line in inF:
    splitL = line.strip().split(delim)
    res.append((splitL[idCol], float(splitL[gainCol])))
  return res


def ShowDetails(catalog, gains, nToDo=-1, outF=sys.stdout, idCol=0, gainCol=1, outDelim=','):
  """
     gains should be a sequence of sequences.  The idCol entry of each
     sub-sequence should be a catalog ID.  _ProcessGainsData()_ provides
     suitable input.

    """
  if nToDo < 0:
    nToDo = len(gains)
  for i in range(nToDo):
    id_ = int(gains[i][idCol])
    gain = float(gains[i][gainCol])
    descr = catalog.GetFragDescription(id_)
    if descr:
      outF.write('%s\n' % (outDelim.join((str(id_), descr, str(gain)))))


def SupplierFromDetails(details):
  from rdkit.VLib.NodeLib.DbMolSupply import DbMolSupplyNode
  from rdkit.VLib.NodeLib.SmilesSupply import SmilesSupplyNode

  if details.dbName:
    conn = DbConnect(details.dbName, details.tableName)
    suppl = DbMolSupplyNode(conn.GetData())
  else:
    suppl = SmilesSupplyNode(details.inFileName, delim=details.delim, nameColumn=details.nameCol,
                             smilesColumn=details.smiCol, titleLine=details.hasTitle)
    if isinstance(details.actCol, int):
      suppl.reset()
      m = next(suppl)
      actName = m.GetPropNames()[details.actCol]
      details.actCol = actName
    if isinstance(details.nameCol, int):
      suppl.reset()
      m = next(suppl)
      nameName = m.GetPropNames()[details.nameCol]
      details.nameCol = nameName
      suppl.reset()
  if isinstance(details.actCol, int):
    suppl.reset()
    m = next(suppl)
    actName = m.GetPropNames()[details.actCol]
    details.actCol = actName
  if isinstance(details.nameCol, int):
    suppl.reset()
    m = next(suppl)
    nameName = m.GetPropNames()[details.nameCol]
    details.nameCol = nameName
    suppl.reset()
  return suppl


def Usage():
  print("This is BuildFragmentCatalog")
  print('usage error')
  # print(__doc__)
  sys.exit(-1)


class RunDetails(object):
  numMols = -1
  doBuild = 0
  doSigs = 0
  doScore = 0
  doGains = 0
  doDetails = 0
  catalogName = None
  onBitsName = None
  scoresName = None
  gainsName = None
  dbName = ''
  tableName = None
  detailsName = None
  inFileName = None
  fpName = None
  minPath = 2
  maxPath = 6
  smiCol = 1
  actCol = -1
  nameCol = -1
  hasTitle = 1
  nActs = 2
  nBits = -1
  delim = ','
  biasList = None
  topN = -1


def ParseArgs(details):
  import getopt
  try:
    args, extras = getopt.getopt(sys.argv[1:], 'n:d:cst', [
      'catalog=', 'onbits=', 'scoresFile=', 'gainsFile=', 'detailsFile=', 'fpFile=', 'minPath=',
      'maxPath=', 'smiCol=', 'actCol=', 'nameCol=', 'nActs=', 'nBits=', 'biasList=', 'topN=',
      'build', 'sigs', 'gains', 'details', 'score', 'noTitle'
    ])
  except Exception:
    sys.stderr.write('Error parsing command line:\n')
    import traceback
    traceback.print_exc()
    Usage()
  for arg, val in args:
    if arg == '-n':
      details.numMols = int(val)
    elif arg == '-c':
      details.delim = ','
    elif arg == '-s':
      details.delim = ' '
    elif arg == '-t':
      details.delim = '\t'
    elif arg == '-d':
      details.dbName = val
    elif arg == '--build':
      details.doBuild = 1
    elif arg == '--score':
      details.doScore = 1
    elif arg == '--gains':
      details.doGains = 1
    elif arg == '--sigs':
      details.doSigs = 1
    elif arg == '-details':
      details.doDetails = 1
    elif arg == '--catalog':
      details.catalogName = val
    elif arg == '--onbits':
      details.onBitsName = val
    elif arg == '--scoresFile':
      details.scoresName = val
    elif arg == '--gainsFile':
      details.gainsName = val
    elif arg == '--detailsFile':
      details.detailsName = val
    elif arg == '--fpFile':
      details.fpName = val
    elif arg == '--minPath':
      details.minPath = int(val)
    elif arg == '--maxPath':
      details.maxPath = int(val)
    elif arg == '--smiCol':
      try:
        details.smiCol = int(val)
      except ValueError:
        details.smiCol = val
    elif arg == '--actCol':
      try:
        details.actCol = int(val)
      except ValueError:
        details.actCol = val
    elif arg == '--nameCol':
      try:
        details.nameCol = int(val)
      except ValueError:
        details.nameCol = val
    elif arg == '--nActs':
      details.nActs = int(val)
    elif arg == '--nBits':
      details.nBits = int(val)
    elif arg == '--noTitle':
      details.hasTitle = 0
    elif arg == '--biasList':
      details.biasList = tuple(eval(val))
    elif arg == '--topN':
      details.topN = int(val)
    elif arg == '-h':
      Usage()
      sys.exit(0)
    else:
      Usage()
  if len(extras):
    if details.dbName:
      details.tableName = extras[0]
    else:
      details.inFileName = extras[0]
  else:
    Usage()


if __name__ == '__main__':
  import time
  details = RunDetails()
  ParseArgs(details)
  from io import StringIO
  suppl = SupplierFromDetails(details)

  cat = None
  obls = None
  if details.doBuild:
    if not suppl:
      message("We require inData to generate a catalog\n")
      sys.exit(-2)
    message("Building catalog\n")
    t1 = time.time()
    cat = BuildCatalog(suppl, maxPts=details.numMols, minPath=details.minPath,
                       maxPath=details.maxPath)
    t2 = time.time()
    message("\tThat took %.2f seconds.\n" % (t2 - t1))
    if details.catalogName:
      message("Dumping catalog data\n")
      pickle.dump(cat, open(details.catalogName, 'wb+'))
  elif details.catalogName:
    message("Loading catalog\n")
    cat = pickle.load(open(details.catalogName, 'rb'))
    if details.onBitsName:
      try:
        obls = pickle.load(open(details.onBitsName, 'rb'))
      except Exception:
        obls = None
      else:
        if len(obls) < (inD.count('\n') - 1):
          obls = None
  scores = None
  if details.doScore:
    if not suppl:
      message("We require inData to score molecules\n")
      sys.exit(-2)
    if not cat:
      message("We require a catalog to score molecules\n")
      sys.exit(-2)
    message("Scoring compounds\n")
    if not obls or len(obls) < details.numMols:
      scores, obls = ScoreMolecules(suppl, cat, maxPts=details.numMols, actName=details.actCol,
                                    nActs=details.nActs)
      if details.scoresName:
        pickle.dump(scores, open(details.scoresName, 'wb+'))
      if details.onBitsName:
        pickle.dump(obls, open(details.onBitsName, 'wb+'))
    else:
      scores = ScoreFromLists(obls, suppl, cat, maxPts=details.numMols, actName=details.actCol,
                              nActs=details.nActs)
  elif details.scoresName:
    scores = pickle.load(open(details.scoresName, 'rb'))

  if details.fpName and os.path.exists(details.fpName) and not details.doSigs:
    message("Reading fingerprints from file.\n")
    fps = pickle.load(open(details.fpName, 'rb'))
  else:
    fps = []
  gains = None
  if details.doGains:
    if not suppl:
      message("We require inData to calculate gains\n")
      sys.exit(-2)
    if not (cat or fps):
      message("We require either a catalog or fingerprints to calculate gains\n")
      sys.exit(-2)
    message("Calculating Gains\n")
    t1 = time.time()
    if details.fpName:
      collectFps = 1
    else:
      collectFps = 0
    if not fps:
      gains, fps = CalcGains(suppl, cat, topN=details.topN, actName=details.actCol,
                             nActs=details.nActs, biasList=details.biasList, collectFps=collectFps)
      if details.fpName:
        message("Writing fingerprint file.\n")
        tmpF = open(details.fpName, 'wb+')
        pickle.dump(fps, tmpF, 1)
        tmpF.close()
    else:
      gains = CalcGainsFromFps(suppl, fps, topN=details.topN, actName=details.actCol,
                               nActs=details.nActs, biasList=details.biasList)
    t2 = time.time()
    message("\tThat took %.2f seconds.\n" % (t2 - t1))
    if details.gainsName:
      outF = open(details.gainsName, 'w+')
      OutputGainsData(outF, gains, cat, nActs=details.nActs)
  else:
    if details.gainsName:
      inF = open(details.gainsName, 'r')
      gains = ProcessGainsData(inF)

  if details.doDetails:
    if not cat:
      message("We require a catalog to get details\n")
      sys.exit(-2)
    if not gains:
      message("We require gains data to get details\n")
      sys.exit(-2)
    io = StringIO()
    io.write('id,SMILES,gain\n')
    ShowDetails(cat, gains, nToDo=details.nBits, outF=io)
    if details.detailsName:
      open(details.detailsName, 'w+').write(io.getvalue())
    else:
      sys.stderr.write(io.getvalue())
