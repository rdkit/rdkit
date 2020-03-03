# $Id$
#
#  Copyright (c) 2007, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met: 
#
#     * Redistributions of source code must retain the above copyright 
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following 
#       disclaimer in the documentation and/or other materials provided 
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
#       nor the names of its contributors may be used to endorse or promote 
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Created by Greg Landrum, July 2007
_version = "0.13.0"

_description = """
  NOTES:

    - the property names for the database are the union of those for
      all molecules.

    - missing property values will be set to 'N/A', though this can be
      changed with the --missingPropertyVal argument.
    
    - The property names may be altered on loading the database.  Any
      non-alphanumeric character in a property name will be replaced
      with '_'. e.g. "Gold.Goldscore.Constraint.Score" becomes
      "Gold_Goldscore_Constraint_Score".  This is important to know
      when querying.

    - Property names are not case sensitive in the database; this may
      cause some problems if they are case sensitive in the sd file.

      
"""
import argparse

from rdkit import RDConfig
from rdkit import Chem
from rdkit.Dbase.DbConnection import DbConnect
from rdkit.Dbase import DbModule
from rdkit.RDLogger import logger
from rdkit.Chem.MolDb import Loader

logger = logger()
import sys, os
import io
import pickle
from rdkit.Chem.MolDb.FingerprintUtils import BuildSigFactory, LayeredOptions
from rdkit.Chem.MolDb import FingerprintUtils

# ---- ---- ---- ----  ---- ---- ---- ----  ---- ---- ---- ----  ---- ---- ---- ----

def initParser():
  """ Initialize the command line parser """ 
  parser = argparse.ArgumentParser(usage='CreateDB [optional arguments] <filename>',
                                   description=_description,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)

  parser.add_argument('filename', nargs='?', help='File containg molecules to load into database')
  parser.add_argument('--version', action='version', version='%(prog)s ' + _version)
  
  parser.add_argument('--outDir', '--dbDir', default='', help='name of the output directory')
  parser.add_argument('--molDbName', default='Compounds.sqlt', help='name of the molecule database')
  parser.add_argument('--molIdName', default='compound_id', help='name of the database key column')
  parser.add_argument('--regName', default='molecules', help='name of the molecular registry table')
  parser.add_argument('--pairDbName', default='AtomPairs.sqlt', help='name of the atom pairs database')
  parser.add_argument('--pairTableName', default='atompairs', help='name of the atom pairs table')
  parser.add_argument('--fpDbName', default='Fingerprints.sqlt',
                    help='name of the 2D fingerprints database')
  parser.add_argument('--fpTableName', default='rdkitfps', help='name of the 2D fingerprints table')
  parser.add_argument('--layeredTableName', default='layeredfps',
                    help='name of the layered fingerprints table')
  parser.add_argument('--descrDbName', default='Descriptors.sqlt',
                    help='name of the descriptor database')
  parser.add_argument('--descrTableName', default='descriptors_v1', help='name of the descriptor table')
  parser.add_argument('--descriptorCalcFilename', default=os.path.join(RDConfig.RDBaseDir, 'Projects',
                                                                     'DbCLI', 'moe_like.dsc'),
                    help='name of the file containing the descriptor calculator')
  parser.add_argument('--errFilename', default='loadErrors.txt',
                    help='name of the file to contain information about molecules that fail to load')
  parser.add_argument('--noPairs', default=True, dest='doPairs', action='store_false',
                    help='skip calculating atom pairs')
  parser.add_argument('--noFingerprints', default=True, dest='doFingerprints', action='store_false',
                    help='skip calculating 2D fingerprints')
  parser.add_argument('--noLayeredFps', default=True, dest='doLayered', action='store_false',
                    help='skip calculating layered fingerprints')
  parser.add_argument('--noDescriptors', default=True, dest='doDescriptors', action='store_false',
                    help='skip calculating descriptors')
  parser.add_argument('--noProps', default=False, dest='skipProps', action='store_true',
                    help="don't include molecular properties in the database")
  parser.add_argument('--noSmiles', default=False, dest='skipSmiles', action='store_true',
                    help="don't include SMILES in the database (can make loading somewhat faster)")
  parser.add_argument('--maxRowsCached', default=-1,
                    help="maximum number of rows to cache before doing a database commit")
  
  parser.add_argument('--silent', default=False, action='store_true',
                    help='do not provide status messages')
  
  parser.add_argument('--molFormat', default='', choices=('smiles', 'sdf', ''),
                    help='specify the format of the input file')
  parser.add_argument(
    '--nameProp', default='_Name',
    help='specify the SD property to be used for the molecule names. Default is to use the mol block name')
  parser.add_argument(
    '--missingPropertyVal', default='N/A',
    help='value to insert in the database if a property value is missing. Default is %(default)s.')
  parser.add_argument('--addProps', default=False, action='store_true',
                    help='add computed properties to the output')
  parser.add_argument('--noExtras', default=False, action='store_true',
                    help='skip all non-molecule databases')
  parser.add_argument('--skipLoad', '--skipMols', action="store_false", dest='loadMols', default=True,
                    help='skip the molecule loading (assumes mol db already exists)')
  parser.add_argument('--updateDb', '--update', default=False, action='store_true',
                    help='add to an existing database')
  parser.add_argument('--doPharm2D', default=False, action='store_true',
                    help='skip calculating Pharm2D fingerprints')
  parser.add_argument('--pharm2DTableName', default='pharm2dfps',
                    help='name of the Pharm2D fingerprints table')
  parser.add_argument('--fdefFile', '--fdef',
                    default=os.path.join(RDConfig.RDDataDir, 'Novartis1.fdef'),
                    help='provide the name of the fdef file to use for 2d pharmacophores')
  parser.add_argument('--doGobbi2D', default=False, action='store_true',
                    help='skip calculating Gobbi 2D fingerprints')
  parser.add_argument('--gobbi2DTableName', default='gobbi2dfps',
                    help='name of the Gobbi 2D fingerprints table')
  
  parser.add_argument('--noMorganFps', '--noCircularFps', default=True, dest='doMorganFps',
                    action='store_false', help='skip calculating Morgan (circular) fingerprints')
  parser.add_argument('--morganFpTableName', default='morganfps',
                    help='name of the Morgan fingerprints table')
  
  parser.add_argument('--delimiter', '--delim', default=' ', help='the delimiter in the input file')
  parser.add_argument('--titleLine', default=False, action='store_true',
                    help='the input file contains a title line')
  parser.add_argument('--smilesColumn', '--smilesCol', default=0, type=int,
                    help='the column index with smiles')
  parser.add_argument('--nameColumn', '--nameCol', default=1, type=int,
                    help='the column index with mol names')
  return parser


def CreateDb(options, dataFilename='', supplier=None):
  if not dataFilename and supplier is None:
    raise ValueError('Please provide either a data filename or a supplier')

  if options.errFilename:
    errFile = open(os.path.join(options.outDir, options.errFilename), 'w+')
  else:
    errFile = None

  if options.noExtras:
    options.doPairs = False
    options.doDescriptors = False
    options.doFingerprints = False
    options.doPharm2D = False
    options.doGobbi2D = False
    options.doLayered = False
    options.doMorganFps = False

  if options.loadMols:
    if supplier is None:
      if not options.molFormat:
        ext = os.path.splitext(dataFilename)[-1].lower()
        if ext == '.sdf':
          options.molFormat = 'sdf'
        elif ext in ('.smi', '.smiles', '.txt', '.csv'):
          options.molFormat = 'smiles'
          if not options.delimiter:
            # guess the delimiter
            import csv
            sniffer = csv.Sniffer()
            dlct = sniffer.sniff(open(dataFilename, 'r').read(2000))
            options.delimiter = dlct.delimiter
            if not options.silent:
              logger.info(
                'Guessing that delimiter is %s. Use --delimiter argument if this is wrong.' %
                repr(options.delimiter))

        if not options.silent:
          logger.info('Guessing that mol format is %s. Use --molFormat argument if this is wrong.' %
                      repr(options.molFormat))
      if options.molFormat == 'smiles':
        if options.delimiter == '\\t':
          options.delimiter = '\t'
        supplier = Chem.SmilesMolSupplier(
          dataFilename, titleLine=options.titleLine, delimiter=options.delimiter,
          smilesColumn=options.smilesColumn, nameColumn=options.nameColumn)
      else:
        supplier = Chem.SDMolSupplier(dataFilename)
    if not options.silent:
      logger.info('Reading molecules and constructing molecular database.')
    Loader.LoadDb(supplier, os.path.join(options.outDir, options.molDbName), errorsTo=errFile,
                  regName=options.regName, nameCol=options.molIdName, skipProps=options.skipProps,
                  defaultVal=options.missingPropertyVal, addComputedProps=options.addProps,
                  uniqNames=True, skipSmiles=options.skipSmiles,
                  maxRowsCached=int(options.maxRowsCached), silent=options.silent,
                  nameProp=options.nameProp, lazySupplier=int(options.maxRowsCached) > 0,
                  startAnew=not options.updateDb)

  if options.doPairs:
    pairConn = DbConnect(os.path.join(options.outDir, options.pairDbName))
    pairCurs = pairConn.GetCursor()
    try:
      pairCurs.execute('drop table %s' % (options.pairTableName))
    except Exception:
      pass
    pairCurs.execute(
      'create table %s (guid integer not null primary key,%s varchar not null unique,atompairfp blob,torsionfp blob)'
      % (options.pairTableName, options.molIdName))

  if options.doFingerprints or options.doPharm2D or options.doGobbi2D or options.doLayered:
    fpConn = DbConnect(os.path.join(options.outDir, options.fpDbName))
    fpCurs = fpConn.GetCursor()
    try:
      fpCurs.execute('drop table %s' % (options.fpTableName))
    except Exception:
      pass
    try:
      fpCurs.execute('drop table %s' % (options.pharm2DTableName))
    except Exception:
      pass
    try:
      fpCurs.execute('drop table %s' % (options.gobbi2DTableName))
    except Exception:
      pass
    try:
      fpCurs.execute('drop table %s' % (options.layeredTableName))
    except Exception:
      pass

    if options.doFingerprints:
      fpCurs.execute(
        'create table %s (guid integer not null primary key,%s varchar not null unique,rdkfp blob)'
        % (options.fpTableName, options.molIdName))
    if options.doLayered:
      layeredQs = ','.join('?' * LayeredOptions.nWords)
      colDefs = ','.join(['Col_%d integer' % (x + 1) for x in range(LayeredOptions.nWords)])
      fpCurs.execute(
        'create table %s (guid integer not null primary key,%s varchar not null unique,%s)' % (
          options.layeredTableName, options.molIdName, colDefs))

    if options.doPharm2D:
      fpCurs.execute(
        'create table %s (guid integer not null primary key,%s varchar not null unique,pharm2dfp blob)'
        % (options.pharm2DTableName, options.molIdName))
      sigFactory = BuildSigFactory(options)
    if options.doGobbi2D:
      fpCurs.execute(
        'create table %s (guid integer not null primary key,%s varchar not null unique,gobbi2dfp blob)'
        % (options.gobbi2DTableName, options.molIdName))
      from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D

  if options.doMorganFps:
    fpConn = DbConnect(os.path.join(options.outDir, options.fpDbName))
    fpCurs = fpConn.GetCursor()
    try:
      fpCurs.execute('drop table %s' % (options.morganFpTableName))
    except Exception:
      pass
    fpCurs.execute(
      'create table %s (guid integer not null primary key,%s varchar not null unique,morganfp blob)'
      % (options.morganFpTableName, options.molIdName))

  if options.doDescriptors:
    descrConn = DbConnect(os.path.join(options.outDir, options.descrDbName))
    with open(options.descriptorCalcFilename, 'r') as inTF:
      buf = inTF.read().replace('\r\n', '\n').encode('utf-8')
      inTF.close()
    calc = pickle.load(io.BytesIO(buf))
    nms = [x for x in calc.GetDescriptorNames()]
    descrCurs = descrConn.GetCursor()
    descrs = ['guid integer not null primary key', '%s varchar not null unique' % options.molIdName]
    descrs.extend(['%s float' % x for x in nms])
    try:
      descrCurs.execute('drop table %s' % (options.descrTableName))
    except Exception:
      pass
    descrCurs.execute('create table %s (%s)' % (options.descrTableName, ','.join(descrs)))
    descrQuery = ','.join([DbModule.placeHolder] * len(descrs))
  pairRows = []
  fpRows = []
  layeredRows = []
  descrRows = []
  pharm2DRows = []
  gobbi2DRows = []
  morganRows = []

  if not options.silent:
    logger.info('Generating fingerprints and descriptors:')
  molConn = DbConnect(os.path.join(options.outDir, options.molDbName))
  molCurs = molConn.GetCursor()
  if not options.skipSmiles:
    molCurs.execute('select guid,%s,smiles,molpkl from %s' % (options.molIdName, options.regName))
  else:
    molCurs.execute('select guid,%s,molpkl from %s' % (options.molIdName, options.regName))
  i = 0
  while 1:
    try:
      tpl = molCurs.fetchone()
      molGuid = tpl[0]
      molId = tpl[1]
      pkl = tpl[-1]
      i += 1
    except Exception:
      break
    if isinstance(pkl, (bytes, str)):
      mol = Chem.Mol(pkl)
    else:
      mol = Chem.Mol(str(pkl))
    if not mol:
      continue

    if options.doPairs:
      pairs = FingerprintUtils.BuildAtomPairFP(mol)
      torsions = FingerprintUtils.BuildTorsionsFP(mol)
      pkl1 = DbModule.binaryHolder(pairs.ToBinary())
      pkl2 = DbModule.binaryHolder(torsions.ToBinary())
      row = (molGuid, molId, pkl1, pkl2)
      pairRows.append(row)
    if options.doFingerprints:
      fp2 = FingerprintUtils.BuildRDKitFP(mol)
      pkl = DbModule.binaryHolder(fp2.ToBinary())
      row = (molGuid, molId, pkl)
      fpRows.append(row)
    if options.doLayered:
      words = LayeredOptions.GetWords(mol)
      row = [molGuid, molId] + words
      layeredRows.append(row)
    if options.doDescriptors:
      descrs = calc.CalcDescriptors(mol)
      row = [molGuid, molId]
      row.extend(descrs)
      descrRows.append(row)
    if options.doPharm2D:
      FingerprintUtils.sigFactory = sigFactory
      fp = FingerprintUtils.BuildPharm2DFP(mol)
      pkl = DbModule.binaryHolder(fp.ToBinary())
      row = (molGuid, molId, pkl)
      pharm2DRows.append(row)
    if options.doGobbi2D:
      FingerprintUtils.sigFactory = Gobbi_Pharm2D.factory
      fp = FingerprintUtils.BuildPharm2DFP(mol)
      pkl = DbModule.binaryHolder(fp.ToBinary())
      row = (molGuid, molId, pkl)
      gobbi2DRows.append(row)
    if options.doMorganFps:
      morgan = FingerprintUtils.BuildMorganFP(mol)
      pkl = DbModule.binaryHolder(morgan.ToBinary())
      row = (molGuid, molId, pkl)
      morganRows.append(row)

    if not i % 500:
      if len(pairRows):
        pairCurs.executemany('insert into %s values (?,?,?,?)' % options.pairTableName, pairRows)
        pairRows = []
        pairConn.Commit()
      if len(fpRows):
        fpCurs.executemany('insert into %s values (?,?,?)' % options.fpTableName, fpRows)
        fpRows = []
        fpConn.Commit()
      if len(layeredRows):
        fpCurs.executemany('insert into %s values (?,?,%s)' % (options.layeredTableName, layeredQs),
                           layeredRows)
        layeredRows = []
        fpConn.Commit()
      if len(descrRows):
        descrCurs.executemany('insert into %s values (%s)' % (options.descrTableName, descrQuery),
                              descrRows)
        descrRows = []
        descrConn.Commit()
      if len(pharm2DRows):
        fpCurs.executemany('insert into %s values (?,?,?)' % options.pharm2DTableName, pharm2DRows)
        pharm2DRows = []
        fpConn.Commit()
      if len(gobbi2DRows):
        fpCurs.executemany('insert into %s values (?,?,?)' % options.gobbi2DTableName, gobbi2DRows)
        gobbi2DRows = []
        fpConn.Commit()
      if len(morganRows):
        fpCurs.executemany('insert into %s values (?,?,?)' % options.morganFpTableName, morganRows)
        morganRows = []
        fpConn.Commit()

    if not options.silent and not i % 500:
      logger.info('  Done: %d' % (i))

  if len(pairRows):
    pairCurs.executemany('insert into %s values (?,?,?,?)' % options.pairTableName, pairRows)
    pairRows = []
    pairConn.Commit()
  if len(fpRows):
    fpCurs.executemany('insert into %s values (?,?,?)' % options.fpTableName, fpRows)
    fpRows = []
    fpConn.Commit()
  if len(layeredRows):
    fpCurs.executemany('insert into %s values (?,?,%s)' % (options.layeredTableName, layeredQs),
                       layeredRows)
    layeredRows = []
    fpConn.Commit()
  if len(descrRows):
    descrCurs.executemany('insert into %s values (%s)' % (options.descrTableName, descrQuery),
                          descrRows)
    descrRows = []
    descrConn.Commit()
  if len(pharm2DRows):
    fpCurs.executemany('insert into %s values (?,?,?)' % options.pharm2DTableName, pharm2DRows)
    pharm2DRows = []
    fpConn.Commit()
  if len(gobbi2DRows):
    fpCurs.executemany('insert into %s values (?,?,?)' % options.gobbi2DTableName, gobbi2DRows)
    gobbi2DRows = []
    fpConn.Commit()
  if len(morganRows):
    fpCurs.executemany('insert into %s values (?,?,?)' % options.morganFpTableName, morganRows)
    morganRows = []
    fpConn.Commit()

  if not options.silent:
    logger.info('Finished.')


if __name__ == '__main__':
  parser = initParser()
  options = parser.parse_args()
  if options.loadMols:
    if options.filename is None:
      parser.error('please provide a filename argument')
    dataFilename = options.filename
    try:
      dataFile = open(dataFilename, 'r')
    except IOError:
      logger.error('input file %s does not exist' % (dataFilename))
      sys.exit(0)
    dataFile = None

  if not options.outDir:
    prefix = os.path.splitext(dataFilename)[0]
    options.outDir = prefix

  if not os.path.exists(options.outDir):
    try:
      os.mkdir(options.outDir)
    except Exception:
      logger.error('could not create output directory %s' % options.outDir)
      sys.exit(1)

  if 1:
    CreateDb(options, dataFilename)
  else:
    import cProfile
    cProfile.run("CreateDb(options,dataFilename)", "create.prof")
    import pstats
    p = pstats.Stats('create.prof')
    p.strip_dirs().sort_stats('cumulative').print_stats(25)
