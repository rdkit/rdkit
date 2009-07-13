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
_version = "0.11.0"

_usage="""
 CreateDb [optional arguments] <filename>

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
from rdkit import RDConfig
from rdkit import Chem
from rdkit.Dbase.DbConnection import DbConnect
from rdkit.Dbase import DbModule
from rdkit.RDLogger import logger
from rdkit.Chem.MolDb import Loader

from rdkit.Chem.Descriptors import rdMolDescriptors
from rdkit import DataStructs

logger = logger()
import cPickle,sys,os

from rdkit.Chem.MolDb.FingerprintUtils import BuildSigFactory,LayeredOptions
from rdkit.Chem.MolDb import FingerprintUtils

# ---- ---- ---- ----  ---- ---- ---- ----  ---- ---- ---- ----  ---- ---- ---- ---- 
from optparse import OptionParser
parser=OptionParser(_usage,version='%prog '+_version)
parser.add_option('--outDir','--dbDir',default='',
                  help='name of the output directory')
parser.add_option('--molDbName',default='Compounds.sqlt',
                  help='name of the molecule database')
parser.add_option('--molIdName',default='compound_id',
                  help='name of the database key column')
parser.add_option('--regName',default='molecules',
                  help='name of the molecular registry table')
parser.add_option('--pairDbName',default='AtomPairs.sqlt',
                  help='name of the atom pairs database')
parser.add_option('--pairTableName',default='atompairs',
                  help='name of the atom pairs table')
parser.add_option('--fpDbName',default='Fingerprints.sqlt',
                  help='name of the 2D fingerprints database')
parser.add_option('--fpTableName',default='rdkitfps',
                  help='name of the 2D fingerprints table')
parser.add_option('--layeredTableName',default='layeredfps',
                  help='name of the layered fingerprints table')
parser.add_option('--descrDbName',default='Descriptors.sqlt',
                  help='name of the descriptor database')
parser.add_option('--descrTableName',default='descriptors_v1',
                  help='name of the descriptor table')
parser.add_option('--descriptorCalcFilename',default=os.path.join(RDConfig.RDBaseDir,'Projects',
                                                                  'DbCLI','moe_like.dsc'),
                  help='name of the file containing the descriptor calculator')
parser.add_option('--errFilename',default='loadErrors.txt',
                  help='name of the file to contain information about molecules that fail to load')
parser.add_option('--noPairs',default=True,dest='doPairs',action='store_false',
                  help='skip calculating atom pairs')
parser.add_option('--noFingerprints',default=True,dest='doFingerprints',action='store_false',
                  help='skip calculating 2D fingerprints')
parser.add_option('--noLayeredFps',default=True,dest='doLayered',action='store_false',
                  help='skip calculating layered fingerprints')
parser.add_option('--noDescriptors',default=True,dest='doDescriptors',action='store_false',
                  help='skip calculating descriptors')
parser.add_option('--noProps',default=False,dest='skipProps',action='store_true',
                  help="don't include molecular properties in the database")
parser.add_option('--noSmiles',default=False,dest='skipSmiles',action='store_true',
                  help="don't include SMILES in the database (can make loading somewhat faster)")
parser.add_option('--maxRowsCached',default=-1,
                  help="maximum number of rows to cache before doing a database commit")

parser.add_option('--silent',default=False,action='store_true',
                  help='do not provide status messages')

parser.add_option('--molFormat',default='',choices=('smiles','sdf',''),
                  help='specify the format of the input file')
parser.add_option('--nameProp',default='_Name',
                  help='specify the SD property to be used for the molecule names. Default is to use the mol block name')
parser.add_option('--missingPropertyVal',default='N/A',
                  help='value to insert in the database if a property value is missing. Default is %default.')
parser.add_option('--addProps',default=False,action='store_true',
                  help='add computed properties to the output')
parser.add_option('--noExtras',default=False,action='store_true',
                  help='skip all non-molecule databases')
parser.add_option('--skipLoad','--skipMols',action="store_false",dest='loadMols',default=True,
                  help='skip the molecule loading (assumes mol db already exists)')

parser.add_option('--doPharm2D',default=False,
                  action='store_true',
                  help='skip calculating Pharm2D fingerprints')
parser.add_option('--pharm2DTableName',default='pharm2dfps',
                  help='name of the Pharm2D fingerprints table')
parser.add_option('--fdefFile','--fdef',
                  default=os.path.join(RDConfig.RDDataDir,'Novartis1.fdef'),
                  help='provide the name of the fdef file to use for 2d pharmacophores')
parser.add_option('--doGobbi2D',default=False,
                  action='store_true',
                  help='skip calculating Gobbi 2D fingerprints')
parser.add_option('--gobbi2DTableName',default='gobbi2dfps',
                  help='name of the Gobbi 2D fingerprints table')

parser.add_option('--noMorganFps','--noCircularFps',default=True,dest='doMorganFps',action='store_false',
                  help='skip calculating Morgan (circular) fingerprints')
parser.add_option('--morganFpTableName',default='morganfps',
                  help='name of the Morgan fingerprints table')

parser.add_option('--delimiter','--delim',default=' ',
                  help='the delimiter in the input file')
parser.add_option('--titleLine',default=False,action='store_true',
                  help='the input file contains a title line')
parser.add_option('--smilesColumn','--smilesCol',default=0,type='int',
                  help='the column index with smiles')
parser.add_option('--nameColumn','--nameCol',default=1,type='int',
                  help='the column index with mol names')



def run():
  options,args = parser.parse_args()
  if options.loadMols:
    if len(args)!=1:
      parser.error('please provide a filename argument')
    dataFilename = args[0]
    try:
      dataFile = file(dataFilename,'r')
    except IOError:
      logger.error('input file %s does not exist'%(dataFilename))
      sys.exit(0)
    dataFile=None

  if not options.outDir:
    prefix = os.path.splitext(dataFilename)[0]
    options.outDir=prefix

  if not os.path.exists(options.outDir):
    try:
      os.mkdir(options.outDir)
    except: 
      logger.error('could not create output directory %s'%options.outDir)
      sys.exit(1)
  errFile=file(os.path.join(options.outDir,options.errFilename),'w+')

  if options.noExtras:
    options.doPairs=False
    options.doDescriptors=False
    options.doFingerprints=False
    options.doPharm2D=False
    options.doGobbi2D=False
    options.doLayered=False
    options.doMorganFps=False

  if options.loadMols:
    if not options.molFormat:
      ext = os.path.splitext(dataFilename)[-1].lower()
      if ext=='.sdf':
        options.molFormat='sdf'
      elif ext in ('.smi','.smiles','.txt','.csv'):
        options.molFormat='smiles'
        if not options.delimiter:
          # guess the delimiter
          import csv
          sniffer = csv.Sniffer()
          dlct=sniffer.sniff(file(dataFilename,'r').read(2000))
          options.delimiter=dlct.delimiter
          if not options.silent:
            logger.info('Guessing that delimiter is %s. Use --delimiter argument if this is wrong.'%repr(options.delimiter))

      if not options.silent:
        logger.info('Guessing that mol format is %s. Use --molFormat argument if this is wrong.'%repr(options.molFormat))  
    if options.molFormat=='smiles':
      if options.delimiter=='\\t': options.delimiter='\t'
      supplier=Chem.SmilesMolSupplier(dataFilename,
                                      titleLine=options.titleLine,
                                      delimiter=options.delimiter,
                                      smilesColumn=options.smilesColumn,
                                      nameColumn=options.nameColumn
                                      )
    else:
      supplier = Chem.SDMolSupplier(dataFilename)


    if not options.silent: logger.info('Reading molecules and constructing molecular database.')
    Loader.LoadDb(supplier,os.path.join(options.outDir,options.molDbName),
                  errorsTo=errFile,regName=options.regName,nameCol=options.molIdName,
                  skipProps=options.skipProps,defaultVal=options.missingPropertyVal,
                  addComputedProps=options.addProps,uniqNames=True,
                  skipSmiles=options.skipSmiles,maxRowsCached=int(options.maxRowsCached),
                  silent=options.silent,nameProp=options.nameProp,
                  lazySupplier=int(options.maxRowsCached)>0)
  if options.doPairs:
    from rdkit.Chem.AtomPairs import Pairs,Torsions
    pairConn = DbConnect(os.path.join(options.outDir,options.pairDbName))
    pairCurs = pairConn.GetCursor()
    try:
      pairCurs.execute('drop table %s'%(options.pairTableName))
    except:
      pass
    pairCurs.execute('create table %s (%s varchar not null primary key,atompairfp blob,torsionfp blob)'%(options.pairTableName,
                                                                                                         options.molIdName))

  if options.doFingerprints or options.doPharm2D or options.doGobbi2D or options.doLayered:
    fpConn = DbConnect(os.path.join(options.outDir,options.fpDbName))
    fpCurs=fpConn.GetCursor()
    try:
      fpCurs.execute('drop table %s'%(options.fpTableName))
    except:
      pass
    try:
      fpCurs.execute('drop table %s'%(options.pharm2DTableName))
    except:
      pass
    try:
      fpCurs.execute('drop table %s'%(options.gobbi2DTableName))
    except:
      pass
    try:
      fpCurs.execute('drop table %s'%(options.layeredTableName))
    except:
      pass

    if options.doFingerprints:
      fpCurs.execute('create table %s (%s varchar not null primary key,rdkfp blob)'%(options.fpTableName,
                                                                                     options.molIdName))
    if options.doLayered:
      colDefs=','.join(['Col_%d integer'%(x+1) for x in range(LayeredOptions.nWords)])
      fpCurs.execute('create table %s (%s varchar not null primary key,%s)'%(options.layeredTableName,
                                                                             options.molIdName,
                                                                             colDefs))
      
    if options.doPharm2D:
      fpCurs.execute('create table %s (%s varchar not null primary key,pharm2dfp blob)'%(options.pharm2DTableName,
                                                                                     options.molIdName))
      sigFactory = BuildSigFactory(options)
    if options.doGobbi2D:
      fpCurs.execute('create table %s (%s varchar not null primary key,gobbi2dfp blob)'%(options.gobbi2DTableName,
                                                                                     options.molIdName))
      from rdkit.Chem.Pharm2D import Generate,Gobbi_Pharm2D

  if options.doMorganFps :
    fpConn = DbConnect(os.path.join(options.outDir,options.fpDbName))
    fpCurs=fpConn.GetCursor()
    try:
      fpCurs.execute('drop table %s'%(options.morganFpTableName))
    except:
      pass
    fpCurs.execute('create table %s (%s varchar not null primary key,morganfp blob)'%(options.morganFpTableName,
                                                                                        options.molIdName))
  morganRows = []

  if options.doDescriptors:
    descrConn=DbConnect(os.path.join(options.outDir,options.descrDbName))
    calc = cPickle.load(file(options.descriptorCalcFilename,'rb'))
    nms = [x for x in calc.GetDescriptorNames()]
    descrCurs = descrConn.GetCursor()
    descrs = ['%s varchar not null primary key'%options.molIdName]
    descrs.extend(['%s float'%x for x in nms])
    try:
      descrCurs.execute('drop table %s'%(options.descrTableName))
    except:
      pass
    descrCurs.execute('create table %s (%s)'%(options.descrTableName,','.join(descrs)))
    descrQuery=','.join([DbModule.placeHolder]*len(descrs))
  pairRows = []
  fpRows = []
  layeredRows = []
  descrRows = []
  pharm2DRows=[]
  gobbi2DRows=[]

  if not options.silent: logger.info('Generating fingerprints and descriptors:')
  molConn = DbConnect(os.path.join(options.outDir,options.molDbName))
  molCurs = molConn.GetCursor()
  if not options.skipSmiles:
    molCurs.execute('select %s,smiles,molpkl from %s'%(options.molIdName,options.regName))
  else:
    molCurs.execute('select %s,molpkl from %s'%(options.molIdName,options.regName))
  i=0

  while 1:
    try:
      tpl = molCurs.fetchone()
      id = tpl[0]
      pkl = tpl[-1]
      i+=1
    except:
      break
    mol = Chem.Mol(str(pkl))
    if not mol: continue
     
    if options.doPairs:
      pairs = FingerprintUtils.BuildAtomPairFP(mol)
      torsions = FingerprintUtils.BuildTorsionsFP(mol)
      pkl1 = DbModule.binaryHolder(pairs.ToBinary())
      pkl2 = DbModule.binaryHolder(torsions.ToBinary())
      row = [id,pkl1,pkl2]
      pairRows.append(row)
      if len(pairRows)>=500:
        pairCurs.executemany('insert into %s values (?,?,?)'%options.pairTableName,
                             pairRows)
        pairRows = []
        pairConn.Commit()
  
    if options.doFingerprints:
      fp2 = FingerprintUtils.BuildRDKitFP(mol)
      pkl = DbModule.binaryHolder(fp2.ToBinary())
      row = [id,pkl]
      fpRows.append(row)
      if len(fpRows)>=500:
        fpCurs.executemany('insert into %s values (?,?)'%options.fpTableName,
                           fpRows)
        fpRows = []
        fpConn.Commit()

    if options.doLayered:
      words = LayeredOptions.GetWords(mol)
      row = [id]+words
      layeredRows.append(row)
      qs = ','.join('?'*LayeredOptions.nWords)
      if len(fpRows)>=500:
        fpCurs.executemany('insert into %s values (?,%s)'%(options.layeredTableName,qs),
                           layeredRows)
        layeredRows = []
        fpConn.Commit()

    if options.doDescriptors:
      descrs= calc.CalcDescriptors(mol)
      row = [id]
      row.extend(descrs)
      descrRows.append(row)
      if len(descrRows)>=500:
        descrCurs.executemany('insert into %s values (%s)'%(options.descrTableName,descrQuery),
                              descrRows)
        descrRows = []
        descrConn.Commit()

    if options.doPharm2D:
      FingerprintUtils.sigFactory=sigFactory
      fp= FingerprintUtils.BuildPharm2DFP(mol)
      pkl = DbModule.binaryHolder(fp.ToBinary())
      row = (id,pkl)
      pharm2DRows.append(row)
      if len(pharm2DRows)>=500:
        fpCurs.executemany('insert into %s values (?,?)'%options.pharm2DTableName,
                           pharm2DRows)
        pharm2DRows = []
        fpConn.Commit()
    if options.doGobbi2D:
      FingerprintUtils.sigFactory=Gobbi_Pharm2D.factory
      fp= FingerprintUtils.BuildPharm2DFP(mol)
      pkl = DbModule.binaryHolder(fp.ToBinary())
      row = (id,pkl)
      gobbi2DRows.append(row)
      if len(gobbi2DRows)>=500:
        fpCurs.executemany('insert into %s values (?,?)'%options.gobbi2DTableName,
                           gobbi2DRows)
        gobbi2DRows = []
        fpConn.Commit()

    if options.doMorganFps:
      morgan = FingerprintUtils.BuildMorganFP(mol)
      pkl1 = DbModule.binaryHolder(morgan.ToBinary())
      row = [id,pkl1]
      morganRows.append(row)
      if len(morganRows)>=500:
        fpCurs.executemany('insert into %s values (?,?)'%options.morganFpTableName,
                             morganRows)
        morganRows = []
        fpConn.Commit()
        
    if not options.silent and not i%500: 
      logger.info('  Done: %d'%(i))

  if len(pairRows):
    pairCurs.executemany('insert into %s values (?,?,?)'%options.pairTableName,
                         pairRows)
    pairRows = []
    pairConn.Commit()
  if len(fpRows):
    fpCurs.executemany('insert into %s values (?,?)'%options.fpTableName,
                       fpRows)
    fpRows = []
    fpConn.Commit()
  if len(layeredRows):
    qs = ','.join('?'*LayeredOptions.nWords)
    fpCurs.executemany('insert into %s values (?,%s)'%(options.layeredTableName,qs),
                       layeredRows)
    layeredRows = []
    fpConn.Commit()
  if len(descrRows):
    descrCurs.executemany('insert into %s values (%s)'%(options.descrTableName,descrQuery),
                          descrRows)
    descrRows = []
    descrConn.Commit()
  if len(pharm2DRows):
    fpCurs.executemany('insert into %s values (?,?)'%options.pharm2DTableName,
                       pharm2DRows)
    pharm2DRows = []
    fpConn.Commit()
  if len(gobbi2DRows):
    fpCurs.executemany('insert into %s values (?,?)'%options.gobbi2DTableName,
                       gobbi2DRows)
    gobbi2DRows = []
    fpConn.Commit()

  if len(morganRows):
    fpCurs.executemany('insert into %s values (?,?)'%options.morganFpTableName,
                       morganRows)
    morganRows = []
    fpConn.Commit()
    
  if not options.silent:
    logger.info('Finished.')

if __name__=='__main__':
  if 0:
    import cProfile
    cProfile.run('run()','fooprof')
    import pstats
    p = pstats.Stats('fooprof')
    p.sort_stats('cumulative').print_stats(25)
  else:
    run()
