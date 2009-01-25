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
_version = "0.7.3"
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
from pyRDKit import RDConfig
from pyRDKit import Chem
from pyRDKit.Dbase.DbConnection import DbConnect
from pyRDKit.Dbase import DbModule
from pyRDKit.RDLogger import logger
from pyRDKit.Chem.MolDb import Loader

logger = logger()
import cPickle,sys,os

# ---- ---- ---- ----  ---- ---- ---- ----  ---- ---- ---- ----  ---- ---- ---- ---- 
from optparse import OptionParser
parser=OptionParser(_usage,version='%prog '+_version)
parser.add_option('--outDir','--dbDir',default='.',
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
parser.add_option('--noOldFingerprints',default=True,dest='doOldFps',action='store_false',
                  help='skip calculating 2D fingerprints using the old fingerprinter')
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

parser.add_option('--molFormat',default='smiles',choices=('smiles','sdf'),
                  help='specify the format of the input file')
parser.add_option('--nameProp',default='_Name',
                  help='specify the SD property to be used for the molecule names. Default is to use the mol block name')
parser.add_option('--missingPropertyVal',default='N/A',
                  help='value to insert in the database if a property value is missing. Default is %default.')
parser.add_option('--addProps',default=False,action='store_true',
                  help='add computed properties to the output')
parser.add_option('--noExtras',default=False,action='store_true',
                  help='skip all non-molecule databases')

parser.add_option('--delimiter','--delim',default=' ',
                  help='the delimiter in the input file')
parser.add_option('--titleLine',default=False,action='store_true',
                  help='the input file contains a title line')

if __name__=='__main__':
  options,args = parser.parse_args()
  if len(args)!=1:
    parser.error('please provide a filename argument')

  if not os.path.exists(options.outDir):
    try:
      os.mkdir(options.outDir)
    except: 
      logger.error('could not create output directory %s'%options.outDir)
      sys.exit(1)
  molConn = DbConnect(os.path.join(options.outDir,options.molDbName))
  dataFilename = args[0]
  dataFile = file(dataFilename,'r')
  dataFile=None
  errFile=file(os.path.join(options.outDir,options.errFilename),'w+')
  
  if options.molFormat=='smiles':
    supplier=Chem.SmilesMolSupplier(dataFilename,
                                    titleLine=options.titleLine,
                                    delimiter=options.delimiter)
  else:
    supplier = Chem.SDMolSupplier(dataFilename)

  if options.noExtras:
    options.doPairs=False
    options.doDescriptors=False
    options.doFingerprints=False

  if not options.silent: logger.info('Reading molecules and constructing molecular database.')
  Loader.LoadDb(supplier,os.path.join(options.outDir,options.molDbName),
                errorsTo=errFile,regName=options.regName,nameCol=options.molIdName,
                skipProps=options.skipProps,defaultVal=options.missingPropertyVal,
                addComputedProps=options.addProps,uniqNames=True,
                skipSmiles=options.skipSmiles,maxRowsCached=int(options.maxRowsCached),
                silent=options.silent,nameProp=options.nameProp,lazySupplier=False)
  if options.doPairs:
    from pyRDKit.Chem.AtomPairs import Pairs,Torsions
    pairConn = DbConnect(os.path.join(options.outDir,options.pairDbName))
    pairCurs = pairConn.GetCursor()
    try:
      pairCurs.execute('drop table %s'%(options.pairTableName))
    except:
      pass
    pairCurs.execute('create table %s (%s varchar not null primary key,atompairfp blob,torsionfp blob)'%(options.pairTableName,
                                                                                                         options.molIdName))

  if options.doFingerprints:
    fpConn = DbConnect(os.path.join(options.outDir,options.fpDbName))
    fpCurs=fpConn.GetCursor()
    try:
      fpCurs.execute('drop table %s'%(options.fpTableName))
    except:
      pass
    fpCurs.execute('create table %s (%s varchar not null primary key,autofragmentfp blob,rdkfp blob)'%(options.fpTableName,
                                                                                                       options.molIdName))
    from pyRDKit.Chem.Fingerprints import FingerprintMols
    details = FingerprintMols.FingerprinterDetails()
    fpArgs = details.__dict__
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
  descrRows = []

  if not options.silent: logger.info('Generating fingerprints and descriptors:')
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
      pairs = Pairs.GetAtomPairFingerprintAsIntVect(mol)
      torsions = Torsions.GetTopologicalTorsionFingerprintAsIntVect(mol)
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
      if options.doOldFps:
        fp = FingerprintMols.FingerprintMol(mol,**fpArgs)
        pkl1 = DbModule.binaryHolder(fp.ToBinary())
      else:
        pkl1 = ''
      fp2 = Chem.RDKFingerprint(mol)
      pkl2 = DbModule.binaryHolder(fp2.ToBinary())
      row = [id,pkl1,pkl2]
      fpRows.append(row)
      if len(fpRows)>=500:
        fpCurs.executemany('insert into %s values (?,?,?)'%options.fpTableName,
                           fpRows)
        fpRows = []
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

    if not options.silent and not i%500: 
      logger.info('  Done: %d'%(i))

  if len(pairRows):
    pairCurs.executemany('insert into %s values (?,?,?)'%options.pairTableName,
                         pairRows)
    pairRows = []
    pairConn.Commit()
  if len(fpRows):
    fpCurs.executemany('insert into %s values (?,?,?)'%options.fpTableName,
                       fpRows)
    fpRows = []
    fpConn.Commit()
  if len(descrRows):
    descrCurs.executemany('insert into %s values (%s)'%(options.descrTableName,descrQuery),
                          descrRows)
    descrRows = []
    descrConn.Commit()
    
                




