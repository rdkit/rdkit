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
#  Created by Greg Landrum, July 2007
#
_version = "0.7.2"
_usage="""
 SearchDb [optional arguments] <sdfilename>

     The sd filename argument can be either an SD file or an MDL mol 
     file.
     

  NOTES:

    - The property names may have been altered on loading the
      database.  Any non-alphanumeric character in a property name
      will be replaced with '_'. e.g."Gold.Goldscore.Constraint.Score" becomes
      "Gold_Goldscore_Constraint_Score".

    - Property names are not case sensitive in the database.

 """
import RDConfig
from Dbase.DbConnection import DbConnect

from RDLogger import logger
logger=logger()
import cPickle,sets


# ----------------------------------------
# ATOM PAIRS
import DataStructs
from DataStructs import SparseIntVect
def DepickleIntVectFP(pkl):
  try:
    fp = cPickle.loads(str(pkl))
  except:
    fp = DataStructs.IntSparseIntVect(str(pkl))
  fp._sumCache = fp.GetTotalVal()
  return fp
def DepickleLongIntVectFP(pkl):
  try:
    fp = cPickle.loads(str(pkl))
  except:
    fp = DataStructs.LongSparseIntVect(str(pkl))
  fp._sumCache = fp.GetTotalVal()
  return fp
def BuildAtomPairFP(mol):
  from Chem.AtomPairs import Pairs
  fp=Pairs.GetAtomPairFingerprintAsIntVect(mol)
  fp._sumCache = fp.GetTotalVal()
  return fp
def BuildTorsionsFP(mol):
  from Chem.AtomPairs import Torsions
  fp=Torsions.GetTopologicalTorsionFingerprintAsIntVect(mol)
  fp._sumCache = fp.GetTotalVal()
  return fp

# ----------------------------------------
# RDKit topological fingerprints:
import DataStructs
def DepickleRDKitFP(pkl):
  fp = DataStructs.ExplicitBitVect(str(pkl))
  return fp
def BuildRDKitFP(mol):
  fp=Chem.RDKFingerprint(mol)
  return fp
def Depickle2DFP(pkl):
  fp = DataStructs.ExplicitBitVect(str(pkl))
  return fp
def Build2DFP(mol):
  from Chem.Fingerprints.FingerprintMols import FingerprintMol
  fp=FingerprintMol(mol)
  return fp

def GetNeighborLists(probes,topN,cursor,
                     simMetric=DataStructs.DiceSimilarity,
                     fpDepickler=DepickleIntVectFP,
                     silent=False):
  probeFps = [x[1] for x in probes]
  validProbes = [x for x in range(len(probeFps)) if probeFps[x] is not None]
  validFps=[probeFps[x] for x in validProbes]
  from DataStructs.TopNContainer import TopNContainer
  nbrLists = [TopNContainer(topN) for x in range(len(probeFps))]

  row = curs.fetchone()
  nDone=0
  while row:
    nDone+=1
    if not silent and not nDone%1000: logger.info('  searched %d rows'%nDone)
    nm,pkl = row
    fp=fpDepickler(pkl)
    if(simMetric==DataStructs.DiceSimilarity):
      scores = DataStructs.BulkDiceSimilarity(fp,validFps)
      for i,score in enumerate(scores):
        nbrLists[validProbes[i]].Insert(score,nm)
    else:
      for i in range(len(probeFps)):
        pfp = probeFps[i]
        if pfp is not None:
          if simMetric==DataStructs.DiceSimilarity:
            score = simMetric(probeFps[i],fp,
                              bounds=nbrLists[i].best[0])
          else:
            score = simMetric(probeFps[i],fp)

          nbrLists[i].Insert(score,nm)
    row = curs.fetchone()
  return nbrLists

def GetMolsFromSmilesFile(dataFilename,errFile,nameProp):
  dataFile=file(dataFilename,'r')
  for idx,line in enumerate(dataFile):
    try:
      smi,nm = line.strip().split(' ')
    except:
      continue
    try:
      m = Chem.MolFromSmiles(smi)
    except:
      m=None
    if not m and errFile:
      print >>errFile,idx,nm,smi
      continue
    yield (nm,smi,m)

def GetMolsFromSDFile(dataFilename,errFile,nameProp):
  suppl = Chem.SDMolSupplier(dataFilename)

  for idx,m in enumerate(suppl):
    if not m and errFile:
      if hasattr(suppl,'GetItemText'):
        d = suppl.GetItemText(idx)
        errFile.write(d)
      else:
        logger.warning('full error file support not complete')
      continue
    smi = Chem.MolToSmiles(m,True)
    if m.HasProp(nameProp):
      nm = m.GetProp(nameProp)
      if not nm:
        logger.warning('molecule found with empty name property')
    else:
      nm = 'Mol_%d'%(idx+1)
    yield nm,smi,m


# ---- ---- ---- ----  ---- ---- ---- ----  ---- ---- ---- ----  ---- ---- ---- ---- 
import os
from optparse import OptionParser
parser=OptionParser(_usage,version='%prog '+_version)
parser.add_option('--dbDir',default='.',
                  help='name of the directory containing the database information. The default is %default')
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
parser.add_option('--pairColName',default='atompairfp',
                  help='name of the atom pair column')
parser.add_option('--torsionsDbName',default='AtomPairs.sqlt',
                  help='name of the topological torsions database (usually the same as the atom pairs database)')
parser.add_option('--torsionsTableName',default='atompairs',
                  help='name of the topological torsions table (usually the same as the atom pairs table)')
parser.add_option('--torsionsColName',default='torsionfp',
                  help='name of the atom pair column')
parser.add_option('--fpDbName',default='Fingerprints.sqlt',
                  help='name of the 2D fingerprints database')
parser.add_option('--fpTableName',default='rdkitfps',
                  help='name of the 2D fingerprints table')
parser.add_option('--fpColName',default='',
                  help='name of the 2D fingerprint column, a sensible default is used')
parser.add_option('--descrDbName',default='Descriptors.sqlt',
                  help='name of the descriptor database')
parser.add_option('--descrTableName',default='descriptors_v1',
                  help='name of the descriptor table')
parser.add_option('--descriptorCalcFilename',default=os.path.join(RDConfig.RDBaseDir,'Projects',
                                                                  'QuickMolDB','moe_like.dsc'),
                  help='name of the file containing the descriptor calculator')
parser.add_option('--similarityType',default='RDK',choices=['RDK','2D','AtomPairs','TopologicalTorsions'],
                  help='Choose the type of similarity to use, possible values: 2D, AtomPairs, TopologicalTorsions. The default is %default')
parser.add_option('--outputDelim',default=',',
                  help='the delimiter for the output file. The default is %default')
parser.add_option('--topN',default=20,type='int',
                  help='the number of neighbors to keep for each query compound. The default is %default')

parser.add_option('--outF','--outFile',default='-',
                  help='The name of the output file. The default is the console (stdout).')

parser.add_option('--transpose',default=False,action="store_true",
                  help='print the results out in a transposed form: e.g. neighbors in rows and probe compounds in columns')

parser.add_option('--molFormat',default='sdf',choices=('smiles','sdf'),
                  help='specify the format of the input file')
parser.add_option('--nameProp',default='_Name',
                  help='specify the SD property to be used for the molecule names. Default is to use the mol block name')

parser.add_option('--smartsQuery','--smarts','--sma',default='',
                  help='provide a SMARTS to be used as a substructure query')
parser.add_option('--negateQuery','--negate',default=False,action='store_true',
                  help='negate the results of the smarts query.')
parser.add_option('--propQuery','--query','-q',default='',
                  help='provide a property query (see the NOTE about property names)')

parser.add_option('--sdfOut','--sdOut',default='',
                  help='export an SD file with the matching molecules')
parser.add_option('--smilesOut','--smiOut',default='',
                  help='export a smiles file with the matching molecules')
parser.add_option('--nonchiralSmiles',dest='chiralSmiles',default=True,action='store_false',
                  help='do not use chiral SMILES in the output')
parser.add_option('--silent',default=False,action='store_true',
                  help='Do not generate status messages.')

if __name__=='__main__':
  import sys,getopt,time
  import Chem
  
  options,args = parser.parse_args()
  if len(args)!=1 and not (options.smartsQuery or options.propQuery):
    parser.error('please either provide a query filename argument or do a data or smarts query')

  if options.similarityType=='AtomPairs':
    fpBuilder=BuildAtomPairFP
    fpDepickler=DepickleIntVectFP
    simMetric=DataStructs.DiceSimilarity
    dbName = os.path.join(options.dbDir,options.pairDbName)
    fpTableName = options.pairTableName
    fpColName = options.pairColName
  elif options.similarityType=='TopologicalTorsions':
    fpBuilder=BuildTorsionsFP
    fpDepickler=DepickleLongIntVectFP
    simMetric=DataStructs.DiceSimilarity
    dbName = os.path.join(options.dbDir,options.torsionsDbName)
    fpTableName = options.torsionsTableName
    fpColName = options.torsionsColName
  elif options.similarityType=='RDK':
    fpBuilder=BuildRDKitFP
    fpDepickler=DepickleRDKitFP
    simMetric=DataStructs.FingerprintSimilarity
    dbName = os.path.join(options.dbDir,options.fpDbName)
    fpTableName = options.fpTableName
    if not options.fpColName:
      options.fpColName='rdkfp'
    fpColName = options.fpColName
  elif options.similarityType=='2D':
    fpBuilder=Build2DFP
    fpDepickler=Depickle2DFP
    simMetric=DataStructs.FingerprintSimilarity
    dbName = os.path.join(options.dbDir,options.fpDbName)
    fpTableName = options.fpTableName
    if not options.fpColName:
      options.fpColName='autofragmentfp'
    fpColName = options.fpColName

  if options.smartsQuery:
    try:
      options.smartsQueryMol = Chem.MolFromSmarts(options.smartsQuery)
    except:
      logger.error('could not build query molecule from smarts "%s"'%options.smartsQuery)
      sys.exit(-1)

  if options.outF=='-':
    outF=sys.stdout
  else:
    outF = file(options.outF,'w+')
  
  molsOut=False
  if options.sdfOut:
    molsOut=True
    if options.sdfOut=='-':
      sdfOut=sys.stdout
    else:
      sdfOut = file(options.sdfOut,'w+')
  else:
    sdfOut=None
  if options.smilesOut:
    molsOut=True
    if options.smilesOut=='-':
      smilesOut=sys.stdout
    else:
      smilesOut = file(options.smilesOut,'w+')
  else:
    smilesOut=None
      
  if len(args):
    queryFilename=args[0]
    try:
      tmpF = file(queryFilename,'r')
    except IOError:
      logger.error('could not open query file %s'%queryFilename)
      sys.exit(1)

    if options.molFormat=='smiles':
      func=GetMolsFromSmilesFile
    elif options.molFormat=='sdf':
      func=GetMolsFromSDFile

    if not options.silent:
      msg='Reading query molecules'
      if fpBuilder: msg+=' and generating fingerprints'
      logger.info(msg)
    probes=[]
    i=0
    nms=[]
    for nm,smi,mol in func(queryFilename,None,options.nameProp):
      i+=1
      nms.append(nm)
      if not mol:
        logger.error('query molecule %d could not be built'%(i))
        probes.append((None,None))
        continue
      if fpBuilder:
        probes.append((mol,fpBuilder(mol)))
      else:
        probes.append((mol,None))
      if not options.silent and not i%1000:
        logger.info("  done %d"%i)
  else:
    queryFilename=None
    probes=None
    
  conn=None
  idName = options.molIdName
  ids=None
  if options.propQuery or options.smartsQuery:
    molDbName = os.path.join(options.dbDir,options.molDbName)
    conn = DbConnect(molDbName)
    curs = conn.GetCursor()
    if options.smartsQuery:
      if not options.silent: logger.info('Doing substructure query')
      if options.propQuery:
        where='where %s'%options.propQuery
      else:
        where=''
      curs.execute('select count(*) from molecules %(where)s'%locals())
      nToDo = curs.fetchone()[0]
      curs.execute('select %(idName)s,molpkl from molecules %(where)s'%locals())
      row=curs.fetchone()
      nDone=0
      ids=[]
      while row:
        id,molpkl = row
        m = Chem.Mol(str(molpkl))
        matched=m.HasSubstructMatch(options.smartsQueryMol)
        if options.negateQuery:
          matched = not matched
        if matched:
          ids.append(id)
        nDone+=1
        if not options.silent and not nDone%500:
          logger.info('  searched %d (of %d) molecules; %d hits so far'%(nDone,nToDo,len(ids)))
        row=curs.fetchone()
    elif options.propQuery:
      if not options.silent: logger.info('Doing property query')
      propQuery=options.propQuery.split(';')[0]
      curs.execute('select %(idName)s from molecules where %(propQuery)s'%locals())
      ids = [str(x[0]) for x in curs.fetchall()]
    if not options.silent: logger.info('Found %d molecules matching the query'%(len(ids)))

  t1=time.time()
  if probes:
    if not options.silent: logger.info('Finding Neighbors')
    conn = DbConnect(dbName)
    curs = conn.GetCursor()
    if ids:
      ids = [(x,) for x in ids]
      curs.execute('create temporary table _tmpTbl (%(idName)s varchar)'%locals())
      curs.executemany('insert into _tmpTbl values (?)',ids)
      join='join  _tmpTbl using (%(idName)s)'%locals()
    else:
      join=''

    curs.execute('select %(idName)s,%(fpColName)s from %(fpTableName)s %(join)s'%locals())
    topNLists = GetNeighborLists(probes,options.topN,curs,simMetric=simMetric,
                                 fpDepickler=fpDepickler)

    uniqIds=sets.Set()
    nbrLists = {}
    for i,nm in enumerate(nms):
      topNLists[i].reverse()
      scores=topNLists[i].GetPts()
      nbrNames = topNLists[i].GetExtras()
      nbrs = []
      for j,nbrNm in enumerate(nbrNames):
        if nbrNm is None:
          break
        else:
          uniqIds.add(nbrNm)
          nbrs.append((nbrNm,scores[j]))
      nbrLists[(i,nm)] = nbrs
    t2=time.time()
    if not options.silent: logger.info('The search took %.1f seconds'%(t2-t1))
    
    if not options.silent: logger.info('Creating output')
    ks = nbrLists.keys()
    ks.sort()
    if not options.transpose:
      for i,nm in ks:
        nbrs= nbrLists[(i,nm)]
        nbrTxt=options.outputDelim.join([nm]+['%s%s%.3f'%(x,options.outputDelim,y) for x,y in nbrs])
        print >>outF,nbrTxt
    else:
      labels = ['%s%sSimilarity'%(x[1],options.outputDelim) for x in ks]
      print >>outF,options.outputDelim.join(labels)
      for i in range(options.topN):
        outL = []
        for idx,nm in ks:
          nbr = nbrLists[(idx,nm)][i]
          outL.append(nbr[0])
          outL.append('%.3f'%nbr[1])
        print >>outF,options.outputDelim.join(outL)
    ids = list(uniqIds)
  else:
    if not options.silent: logger.info('Creating output')
    print >>outF,'\n'.join(ids)
  if molsOut and ids:
    molDbName = os.path.join(options.dbDir,options.molDbName)
    conn = DbConnect(molDbName)
    cns = conn.GetColumnNames('molecules')
    curs = conn.GetCursor()
    ids = [(x,) for x in ids]
    curs.execute('create temporary table _tmpTbl (%(idName)s varchar)'%locals())
    curs.executemany('insert into _tmpTbl values (?)',ids)
    curs.execute('select * from molecules join _tmpTbl using (%(idName)s)'%locals())
    row=curs.fetchone()
    while row:
      m = Chem.Mol(str(row[-1]))
      if sdfOut:
        m.SetProp('_Name',str(row[0]))
        print >>sdfOut,Chem.MolToMolBlock(m)
        for i in range(1,len(cns)-2):
          pn = cns[i]
          pv = str(row[i])
          print >>sdfOut,'> <%s>\n%s\n'%(pn,pv)
        print >>sdfOut,'$$$$'
      if smilesOut:
        print >>smilesOut,'%s %s'%(Chem.MolToSmiles(m,options.chiralSmiles),str(row[0]))
      row=curs.fetchone()
  if not options.silent: logger.info('Done!')
