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
_version = "0.1.3"
_usage="""
 SearchDb [optional arguments] <sdfilename>

     The sd filename argument can be either an SD file or an MDL mol 
     file.
     
 """
import RDConfig
from Dbase.DbConnection import DbConnect

from RDLogger import logger
logger=logger()
import cPickle


# ----------------------------------------
# ATOM PAIRS
from DataStructs import SparseIntVect
def DepickleIntVectFP(pkl):
  fp = cPickle.loads(str(pkl))
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
  from Chem.Fingerprints.FingerprintMols import FingerprintMol
  fp=FingerprintMol(mol)
  return fp

def GetNeighborLists(probeFps,topN,cursor,
                     simMetric=SparseIntVect.DiceSimilarity,
                     fpDepickler=DepickleIntVectFP,
                     silent=False):

  from DataStructs.TopNContainer import TopNContainer
  nbrLists = [TopNContainer(topN) for x in range(len(probeFps))]

  row = curs.fetchone()
  nDone=0
  while row:
    nDone+=1
    if not silent and not nDone%1000: logger.info('  searched %d rows'%nDone)
    nm,pkl = row
    fp=fpDepickler(pkl)
    for i in range(len(probeFps)):
      pfp = probeFps[i]
      if pfp is not None:
        score = simMetric(probeFps[i],fp)
        nbrLists[i].Insert(score,nm)
    row = curs.fetchone()
  return nbrLists


# ---- ---- ---- ----  ---- ---- ---- ----  ---- ---- ---- ----  ---- ---- ---- ---- 
import os
from optparse import OptionParser
parser=OptionParser(_usage,version='%prog '+_version)
parser.add_option('--dbDir',default='/db/CADD/NOV_Q2_2007/rdk_db',
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
parser.add_option('--fpTableName',default='fingerprints',
                  help='name of the 2D fingerprints table')
parser.add_option('--fpColName',default='autofragmentfp',
                  help='name of the 2D fingerprint column')
parser.add_option('--descrDbName',default='Descriptors.sqlt',
                  help='name of the descriptor database')
parser.add_option('--descrTableName',default='descriptors_v1',
                  help='name of the descriptor table')
parser.add_option('--descriptorCalcFilename',default=os.path.join(RDConfig.RDBaseDir,'Projects',
                                                                  'QuickMolDB','moe_like.dsc'),
                  help='name of the file containing the descriptor calculator')
parser.add_option('--similarityType',default='2D',choices=['2D','AtomPairs','TopologicalTorsions'],
                  help='Choose the type of similarity to use, possible values: 2D, AtomPairs, and TopologicalTorsions. The default is %default')
parser.add_option('--outputDelim',default=',',
                  help='the delimiter for the output file. The default is %default')
parser.add_option('--topN',default=20,type='int',
                  help='the number of neighbors to keep for each query compound. The default is %default')

parser.add_option('--outF',default='-',
                  help='The name of the output file. The default is the console (stdout).')

parser.add_option('--transpose',default=False,action="store_true",
                  help='print the results out in a transposed form: e.g. neighbors in rows and probe compounds in columns')


parser.add_option('--silent',default=False,action='store_true',
                  help='Do not generate status messages.')

if __name__=='__main__':
  import sys,getopt,time
  import Chem
  
  options,args = parser.parse_args()
  if len(args)!=1:
    parser.error('please provide a query filename argument')

  if options.similarityType=='AtomPairs':
    fpBuilder=BuildAtomPairFP
    fpDepickler=DepickleIntVectFP
    #simMetric=SparseIntVect.DiceSimilarity
    simMetric=DataStructs.DiceSimilarity
    dbName = os.path.join(options.dbDir,options.pairDbName)
    fpTableName = options.pairTableName
    fpColName = options.pairColName
  elif options.similarityType=='TopologicalTorsions':
    fpBuilder=BuildTorsionsFP
    fpDepickler=DepickleIntVectFP
    #simMetric=SparseIntVect.DiceSimilarity
    simMetric=DataStructs.DiceSimilarity
    dbName = os.path.join(options.dbDir,options.torsionsDbName)
    fpTableName = options.torsionsTableName
    fpColName = options.torsionsColName
  elif options.similarityType=='2D':
    fpBuilder=BuildRDKitFP
    fpDepickler=DepickleRDKitFP
    simMetric=DataStructs.FingerprintSimilarity
    dbName = os.path.join(options.dbDir,options.fpDbName)
    fpTableName = options.fpTableName
    fpColName = options.fpColName
    
  if options.outF=='-':
    outF=sys.stdout
  else:
    outF = file(options.outF,'w+')
      
  queryFilename=args[0]

  try:
    tmpF = file(queryFilename,'r')
  except IOError:
    logger.error('could not open query file %s'%queryFilename)
    sys.exit(1)
      
  if not options.silent: logger.info('Reading query molecules')
  suppl = Chem.SDMolSupplier(queryFilename)
  queryMols = [x for x in suppl]
  suppl=None

  if not options.silent: logger.info('Generating fingerprints')
  probeFps=[]
  for i in range(len(queryMols)):
    mol = queryMols[i]
    if not mol:
      logger.error('query molecule %d could not be built'%(i+1))
      probeFps.append(None)
      continue
    probeFps.append(fpBuilder(mol))

    
  if not options.silent: logger.info('Finding Neighbors')
  t1=time.time()
  conn = DbConnect(dbName)
  curs = conn.GetCursor()
  idName = options.molIdName
  curs.execute('select %(idName)s,%(fpColName)s from %(fpTableName)s'%locals())
  topNLists = GetNeighborLists(probeFps,options.topN,curs,simMetric=simMetric,
                               fpDepickler=fpDepickler)

  nbrLists = {}
  for i in range(len(queryMols)):
    mol = queryMols[i]
    if not mol: continue
    if mol.HasProp('_Name'):
      nm =mol.GetProp('_Name').strip()
    else:
      nm = 'Query_%d'%(i+1)
    scores=topNLists[i].GetPts()
    nbrNames = topNLists[i].GetExtras()
    nbrs = zip(nbrNames,scores)
    nbrs.reverse()
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
      
  if not options.silent: logger.info('Done!')
