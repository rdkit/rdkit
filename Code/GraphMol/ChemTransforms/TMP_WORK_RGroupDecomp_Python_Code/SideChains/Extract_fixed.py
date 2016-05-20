#  $Id: Extract.py 16238 2015-09-25 17:00:41Z kellebr5 $
#
#  Created by Greg Landrum, Dec. 2008
#
from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import rdqueries
from rdkit.Chem import AllChem
from rdkit.Avalon import pyAvalonTools
import os, re, sys
from collections import defaultdict

from rdkit.RDLogger import logger
logger = logger()

_version = "0.6.0"

# uncomment for lots of debugging
#def debug_print(*a, **kw):
#  print(*a, **kw)

# uncomment for no debugging
#def debug_print(*a, **kw):
#  pass
  
def GetMolDegenPoints(mol,options=None):
  """ 

  >>> GetMolDegenPoints(Chem.MolFromSmiles('Nc1ncccc1')) == set()
  True
  >>> GetMolDegenPoints(Chem.MolFromSmiles('Nc1ccccc1')) == set([2, 3, 5, 6])
  True

  cubane is a nice case:
  >>> GetMolDegenPoints(Chem.MolFromSmiles('NC12C3C4C5C3C1C5C24')) == set([2, 3, 5, 6, 7, 8]) 
  True
  """
  equivPoints = set()
  matches = mol.GetSubstructMatches(mol,uniquify=False)
  for match in matches:
    for idx,mIdx in enumerate(match):
      if idx!=mIdx:
        equivPoints.add(idx)
        equivPoints.add(mIdx)
  return equivPoints

def SymmetrizeSidechains(mols,core,sidechains,options=None):
  """
  indexing of arguments and results:
  [mol][group] -> (attachIdx,smiles)

  >>> smis = ['Nc1c(C)cccc1','Nc1ccccc1C']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> cores = [Chem.MolFromSmiles('Nc1ccccc1')]

  Preparing the results from GetSidechains:
  >>> d= GetSidechains(ms,cores)
  >>> chains = [x[0] for x in d]

  Starting point: non-symmetric chains:
  >>> chains
  [((2, '*C'),), ((6, '*C'),)]

  >>> schains = SymmetrizeSidechains(ms,cores[0],chains)
  >>> len(schains)==len(chains)
  True
  >>> schains
  [((2, '*C'),), ((2, '*C'),)]

  When the core is nonsymmetric, this doesn't do anything:
  >>> smis = ['Nc1c(C)nccc1','Nc1cnccc1C']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> cores = [Chem.MolFromSmiles('Nc1cnccc1')]
  >>> chains = GetSidechains(ms,cores)
  >>> schains = SymmetrizeSidechains(ms,cores[0],chains)
  >>> schains is chains
  True

  A more interesting variation, with two substituents per mol:
  >>> smis = ['Nc1c(C)c(O)ccc1','Nc1cccc(O)c1C']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> cores = [Chem.MolFromSmiles('Nc1ccccc1')]

  >>> d= GetSidechains(ms,cores)
  >>> chains = [x[0] for x in d]
  >>> schains = SymmetrizeSidechains(ms,cores[0],chains)
  >>> len(schains)==len(chains)
  True
  >>> chains
  [((2, '*C'), (3, '*O')), ((5, '*O'), (6, '*C'))]
  >>> schains
  [((2, '*C'), (3, '*O')), ((2, '*C'), (3, '*O'))]


  First "real" test case:
  >>> smis = ['Nc1c(C)cccc1','Nc1ccccc1C','Nc1c(Cl)cccc1C']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> cores = [Chem.MolFromSmiles('Nc1ccccc1')]

  >>> d= GetSidechains(ms,cores)
  >>> chains = [x[0] for x in d]
  >>> schains = SymmetrizeSidechains(ms,cores[0],chains)
  >>> len(schains)==len(chains)
  True
  >>> chains
  [((2, '*C'),), ((6, '*C'),), ((2, '*Cl'), (6, '*C'))]
  >>> schains
  [((2, '*C'),), ((2, '*C'),), ((2, '*C'), (6, '*Cl'))]

  >>> smis = ['Nc1c(C)cccc1','Nc1c(Cl)cccc1C','Nc1ccccc1C']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> cores = [Chem.MolFromSmiles('Nc1ccccc1')]
  >>> d= GetSidechains(ms,cores)
  >>> chains = [x[0] for x in d]
  >>> schains = SymmetrizeSidechains(ms,cores[0],chains)
  >>> len(schains)==len(chains)
  True
  >>> chains
  [((2, '*C'),), ((2, '*Cl'), (6, '*C')), ((6, '*C'),)]
  >>> schains
  [((2, '*C'),), ((2, '*C'), (6, '*Cl')), ((2, '*C'),)]

  This case demonstrates a wart in the current procedure. This is correct
  for the current implementation, but should ideally behave differently:
  >>> smis = ['Nc1c(C)cccc1','Nc1c(C)ccc(Cl)c1','Nc1cccc(Cl)c1']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> cores = [Chem.MolFromSmiles('Nc1ccccc1')]
  >>> d= GetSidechains(ms,cores)
  >>> chains = [x[0] for x in d]
  >>> schains = SymmetrizeSidechains(ms,cores[0],chains)
  >>> len(schains)==len(chains)
  True
  >>> chains
  [((2, '*C'),), ((2, '*C'), (5, '*Cl')), ((5, '*Cl'),)]
  >>> schains
  [((2, '*C'),), ((2, '*C'), (5, '*Cl')), ((3, '*Cl'),)]

  # really should be: [((2, '*C'),), ((2, '*C'), (5, '*Cl')), ((5, '*Cl'),)]

  Robust w.r.t to multiple attachment points in sidechain:
  >>> options,args = parser.parse_args([])
  >>> smis = ['C2C[C]c1ccccc1C2','C2C[N]c1ccccc1C2','NC(Cc1ccc(O)cc1)C(O)=O']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> cores = [Chem.MolFromSmiles('c1ccccc1C')]
  >>> d= GetSidechains(ms,cores,options)
  >>> d
  [(((None, 'core matches multiple times'),),), (((None, 'multiple attachment points in sidechain'),),), (((2, '*O'), (6, '*C(N)C(O)=O')),)]
  >>> chains = [x[0] for x in d]
  >>> schains = SymmetrizeSidechains(ms,cores[0],chains)
  >>> len(schains)==len(chains)
  True
  >>> schains
  [((None, 'core matches multiple times'),), ((None, 'multiple attachment points in sidechain'),), ((2, '*O'), (6, '*C(N)C(O)=O'))]

  """
  if len(mols)!=len(sidechains):
    raise ValueError('sidechains list must be as long as molecules list')
  if not core:
    return sidechains

  matches = core.GetSubstructMatches(core,
                                     uniquify=False,
                                     useChirality=True)
  debug_print("matches", matches)
  if len(matches)<=1:
    return sidechains

  degenPts = GetMolDegenPoints(core)

  # start by ordering the molecules by the number of attachment points on degenerate positions:
  newOrder=[]
  for i in range(len(mols)):
    nDegen=0
    for idx,chain in sidechains[i]:
      if idx in degenPts:
        nDegen+=1
    newOrder.append((nDegen,i))

  #core.SetProp("neworders", repr(newOrder))
  res=[]
  i=0
  # all the mols without chains in degenerate positions:
  while i<len(newOrder) and newOrder[i][0]==0:
    molIdx = newOrder[i][1]
    tmp = list(sidechains[molIdx])
    tmp.sort()
    res.append((molIdx,tuple(tmp)))
    i+=1

  # all the mols with single chains in degenerate positions:
  chainLocCount=defaultdict(lambda:defaultdict(int))
  posChainLists=defaultdict(list)
  while i<len(newOrder) and newOrder[i][0]==1:
    molIdx = newOrder[i][1]
    # assign the degenerate positions the lowest index:
    tChains=[]
    idx=None
    for tIdx,tChain in sidechains[molIdx]:
      if tIdx not in degenPts:
        debug_print (tIdx, "-", tIdx)
        tChains.append((tIdx,tChain))
      else:
        debug_print (tIdx, "-", idx)        
        idx,chain = tIdx,tChain
    assert idx is not None
    minId=idx
    for matchId,match in enumerate(matches):
      if match[idx]<minId:
        minId=match[idx]
    tChains.append((minId,chain))
    tChains.sort()
    res.append((molIdx,tuple(tChains)))
    chainLocCount[chain][minId]+=1
    posChainLists[minId].append(chain)
    i+=1

  # the remaining molecules have more than one chain in a degenerate position
  while i<len(newOrder) :
    count = newOrder[i][0]
    while i<len(newOrder) and newOrder[i][0]==count:
      molIdx = newOrder[i][1]

      localStates=[]
      seen = set()
      for matchId,match in enumerate(matches):
        chainState = []
        for idx,chain in sidechains[molIdx]:
          pIdx = match[idx]
          chainState.append((idx in degenPts,-chainLocCount[chain][pIdx],pIdx,chain))
        chainState.sort()
        localStates.append(chainState)
      localStates.sort()

      winner = localStates[0]
      tChains=[]
      for a,b,idx,chain in winner:
        tChains.append((idx,chain))
        chainLocCount[chain][idx]+=1
        posChainLists[idx].append(chain)
      tChains.sort()
      res.append((molIdx,tuple(tChains)))
      i+=1
    
  res.sort()
  res = [y for x,y in res]
  
  return res

def GetSidechains(mols,cores,options=None):
  """
  indexing of results:
  [mol][core][group] -> (attachIdx,smiles)

  >>> smis = ['Nc1c(C)cccc1','Nc1c(O)cccc1']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> cores = [Chem.MolFromSmiles('Nc1ccccc1')]
  >>> d= GetSidechains(ms,cores)
  >>> len(d)==len(ms)
  True
  >>> len(d[0])==len(cores)
  True
  >>> len(d[0][0])==len(d[1][0])
  True
  >>> idx0,smi0=d[0][0][0]
  >>> idx1,smi1=d[1][0][0]
  >>> idx0==idx1
  True
  >>> smi0
  '*C'
  >>> smi1
  '*O'
  
  This doesn't symmetrize things, so for symmetric scaffolds we aren't
  guaranteed to get the same R labels:
  >>> smis = ['Nc1c(C)cccc1','Nc1ccccc1C']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> cores = [Chem.MolFromSmiles('Nc1ccccc1')]
  >>> d= GetSidechains(ms,cores)
  >>> idx0,smi0=d[0][0][0]
  >>> idx1,smi1=d[1][0][0]
  >>> idx0==idx1
  False
  >>> smi0
  '*C'
  >>> smi1
  '*C'
  
  If the core doesn't match, empty lists are returned:
  >>> smis = ['Nc1c(C)cccc1']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> cores = [Chem.MolFromSmiles('Nc1cnccc1'),Chem.MolFromSmiles('Nc1cncnc1')]
  >>> d= GetSidechains(ms,cores)
  >>> len(d)==len(ms)
  True
  >>> len(d[0])==len(cores)
  True
  >>> d[0][0]
  ((None, 'no core matches'),)
  >>> d[0][1]
  ((None, 'no core matches'),)


  And we're robust w.r.t. bogus cores:
  >>> smis = ['Nc1c(C)cccc1']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> cores = [Chem.MolFromSmiles('Nc1ccccc1'),None]
  >>> d= GetSidechains(ms,cores)
  >>> len(d)==len(ms)
  True
  >>> len(d[0])==len(cores)
  True
  >>> len(d[0][0])
  1
  >>> len(d[0][1])
  2
  >>> d[0][0]
  ((2, '*C'),)
  >>> d[0][1]
  (None, 'bad core')

  and bogus molecules:
  >>> ms = [None,Chem.MolFromSmiles('Nc1c(C)cccc1')]
  >>> cores = [Chem.MolFromSmiles('Nc1cnccc1')]
  >>> d= GetSidechains(ms,cores)
  >>> len(d)==len(ms)
  True
  >>> len(d[0])==len(cores)
  True
  >>> len(d[0][0])
  2
  >>> d[0][0]
  (None, 'bad molecule')

  by default we accept doubly attached sidechains:
  >>> smis = ['C2Nc1c(C2)cccc1']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> cores = [Chem.MolFromSmiles('Cc1c(N)cccc1')]
  >>> d = GetSidechains(ms,cores)
  >>> d[0]
  (((0, '*C*'),),)

  but this can be toggled off with the options construct:
  >>> options,args = parser.parse_args([])
  >>> options.rejectDoubleAttachments=True
  >>> d = GetSidechains(ms,cores,options=options)
  >>> d[0]
  (((None, 'multiple attachment points in sidechain'),),)

  """
  dummyExpr=re.compile(r'\[([0-9]*)\*\]')
  res = []
  for molIdx,mol in enumerate(mols):
    if not mol:
      localRes =[(None,'bad molecule')]*len(cores)
      res.append(tuple(localRes))
      continue
    molMatched=False
    localRes =[]
    for i,core in enumerate(cores):
      if not core:
        localRes.append((None,'bad core'))
        continue
      coreRes=[]
      tmatches = mol.GetSubstructMatches(core,useChirality=True)
      if not tmatches:
        if not options or not options.silent:
          logger.warning("molecule %d did not match core %d"%(molIdx,i))
        coreRes=[(None,'no core matches')]
      elif len(tmatches)>1:
        if not options or not options.silent:
          logger.warning("core %d matches molecule multiple times"%i)
        coreRes=[(None,'core matches multiple times')]
      else:
        tMol = Chem.ReplaceCore(mol,core,False,True) # False => move dummies to sidechains
                                                     # labelByIndex => true
        if tMol:
          molMatched=True
          tSmi = Chem.MolToSmiles(tMol,True)
          tSmi = tSmi.split('.')
          for elem in tSmi:
            matches = dummyExpr.findall(elem)
            if not len(matches):
              if not options or not options.silent:
                logger.warning("sidechain with no attachment point label: %s"%elem)
            else: 
              if len(matches)>1:
                if options and options.rejectDoubleAttachments :
                  coreRes=((None,'multiple attachment points in sidechain'),)
                  break
                else:
                  if not options or not options.silent:
                    logger.warning("sidechain with multiple attachment point labels: %s"%elem)
              idx = matches[0]
              if not idx:
                idx = 0
              else:
                idx = int(idx)
              elem = dummyExpr.sub('*',elem)
              try:
                cSmi = pyAvalonTools.GetCanonSmiles(elem,True)
              except:
                pass
              else:
                elem = cSmi
              coreRes.append((idx,elem))
        
      localRes.append(tuple(coreRes))
    if not molMatched and options and not options.silent:
      logger.warning("molecule %d did not match any cores"%molIdx)

    res.append(tuple(localRes))
  return res

def ParseAvalonCore(filename,options=None):
  from xml.etree.cElementTree import parse
  tree = parse(filename)
  structs = tree.findall('QUERY/*/STRUCTURE')
  if options:
    options.cores=None
  if structs:
    cores=[]
    for struct in structs:
      if struct.attrib['ischecked']=='true':
        txt = str(struct.text)
        if txt[:3]=='\n\t\n':
          txt=txt[2:]
        mol = Chem.MolFromMolBlock(txt)
        if not mol:
          logger.error('Could not parse mol query:')
          print('----\n%s----'%txt, file=sys.stderr)
        else:
          cores.append(mol)
    if options:
      options.cores=cores
    return cores

def GetAvalonData(options,filename,
                  avalon = '/netfs/people/camm/landrgr1/highspeedAvalon/speed_export.sh'):
  raise NotImplementedError('needs to be updated')

def ProcessCoreLabels(cores,options):
  """

  >>> options,args = parser.parse_args([])
  >>> options.labelledCores=1
  >>> core = Chem.MolFromSmiles('c1ccccc1[2*]')
  >>> cores=[core]
  >>> ProcessCoreLabels(cores,options)
  >>> Chem.MolToSmiles(cores[0])
  'c1ccccc1'
  >>> cores[0].GetAtomWithIdx(5).HasProp('_RLabel')
  1
  >>> cores[0].GetAtomWithIdx(5).GetProp('_RLabel')
  '2'
  >>> core = Chem.MolFromSmarts('c1ccccc1[*:2]')
  >>> cores=[core]
  >>> ProcessCoreLabels(cores,options)
  >>> Chem.MolToSmiles(cores[0])
  'c1ccccc1'
  >>> cores[0].GetAtomWithIdx(5).HasProp('_RLabel')
  1
  >>> cores[0].GetAtomWithIdx(5).GetProp('_RLabel')
  '2'

  >>> core = Chem.MolFromSmarts('c1ccc[c,n]c1[*:2]')
  >>> cores=[core]
  >>> ProcessCoreLabels(cores,options)
  >>> Chem.MolToSmiles(cores[0])
  'c1ccccc1'
  
  """
  if not options.labelledCores:
    logger.warning('ProcessCoreLabels called without the labelledCores options set. Returning')
    return
  okSmarts=re.compile(r'\[[0-9]*\*[:,0-9]*\]')
  for coreIdx,core in enumerate(cores):
    labelDict={}
    dummies = [core.GetAtomWithIdx(x[0]) for x in core.GetSubstructMatches(Chem.MolFromSmiles('[*]'))]
    dummyNbrs=[-1]*(len(dummies))
    rmv = []
    for idx,dummy in enumerate(dummies):
      # start by making sure we should be removing this dummy (not all are removed):
      #if not okSmarts.match(dummy.GetSmarts()):
      #  print("skipping")
      #  continue
      
      # use the isotope we read in to set the R label.
      # this is appropriate because the mol file parser sets R group
      # labels the same way.
      # if we come in as SMARTS we use the atom mapping facility
      if not dummy.HasProp('molAtomMapNumber'):
        label=int(dummy.GetIsotope())
      else:
        label = int(dummy.GetProp('molAtomMapNumber'))
      nbrs = dummy.GetNeighbors()
      if not len(nbrs):
        logger.warning('attachment point has no neighbors')
        continue
      if len(nbrs)>1:
        logger.warning('attachment point has more than one neighbor')
        
      nbr = nbrs[0]
      if nbr.HasProp("_RLabel"):
        labels = nbr.GetProp("_RLabel") + "\t"
      else:
        labels = ""
      
      nbr.SetProp('_RLabel',labels + str(label))
      labelDict[nbr.GetIdx()] = labelDict.get(nbr.GetIdx(), [])  + [label]
      rmv.append(dummy.GetIdx())

    em = Chem.RWMol(core)
    rmv.sort(reverse=True)
    for idx in rmv:
      query = rdqueries.AtomNumGreaterQueryAtom(0)
      em.ReplaceAtom(idx, query)
    core = em#.GetMol()

    core._labelDict=labelDict
    cores[coreIdx]=core

def RunDecomposition(mols,options):
  """

  >>> options,args = parser.parse_args([])
  >>> options.labelledCores=1
  >>> options.requireLabels=1
  >>> smis = ['Nc1ccccc1','c1cc(F)ccc1','c1cc(F)c(C)cc1']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> options.cores = [Chem.MolFromSmiles('[1*]c1ccccc1')]
  >>> mols,d= RunDecomposition(ms,options)

  >>> len(d)==len(ms)
  True
  >>> len(d[0])==len(options.cores)
  True
  >>> len(d[0][0])
  1
  >>> d[0][0][0]
  (1, '*N')
  >>> d[1][0][0]
  (1, '*F')
  >>> d[2]
  (None, 'substitution at unmarked position')

  >>> options,args = parser.parse_args([])
  >>> options.labelledCores=1
  >>> options.requireLabels=1
  >>> smis = ['Nc1cc(F)ccc1','c1cc(F)cc(F)c1','c1c(C)cc(C)cc1','c1c(C)c(C)ccc1']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> options.cores = [Chem.MolFromSmiles('[1*]c1cc([2*])ccc1')]
  >>> mols,d= RunDecomposition(ms,options)

  >>> len(d)==len(ms)
  True
  >>> len(d[0])==len(options.cores)
  True
  >>> len(d[0][0])
  2
  >>> d[0][0]
  ((2, '*N'), (1, '*F'))
  >>> d[1][0]
  ((1, '*F'), (2, '*F'))
  >>> d[2][0]
  ((1, '*C'), (2, '*C'))
  >>> d[3]
  (None, 'substitution at unmarked position')

  >>> options,args = parser.parse_args([])
  >>> options.labelledCores=1
  >>> options.requireLabels=1
  >>> smis = ['CCN3C=Nc2c(Nc1cccc(Cl)c1)nc(nc23)N4CCC(O)CC4','COc4ccc(CNC2=NNc3ncnc(Nc1cccc(Cl)c1)c23)cc4Cl']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> options.cores = [Chem.MolFromSmiles('[1*]c1cccc([2*])c1')]
  >>> mols,d= RunDecomposition(ms,options)

  >>> len(d)==len(ms)
  True
  >>> len(d[0])==len(options.cores)
  True
  >>> d[0][0]
  ((2, '*Nc1nc(nc2c1N=CN2CC)N3CCC(O)CC3'), (1, '*Cl'))

  >>> d[1]
  (None, 'core matches multiple times')

  >>> options,args = parser.parse_args([])
  >>> options.labelledCores=1
  >>> options.requireLabels=1
  >>> smis = ['CC1CC1','C1CCC1']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> options.cores = [Chem.MolFromSmiles('C1CC1[1*]')]
  >>> mols,d= RunDecomposition(ms,options)

  >>> len(d)==len(ms)
  True
  >>> len(d[0])==len(options.cores)
  True
  >>> d[0][0]
  ((1, '*C'),)
  >>> d[1]
  (None, 'no core matches')

  >>> options,args = parser.parse_args([])
  >>> options.labelledCores=1
  >>> options.nonLabelledSubstituentHandling='FAIL'
  >>> smis = ['CCC1CC1','C1C(Cl)C1CC']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> options.cores = [Chem.MolFromSmiles('C1CC1C[1*]')]
  >>> mols,d= RunDecomposition(ms,options)
  >>> len(d)==len(ms)
  True
  >>> len(d[0])==len(options.cores)
  True
  >>> d[0][0]
  ((1, '*C'),)
  >>> d[1]
  (None, 'substitution at unmarked position')

  >>> options,args = parser.parse_args([])
  >>> options.labelledCores=1
  >>> options.nonLabelledSubstituentHandling='PASS'
  >>> smis = ['CCC1CC1','C1C(Cl)C1CC']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> options.cores = [Chem.MolFromSmiles('C1CC1C[1*]')]
  >>> mols,d= RunDecomposition(ms,options)
  >>> len(d)==len(ms)
  True
  >>> len(d[0])==len(options.cores)
  True
  >>> d[0][0]
  ((1, '*C'),)
  >>> d[1][0]
  ((101, '*Cl'), (1, '*C'))


  >>> options,args = parser.parse_args([])
  >>> options.labelledCores=1
  >>> options.nonLabelledSubstituentHandling='IGNORE'
  >>> smis = ['CCC1CC1','C1C(Cl)C1CC']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> options.cores = [Chem.MolFromSmiles('C1CC1C[1*]')]
  >>> mols,d= RunDecomposition(ms,options)
  >>> len(d)==len(ms)
  True
  >>> len(d[0])==len(options.cores)
  True
  >>> d[0][0]
  ((1, '*C'),)
  >>> d[1][0]
  ((1, '*C'),)




  """
  if not options.cores:
    logger.error('no core information found')
    sys.exit(1)

  coreCoreMatches=[]
  if options.labelledCores:
    ProcessCoreLabels(options.cores,options)

    for core in options.cores:
      matches = core.GetSubstructMatches(core,uniquify=False)
      #assert matches, "%s %s"%(Chem.MolToSmarts(core), Chem.MolToSmiles(core))
      debug_print ("matches", matches)
      if not matches:
        # pathology happens at times with query features in the core
        matches=(tuple(range(core.GetNumAtoms())),)
      coreCoreMatches.append(matches)
    debug_print( 'CCM:',coreCoreMatches )

    # commented out as part of provisional fix for CIXFW-98
    # if options.symmetrize:
    #   for coreIdx,core in enumerate(options.cores):
    #     equivs = defaultdict(list)
    #     for match in coreCoreMatches[coreIdx]:
    #       for i,j in enumerate(match):
    #         if j>i:
    #           if core.GetAtomWithIdx(i).HasProp('_RLabel') ^ core.GetAtomWithIdx(j).HasProp('_RLabel'):
    #             raise ValueError('when using symmetrization with labelled cores, either all symmetric attachment points must be labelled or all must be unlabelled')

  res = GetSidechains(mols,options.cores,options)
  if options.symmetrize:
    if not options.silent:
      logger.info('Symmetrizing R groups')
    for coreIdx,core in enumerate(options.cores):
      chains = [x[coreIdx] for x in res]
      sChains=SymmetrizeSidechains(mols,core,chains,options=options)
      for i,row in enumerate(res):
        row = list(row)
        row[coreIdx]=sChains[i]
        res[i] = row

  if options.labelledCores:
    if not options.silent:
      logger.info('Renumbering to match input ordering')

    for coreIdx,core in enumerate(options.cores):

      # for symmetric cores where we're labelling, we have to deal with the real PITA
      # of a core like *C1CC1. Since the dummy gets subsumed into its attached atom
      # when we do the replacement, this could match C1CC1F three different ways. The one
      # we get isn't determined by logic, but by atom ordering. We want to present
      # a sensible result, so we have to do some work:
      ccM = coreCoreMatches[coreIdx]
      chains = [x[coreIdx] for x in res]
      for i,row in enumerate(res):
        core_atoms_used_count = {}

        possibleMappings=[0]*len(ccM)
        debug_print ('#------------------:',list(chains[i]))
        chain = chains[i]
        chain = list(chain)
        for cidx,(idx,smi) in enumerate(chain):
          if idx is None:
            res[i] = (idx,smi)
            chain=None
            break
          for midx,mapping in enumerate(ccM):
            if core.GetAtomWithIdx(mapping[idx]).HasProp('_RLabel'):
              #print '    ',cidx,midx,idx,'->',mapping[idx]
              possibleMappings[midx] += 1
        if chain is not None:
          debug_print ('POSS:',possibleMappings)
          best = sorted(zip(possibleMappings,range(len(possibleMappings))))
          #raise ValueError(best)
          best = best[-1][1]
          bestMapping = ccM[best]
          #print 'BEST:',bestMapping

          for cidx,(idx,smi) in enumerate(chain):
            if not core.GetAtomWithIdx(bestMapping[idx]).HasProp('_RLabel'):
              if options.requireLabels or options.nonLabelledSubstituentHandling=='FAIL':
                if not options or not options.silent:
                  logger.warning('skipping molecule %d because it has a substitution at an unmarked position'%i)
                res[i] = (None,'substitution at unmarked position')
                chain=None
                break
              elif options.nonLabelledSubstituentHandling=='IGNORE':
                idx=-1
              else:
                idx+=100
            else:
              count = core_atoms_used_count.get(idx, 0)
              labels = core.GetAtomWithIdx(bestMapping[idx]).GetProp('_RLabel').split("\t")
              if count >= len(labels):
                label = labels[0]
              else:
                label = labels[count]
              core_atoms_used_count[idx] = count+1
                
              idx=int(label)#core.GetAtomWithIdx(bestMapping[idx]).GetProp('_RLabel'))
            if idx>=0:
              chain[cidx]=(idx,smi)
            else:
              chain[cidx]=None
        row = list(row)
        if chain is not None:
          chain = tuple([x for x in chain if x is not None])
          row[coreIdx]=chain
          res[i] = row
  return mols,res

def CreateOutput(mols,res,options):
  """

  repeat a test case we used above:

  >>> from rdkit.six.moves import cStringIO as StringIO
  >>> options,args = parser.parse_args([])
  >>> options.labelledCores=1
  >>> options.requireLabels=1
  >>> smis = ['CCN3C=Nc2c(Nc1cccc(Cl)c1)nc(nc23)N4CCC(O)CC4','COc4ccc(CNC2=NNc3ncnc(Nc1cccc(Cl)c1)c23)cc4Cl']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> options.cores = [Chem.MolFromSmiles('[1*]c1cccc([2*])c1')]
  >>> mols,d= RunDecomposition(ms,options)
  >>> sio = StringIO()
  >>> options.outFile=sio
  >>> options.outputFormat='csv'
  >>> CreateOutput(mols,d,options)
  >>> print(sio.getvalue())
  compound_id,smiles,Core,R1,R2
  Mol_1,CCn1cnc2c1nc(N1CCC(O)CC1)nc2Nc1cccc(Cl)c1,0,*Cl,*Nc1nc(nc2c1N=CN2CC)N3CCC(O)CC3
  Mol_2,COc1ccc(CNc2n[nH]c3ncnc(Nc4cccc(Cl)c4)c23)cc1Cl,ERROR: core matches multiple times
  <BLANKLINE>

  without labels:

  >>> ms.append(Chem.MolFromSmiles('[C]1CCCc2c1cccc2'))
  >>> sio = StringIO()
  >>> options.outFile=sio
  >>> options.labelledCores=0
  >>> options.requireLabels=0
  >>> options.cores = [Chem.MolFromSmiles('c1ccccc1')]
  >>> mols,d= RunDecomposition(ms,options)
  >>> CreateOutput(mols,d,options)
  >>> print(sio.getvalue())
  compound_id,smiles,Core,R1,R5
  Mol_1,CCn1cnc2c1nc(N1CCC(O)CC1)nc2Nc1cccc(Cl)c1,0,*Nc1nc(nc2c1N=CN2CC)N3CCC(O)CC3,*Cl
  Mol_2,COc1ccc(CNc2n[nH]c3ncnc(Nc4cccc(Cl)c4)c23)cc1Cl,ERROR: core matches multiple times
  Mol_3,[C]1CCCc2ccccc21,ERROR: multiple attachment points in sidechain
  <BLANKLINE>
  
  robust r.s.t bogus core when includeCoreSmiles is ON

  >>> sio = StringIO()
  >>> options.outFile=sio
  >>> options.includeCoreSmiles = '1'
  >>> ms.append(Chem.MolFromSmiles('CCC'))
  >>> mols,d= RunDecomposition(ms,options)
  >>> CreateOutput(mols,d,options)
  >>> print(sio.getvalue())
  compound_id,smiles,Core,R1,R5
  Mol_1,CCn1cnc2c1nc(N1CCC(O)CC1)nc2Nc1cccc(Cl)c1,c1ccccc1,*Nc1nc(nc2c1N=CN2CC)N3CCC(O)CC3,*Cl
  Mol_2,COc1ccc(CNc2n[nH]c3ncnc(Nc4cccc(Cl)c4)c23)cc1Cl,ERROR: core matches multiple times
  Mol_3,[C]1CCCc2ccccc21,ERROR: multiple attachment points in sidechain
  Mol_4,CCC,ERROR: no core matches
  <BLANKLINE>


  """
  pns = list(mols[0].GetPropNames())
  oddE = re.compile(r'[^a-zA-Z0-9/]')
  outNs = [oddE.sub('_',x) for x in pns]
  if options.outputFormat=='spotfire':
    heads = ['sid','compound_id','smiles','Core','RIdx','RSmiles']+outNs
    print('|'.join(heads), file=options.outFile)
  else:
    rIndices = set()
    for molRes in res:
      if molRes[0] is None:
        continue
      for row in molRes:
        if not row or row[0] is None: continue
        for rIdx,rSmi in row:
          if rIdx is None: continue   # when core matches multiple times
          rIndices.add(rIdx)
    rIndices = list(rIndices)
    rIndices.sort()
    rMap={}
    for i,idx in enumerate(rIndices): rMap[idx]=i

    heads = ['compound_id','smiles','Core']
    for idx in rIndices:
      if not options.labelledCores:
        idx+=1
      heads.extend(['R%d'%(idx)])
    heads.extend(pns)
    print(','.join(heads), file=options.outFile)

  for i,molRes in enumerate(res):
    molPs= []
    for x in pns:
      try:
        v = mols[i].GetProp(x)
      except KeyError:
        v ='N/A'
      molPs.append(v)
    if options.outputFormat=='spotfire':
      # FIX: not complete anymore
      for coreIdx,row in enumerate(molRes):
        if row is None: continue
        row = list(row)
        row.sort()
        for rIdx,rSmi in row:
          if not options.labelledCores:
            rIdx += 1
          smi = Chem.MolToSmiles(mols[i])
          tRow = [mols[i].GetProp('_Name')+'_'+str(coreIdx)+"_"+str(rIdx),mols[i].GetProp('_Name'),smi]
          if not options.includeCoreSmiles:
            tRow.append(str(coreIdx))
          else:
            tRow.append(Chem.MolToSmiles(options.cores[coreIdx],True))
          tRow += [str(rIdx),rSmi]
          tRow += molPs
          print('|'.join(tRow), file=options.outFile)
    else:
      smi = Chem.MolToSmiles(mols[i])
      if mols[i].HasProp('_Name'):
        nm = mols[i].GetProp('_Name')
      else:
        nm = 'Mol_%d'%(i+1)
      tRow = [nm,smi]
      if molRes[0] is None:
        # this is an error condition
        print(','.join(tRow+['ERROR: '+molRes[1]]), file=options.outFile)
      else:
        for coreIdx,row in enumerate(molRes):
          if not row and not options.keepNonmatchingMols:
            continue
          if row:
            rs = ['*[H]']*len(rIndices)
            if not options.includeCoreSmiles:
              outD = [str(coreIdx)]
            else:
              tm = AllChem.ReplaceSidechains(mols[i],options.cores[coreIdx])
              if tm:
                tm = Chem.DeleteSubstructs(tm,Chem.MolFromSmarts('[#0]'))
                outD = [Chem.MolToSmiles(tm,True)]
              else:
                outD = ['']

            if row[0] is None:
              # bad result, due to for instance multiple attachment points in
              # sidechain
              outD=['']
              rs=['']*len(rIndices)
              print(','.join(tRow+['ERROR: '+row[1]]), file=options.outFile)
              break
            if len(row) == 1 and row[0][0] is None:
              # bad result, due to core matches multiple times
              outD=['']
              rs=['']*len(rIndices)
              print(','.join(tRow+['ERROR: '+row[0][1]]), file=options.outFile)
              break
            for rIdx,rSmi in row:
              # this bit of hackery is because the avalon depictor doesn't like SMILES with asterixes in square brackets:
              rSmi=rSmi.replace('[*]','*')
              rs[rMap[rIdx]]=rSmi
          else:
            outD=['']
            rs=['']*len(rIndices)
          outD.extend(rs)
          print(','.join(tRow+outD+molPs), file=options.outFile)


#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest
  return doctest.testmod(sys.modules["__main__"])

# ---- ---- ---- ----  ---- ---- ---- ----  ---- ---- ---- ----  ---- ---- ---- ---- 
from optparse import OptionParser
parser=OptionParser("Extract.py",version='%prog '+_version)
parser.add_option('--runTests','--test',default=False,action='store_true',
                  help='run internal tests')
parser.add_option('--core',default='',
                  help='specify the core')
parser.add_option('--coreFormat',default='smarts',
                  choices=['smarts','smiles','sdf','avalon','mol'],
                  help='specify the format the core is in')
parser.add_option('--labelledCores',default=False,action='store_true',
                  help='allow cores to have attachment points labelled')
parser.add_option('--requireLabels',default=False,action='store_true',
                  help='only allow substitutions at labelled points')
parser.add_option('--dataFormat',default='other',
                  choices=['other'],
                  help='specify the format the dataset is in')
parser.add_option('--silent',default=False,action='store_true',
                  help='be quiet')
parser.add_option('--outputFormat','--outFormat',default='spotfire',
                  choices=['spotfire','csv'],
                  help='specify the output format')
parser.add_option('--outFile','--outF','-o',default='-',
                  help='outfile name')
parser.add_option('--saveAvalon',default='',
                  help='specify the filename to save the avalon results')
parser.add_option('--symmetrize','--symm','-s',default=False,action='store_true',
                  help='attempt to symmetrize the R group assignments (experimental)')
parser.add_option('--includeCoreSmiles',default=False,action='store_true',
                  help='include the SMILES of the core in the output instead of its index')
parser.add_option('--keepNonmatchingMols',default=False,action='store_true',
                  help='retain molecules that do not match any of the known cores, or that have disallowed substituents, in the output.')
parser.add_option('--rejectDoubleAttachments',default=True,action='store_false',
                  help='do not reject molecules where sidechains end up with two attachment points')
parser.add_option('--nonLabelledSubstituentHandling',default='PASS',choices=['PASS','IGNORE','FAIL'],
                  help='a more flexible way of specifying what to do with substituents in non-labelled positions when labelledCores is true. Options are: PASS the substituents through anyway, IGNORE the substituents, or FAIL the molecule (this is equivalent to having the requireLabels option set).')

if __name__=='__main__':
  options,args = parser.parse_args()

  if options.runTests:
    failed,tried = _test()
    sys.exit(failed)
    
  if len(args)!=1:
    parser.error('please a data filename')

  if not os.path.exists(args[0]):
    logger.error('could not open input file: %s'%args[0])
    sys.exit(1)

  if options.outFile=='-' or not options.outFile:
    options.outFile = sys.stdout
  else:
    try:
      options.outFile = open(options.outFile,'w+')
    except IOError:
      logger.error('could not open output file: %s'%options.outFile)
      sys.exit(-1)

  if not options.core and options.dataFormat=='avalon':
    ParseAvalonCore(args[0],options)
  elif options.core:
    if options.coreFormat=='smarts':
      options.cores=[Chem.MolFromSmarts(options.core)]
    elif options.coreFormat=='smiles':
      options.cores=[Chem.MolFromSmiles(options.core)]
    elif options.coreFormat=='mol':
      if not os.path.exists(options.core):
        logger.error('could not open core file: %s'%options.core)
        sys.exit(1)
      options.cores=[Chem.MolFromMolFile(options.core)]
    elif options.coreFormat=='sdf':
      if not os.path.exists(options.core):
        logger.error('could not open core file: %s'%options.core)
        sys.exit(1)
      suppl = Chem.SDMolSupplier(options.core,sanitize=False,removeHs=False)
      options.cores=[x for x in suppl]
      for i,core in enumerate(options.cores):
        Chem.SanitizeMol(core)
        options.cores[i] = Chem.MergeQueryHs(core)
    elif options.coreFormat=='avalon':
      if not os.path.exists(options.core):
        logger.error('could not open core file: %s'%options.core)
        sys.exit(1)
      ParseAvalonCore(options.core,options)
  if not options.cores:
    logger.error('no core information found')
    sys.exit(1)

  if options.dataFormat=='avalon':
    raise NotImplementedError('direct avalon import not currently supported')
    mols = GetAvalonData(options,args[0])
  else:
    from rdkit.Novartis.Plumbing import Plumbing
    from rdkit.Novartis.Plumbing.Types.Mol import getRDKMol
    mols = Plumbing.LoadMols(args[0])
    mols = [getRDKMol(x) for x in mols]

  mols,res=RunDecomposition(mols,options)
  CreateOutput(mols,res,options)
  

        
          
      
