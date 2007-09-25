#  $Id$
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
""" Implementation of the RECAP algorithm from Lewell et al. JCICS *38* 511-522 (1998)

The published algorithm is implemented more or less without
modification. The results are returned as a hierarchy of nodes instead
of just as a set of fragments. The hope is that this will allow a bit
more flexibility in working with the results.

For example:
>>> m = Chem.MolFromSmiles('C1CC1Oc1ccccc1-c1ncc(OC)cc1')
>>> res = Recap.RecapDecompose(m)
>>> res
<Chem.Recap.RecapHierarchyNode object at 0x00CDB5D0>
>>> res.children.keys()
['[Du]C1CC1', '[Du]c1ccccc1-c1ncc(OC)cc1']
>>> res.GetAllChildren().keys()
['[Du]c1ccccc1[Du]', '[Du]c1ccc(OC)cn1', '[Du]C1CC1', '[Du]c1ccccc1-c1ncc(OC)cc1']

To get the standard set of RECAP results, use GetLeaves():
>>> leaves=res.GetLeaves()
>>> leaves.keys()
['[Du]c1ccccc1[Du]', '[Du]c1ccc(OC)cn1', '[Du]C1CC1']
>>> leaf = leaves['[Du]C1CC1']
>>> leaf.mol
<Chem.rdchem.Mol object at 0x00CBE0F0>


"""
import Chem
from Chem import rdChemReactions as Reactions

# These are the definitions that will be applied to fragment molecules:
reactionDefs = (
  "[C:1](=!@[O:2])!@[N:3]>>[X][C:1]=[O:2].[X][N:3]", # amide

  "[C:1](=!@[O:2])!@[O:3]>>[X][C:1]=[O:2].[O:3][X]", # ester

  "[N:1;D2,D3]!@C(!@=O)!@[N:2,D2,D3]>>[X][N:1].[N:2][X]", # urea

  "[N;!D1](-[*:1])-!@[*:2]>>[X][*:1].[*:2][X]", # amines
  #"[N;!D1](!@[*:1])!@[*:2]>>[X][*:1].[*:2][X]", # amines

  "[#6:1]-!@O-!@[#6:2]>>[#6:1][X].[X][#6:2]", # ether

  "[C:1]=!@[C:2]>>[C:1][X].[X][C:2]", # olefin

  "[n:1]-!@[C:2]>>[n:1][X].[C:2][X]", # aromatic nitrogen - aliphatic carbon

  "O=C-@[N:1]-!@[C:2]>>[N:1][X].[C:2][X]", # lactam nitrogen - aliphatic carbon

  "[c:1]-!@[c:2]>>[c:1][X].[X][c:2]", # aromatic carbon - aromatic carbon

  "[n:1]-!@[c:2]>>[n:1][X].[X][c:2]", # aromatic nitrogen - aromatic carbon *NOTE* this is not part of the standard recap set.

  "[N:1;D2,D3]-!@[S:2](=[O:3])=[O:4]>>[N:1][X].[X][S:2](=[O:3])=[O:4]", # sulphonamide
  )

reactions = tuple([Reactions.ReactionFromSmarts(x) for x in reactionDefs])


class RecapHierarchyNode(object):
  """ This class is used to hold the Recap hiearchy
  """
  mol=None
  children=None
  parents=None
  smiles = None
  def __init__(self,mol):
    self.mol=mol
    self.children = {}
    self.parents = {}

  def GetAllChildren(self):
    " returns a dictionary, keyed by SMILES, of children "
    res = {}
    for smi,child in self.children.iteritems():
      res[smi] = child
      child._gacRecurse(res,terminalOnly=False)
    return res

  def GetLeaves(self):
    " returns a dictionary, keyed by SMILES, of leaf (terminal) nodes "
    res = {}
    for smi,child in self.children.iteritems():
      if not len(child.children):
        res[smi] = child
      else:
        child._gacRecurse(res,terminalOnly=True)
    return res

  def getUltimateParents(self):
    """ returns all the nodes in the hierarchy tree that contain this
        node as a child
    """
    if not self.parents:
      res = [self]
    else:
      res = []
      for p in self.parents.values():
        for uP in p.getUltimateParents():
          if uP not in res:
            res.append(uP)
    return res
  
  def _gacRecurse(self,res,terminalOnly=False):
    for smi,child in self.children.iteritems():
      if not terminalOnly or not len(child.children):
        res[smi] = child
      child._gacRecurse(res,terminalOnly=terminalOnly)

  def _remove(self):
    " removes this entry from the hierarchy (breaks reference loops) "
    for parent in self.parents.values():
      #parent = parent()
      if parent:
	ks = parent.children.keys()[:]
	for k in ks:
	  if parent.children[k] is self:
	    del parent.children[k]

  def __del__(self):
    self._remove()
    self.children={}
    self.parent={}
    self.mol=None


def RecapDecompose(mol,allNodes=None):
  """ returns the recap decomposition for a molecule """
  mSmi = Chem.MolToSmiles(mol,1)

  if allNodes is None:
    allNodes={}
  if allNodes.has_key(mSmi):
    return allNodes[mSmi]

  res = RecapHierarchyNode(mol)
  res.smiles =mSmi
  activePool={mSmi:res}
  allNodes[mSmi]=res
  for rxnIdx,reaction in enumerate(reactions):
    localRes = {}
    #print '>',rxnIdx,len(activePool.keys())
    missed = {}
    while activePool:
      nSmi = activePool.keys()[0]
      node = activePool.pop(nSmi)
      ps = reaction.RunReactants((node.mol,))
      if not ps:
        #print '  !',nSmi
	localRes[nSmi]=node
      else:
        rxnApplied=False
        #print '  .',nSmi
	for prodSeq in ps:
	  seqOk=True
	  # we want to disqualify small fragments, so sort the product sequence by size
	  # and then look for "forbidden" fragments
	  prodSeq = [(prod.GetNumAtoms(onlyHeavy=True),prod) for prod in prodSeq]
	  prodSeq.sort()
	  for nats,prod in prodSeq:
	    pSmi = Chem.MolToSmiles(prod,1)
            # don't forget after replacing dummy atoms to remove any empty
            # branches:
	    if pSmi.replace('[Du]','').replace('()','') in ('','C','CC','CCC'):
	      seqOk=False
	      break
	    prod.pSmi = pSmi
	  if seqOk:
            rxnApplied=True
  	    for nats,prod in prodSeq:
	      pSmi = prod.pSmi
	      pNode = allNodes.get(pSmi,RecapHierarchyNode(prod))
              pNode.smiles=pSmi
              pNode.parents[nSmi]=node
	      node.children[pSmi]=pNode
	      activePool[pSmi] = pNode
	      allNodes[pSmi]=pNode
        if not rxnApplied:
          localRes[nSmi]=node
    activePool=localRes
  return res

      
# ------- ------- ------- ------- ------- ------- ------- -------
# Begin testing code
if __name__=='__main__':
  import unittest
  class TestCase(unittest.TestCase):
    def test1(self):
      m = Chem.MolFromSmiles('C1CC1Oc1ccccc1-c1ncc(OC)cc1')
      res = RecapDecompose(m)
      self.failUnless(res)
      self.failUnless(len(res.children.keys())==2)
      self.failUnless(len(res.GetAllChildren().keys())==4)
      self.failUnless(len(res.GetLeaves().keys())==3)
    def test2(self):
      m = Chem.MolFromSmiles('CCCOCCC')
      res = RecapDecompose(m)
      self.failUnless(res)
      self.failUnless(res.children=={})
    def test3(self):
      allNodes={}
      m = Chem.MolFromSmiles('c1ccccc1-c1ncccc1')
      res = RecapDecompose(m,allNodes=allNodes)
      self.failUnless(res)
      self.failUnless(len(res.children.keys())==2)
      self.failUnless(len(allNodes.keys())==3)

      m = Chem.MolFromSmiles('COc1ccccc1-c1ncccc1')
      res = RecapDecompose(m,allNodes=allNodes)
      self.failUnless(res)
      self.failUnless(len(res.children.keys())==2)
      # we get two more nodes from that:
      self.failUnless(len(allNodes.keys())==5)
      self.failUnless(allNodes.has_key('[Du]c1ccccc1OC'))
      self.failUnless(allNodes.has_key('[Du]c1ccccc1'))
      
      m = Chem.MolFromSmiles('C1CC1Oc1ccccc1-c1ncccc1')
      res = RecapDecompose(m,allNodes=allNodes)
      self.failUnless(res)
      self.failUnless(len(res.children.keys())==2)
      self.failUnless(len(allNodes.keys())==9)
      
    def testSFNetIssue1801871(self):
      m = Chem.MolFromSmiles('c1ccccc1OC(Oc1ccccc1)Oc1ccccc1')
      res = RecapDecompose(m)
      self.failUnless(res)
      self.failUnless(len(res.GetLeaves())==2)
      ks = res.GetLeaves().keys()
      self.failIf('[Du]C([Du])[Du]' in ks)
      self.failUnless('[Du]c1ccccc1' in ks)
      self.failUnless('[Du]C([Du])Oc1ccccc1' in ks)
      


  unittest.main()

