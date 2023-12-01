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
>>> from rdkit import Chem
>>> from rdkit.Chem import Recap
>>> m = Chem.MolFromSmiles('C1CC1Oc1ccccc1-c1ncc(OC)cc1')
>>> res = Recap.RecapDecompose(m)
>>> res
<...Chem.Recap.RecapHierarchyNode object at ...>
>>> sorted(res.children.keys())
['*C1CC1', '*c1ccc(OC)cn1', '*c1ccccc1-c1ccc(OC)cn1', '*c1ccccc1OC1CC1']
>>> sorted(res.GetAllChildren().keys())
['*C1CC1', '*c1ccc(OC)cn1', '*c1ccccc1*', '*c1ccccc1-c1ccc(OC)cn1', '*c1ccccc1OC1CC1']

To get the standard set of RECAP results, use GetLeaves():
>>> leaves=res.GetLeaves()
>>> sorted(leaves.keys())
['*C1CC1', '*c1ccc(OC)cn1', '*c1ccccc1*']
>>> leaf = leaves['*C1CC1']
>>> leaf.mol
<...Chem.rdchem.Mol object at ...>


"""
import sys
import weakref

from rdkit import Chem
from rdkit.Chem import rdChemReactions as Reactions

# These are the definitions that will be applied to fragment molecules:
reactionDefs = (
  "[#7;+0;D2,D3:1]!@C(!@=O)!@[#7;+0;D2,D3:2]>>*[#7:1].[#7:2]*",  # urea
  "[C;!$(C([#7])[#7]):1](=!@[O:2])!@[#7;+0;!D1:3]>>*[C:1]=[O:2].*[#7:3]",  # amide
  "[C:1](=!@[O:2])!@[O;+0:3]>>*[C:1]=[O:2].[O:3]*",  # ester
  "[N;!D1;+0;!$(N-C=[#7,#8,#15,#16])](-!@[*:1])-!@[*:2]>>*[*:1].[*:2]*",  # amines
  # "[N;!D1](!@[*:1])!@[*:2]>>*[*:1].[*:2]*", # amines

  # again: what about aromatics?
  "[#7;R;D3;+0:1]-!@[*:2]>>*[#7:1].[*:2]*",  # cyclic amines
  "[#6:1]-!@[O;+0]-!@[#6:2]>>[#6:1]*.*[#6:2]",  # ether
  "[C:1]=!@[C:2]>>[C:1]*.*[C:2]",  # olefin
  "[n;+0:1]-!@[C:2]>>[n:1]*.[C:2]*",  # aromatic nitrogen - aliphatic carbon
  "[O:3]=[C:4]-@[N;+0:1]-!@[C:2]>>[O:3]=[C:4]-[N:1]*.[C:2]*",  # lactam nitrogen - aliphatic carbon
  "[c:1]-!@[c:2]>>[c:1]*.*[c:2]",  # aromatic carbon - aromatic carbon
  # aromatic nitrogen - aromatic carbon *NOTE* this is not part of the standard recap set.
  "[n;+0:1]-!@[c:2]>>[n:1]*.*[c:2]",
  "[#7;+0;D2,D3:1]-!@[S:2](=[O:3])=[O:4]>>[#7:1]*.*[S:2](=[O:3])=[O:4]",  # sulphonamide
)

reactions = tuple([Reactions.ReactionFromSmarts(x) for x in reactionDefs])


class RecapHierarchyNode(object):
  """ This class is used to hold the Recap hiearchy
    """
  mol = None
  children = None
  parents = None
  smiles = None

  def __init__(self, mol):
    self.mol = mol
    self.children = {}
    self.parents = {}

  def GetAllChildren(self):
    " returns a dictionary, keyed by SMILES, of children "
    res = {}
    for smi, child in self.children.items():
      res[smi] = child
      child._gacRecurse(res, terminalOnly=False)
    return res

  def GetLeaves(self):
    " returns a dictionary, keyed by SMILES, of leaf (terminal) nodes "
    res = {}
    for smi, child in self.children.items():
      if not len(child.children):
        res[smi] = child
      else:
        child._gacRecurse(res, terminalOnly=True)
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

  def _gacRecurse(self, res, terminalOnly=False):
    for smi, child in self.children.items():
      if not terminalOnly or not len(child.children):
        res[smi] = child
      child._gacRecurse(res, terminalOnly=terminalOnly)

  def __del__(self):
    self.children = {}
    self.parents = {}
    self.mol = None


def RecapDecompose(mol, allNodes=None, minFragmentSize=0, onlyUseReactions=None):
  """ returns the recap decomposition for a molecule """
  mSmi = Chem.MolToSmiles(mol, 1)

  if allNodes is None:
    allNodes = {}
  if mSmi in allNodes:
    return allNodes[mSmi]

  res = RecapHierarchyNode(mol)
  res.smiles = mSmi
  activePool = {mSmi: res}
  allNodes[mSmi] = res
  while activePool:
    nSmi = next(iter(activePool))
    node = activePool.pop(nSmi)
    if not node.mol:
      continue
    for rxnIdx, reaction in enumerate(reactions):
      if onlyUseReactions and rxnIdx not in onlyUseReactions:
        continue
      # print '  .',nSmi
      # print '         !!!!',rxnIdx,nSmi,reactionDefs[rxnIdx]
      ps = reaction.RunReactants((node.mol, ))
      # print '    ',len(ps)
      if ps:
        for prodSeq in ps:
          seqOk = True
          # we want to disqualify small fragments, so sort the product sequence by size
          # and then look for "forbidden" fragments
          tSeq = [(prod.GetNumAtoms(onlyExplicit=True), idx) for idx, prod in enumerate(prodSeq)]
          tSeq.sort()
          ts = [(x, prodSeq[y]) for x, y in tSeq]
          prodSeq = ts
          for nats, prod in prodSeq:
            try:
              Chem.SanitizeMol(prod)
            except Exception:
              continue
            pSmi = Chem.MolToSmiles(prod, 1)
            if minFragmentSize > 0:
              nDummies = pSmi.count('*')
              if nats - nDummies < minFragmentSize:
                seqOk = False
                break
            # don't forget after replacing dummy atoms to remove any empty
            # branches:
            elif pSmi.replace('*', '').replace('()', '') in ('', 'C', 'CC', 'CCC'):
              seqOk = False
              break
            prod.pSmi = pSmi
          if seqOk:
            for nats, prod in prodSeq:
              pSmi = prod.pSmi
              # print '\t',nats,pSmi
              if pSmi not in allNodes:
                pNode = RecapHierarchyNode(prod)
                pNode.smiles = pSmi
                pNode.parents[nSmi] = weakref.proxy(node)
                node.children[pSmi] = pNode
                activePool[pSmi] = pNode
                allNodes[pSmi] = pNode
              else:
                pNode = allNodes[pSmi]
                pNode.parents[nSmi] = weakref.proxy(node)
                node.children[pSmi] = pNode
            # print '                >>an:',allNodes.keys()
  return res


# ------- ------- ------- ------- ------- ------- ------- -------
# Begin testing code
import unittest


class TestCase(unittest.TestCase):

  def test1(self):
    m = Chem.MolFromSmiles('C1CC1Oc1ccccc1-c1ncc(OC)cc1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.children.keys()) == 4)
    self.assertTrue(len(res.GetAllChildren().keys()) == 5)
    self.assertTrue(len(res.GetLeaves().keys()) == 3)

  def test2(self):
    m = Chem.MolFromSmiles('CCCOCCC')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(res.children == {})

  def test3(self):
    allNodes = {}
    m = Chem.MolFromSmiles('c1ccccc1-c1ncccc1')
    res = RecapDecompose(m, allNodes=allNodes)
    self.assertTrue(res)
    self.assertTrue(len(res.children.keys()) == 2)
    self.assertTrue(len(allNodes.keys()) == 3)

    m = Chem.MolFromSmiles('COc1ccccc1-c1ncccc1')
    res = RecapDecompose(m, allNodes=allNodes)
    self.assertTrue(res)
    self.assertTrue(len(res.children.keys()) == 2)
    # we get two more nodes from that:
    self.assertTrue(len(allNodes.keys()) == 5)
    self.assertTrue('*c1ccccc1OC' in allNodes)
    self.assertTrue('*c1ccccc1' in allNodes)

    m = Chem.MolFromSmiles('C1CC1Oc1ccccc1-c1ncccc1')
    res = RecapDecompose(m, allNodes=allNodes)
    self.assertTrue(res)
    self.assertTrue(len(res.children.keys()) == 4)
    self.assertTrue(len(allNodes.keys()) == 10)

  def testSFNetIssue1801871(self):
    m = Chem.MolFromSmiles('c1ccccc1OC(Oc1ccccc1)Oc1ccccc1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertFalse('*C(*)*' in ks)
    self.assertTrue('*c1ccccc1' in ks)
    self.assertTrue('*C(*)Oc1ccccc1' in ks)

  def testSFNetIssue1804418(self):
    m = Chem.MolFromSmiles('C1CCCCN1CCCC')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*N1CCCCC1' in ks)
    self.assertTrue('*CCCC' in ks)

  def testMinFragmentSize(self):
    m = Chem.MolFromSmiles('CCCOCCC')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(res.children == {})
    res = RecapDecompose(m, minFragmentSize=3)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 1)
    ks = res.GetLeaves().keys()
    self.assertTrue('*CCC' in ks)

    m = Chem.MolFromSmiles('CCCOCC')
    res = RecapDecompose(m, minFragmentSize=3)
    self.assertTrue(res)
    self.assertTrue(res.children == {})

    m = Chem.MolFromSmiles('CCCOCCOC')
    res = RecapDecompose(m, minFragmentSize=2)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*CCC' in ks)
    ks = res.GetLeaves().keys()
    self.assertTrue('*CCOC' in ks)

  def testAmideRxn(self):
    m = Chem.MolFromSmiles('C1CC1C(=O)NC1OC1')
    res = RecapDecompose(m, onlyUseReactions=[1])
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*C(=O)C1CC1' in ks)
    self.assertTrue('*NC1CO1' in ks)

    m = Chem.MolFromSmiles('C1CC1C(=O)N(C)C1OC1')
    res = RecapDecompose(m, onlyUseReactions=[1])
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*C(=O)C1CC1' in ks)
    self.assertTrue('*N(C)C1CO1' in ks)

    m = Chem.MolFromSmiles('C1CC1C(=O)n1cccc1')
    res = RecapDecompose(m, onlyUseReactions=[1])
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*C(=O)C1CC1' in ks)
    self.assertTrue('*n1cccc1' in ks)

    m = Chem.MolFromSmiles('C1CC1C(=O)CC1OC1')
    res = RecapDecompose(m, onlyUseReactions=[1])
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 0)

    m = Chem.MolFromSmiles('C1CCC(=O)NC1')
    res = RecapDecompose(m, onlyUseReactions=[1])
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 0)

    m = Chem.MolFromSmiles('CC(=O)NC')
    res = RecapDecompose(m, onlyUseReactions=[1])
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()

    m = Chem.MolFromSmiles('CC(=O)N')
    res = RecapDecompose(m, onlyUseReactions=[1])
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 0)

    m = Chem.MolFromSmiles('C(=O)NCCNC(=O)CC')
    res = RecapDecompose(m, onlyUseReactions=[1])
    self.assertTrue(res)
    self.assertTrue(len(res.children) == 4)
    self.assertTrue(len(res.GetLeaves()) == 3)

  def testEsterRxn(self):
    m = Chem.MolFromSmiles('C1CC1C(=O)OC1OC1')
    res = RecapDecompose(m, onlyUseReactions=[2])
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*C(=O)C1CC1' in ks)
    self.assertTrue('*OC1CO1' in ks)

    m = Chem.MolFromSmiles('C1CC1C(=O)CC1OC1')
    res = RecapDecompose(m, onlyUseReactions=[2])
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 0)

    m = Chem.MolFromSmiles('C1CCC(=O)OC1')
    res = RecapDecompose(m, onlyUseReactions=[2])
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 0)

  def testUreaRxn(self):
    m = Chem.MolFromSmiles('C1CC1NC(=O)NC1OC1')
    res = RecapDecompose(m, onlyUseReactions=[0])
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*NC1CC1' in ks)
    self.assertTrue('*NC1CO1' in ks)

    m = Chem.MolFromSmiles('C1CC1NC(=O)N(C)C1OC1')
    res = RecapDecompose(m, onlyUseReactions=[0])
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*NC1CC1' in ks)
    self.assertTrue('*N(C)C1CO1' in ks)

    m = Chem.MolFromSmiles('C1CCNC(=O)NC1C')
    res = RecapDecompose(m, onlyUseReactions=[0])
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 0)

    m = Chem.MolFromSmiles('c1cccn1C(=O)NC1OC1')
    res = RecapDecompose(m, onlyUseReactions=[0])
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*n1cccc1' in ks)
    self.assertTrue('*NC1CO1' in ks)

    m = Chem.MolFromSmiles('c1cccn1C(=O)n1c(C)ccc1')
    res = RecapDecompose(m, onlyUseReactions=[0])
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*n1cccc1C' in ks)

  def testAmineRxn(self):
    m = Chem.MolFromSmiles('C1CC1N(C1NC1)C1OC1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 3)
    ks = res.GetLeaves().keys()
    self.assertTrue('*C1CC1' in ks)
    self.assertTrue('*C1CO1' in ks)
    self.assertTrue('*C1CN1' in ks)

    m = Chem.MolFromSmiles('c1ccccc1N(C1NC1)C1OC1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 3)
    ks = res.GetLeaves().keys()
    self.assertTrue('*c1ccccc1' in ks)
    self.assertTrue('*C1CO1' in ks)
    self.assertTrue('*C1CN1' in ks)

    m = Chem.MolFromSmiles('c1ccccc1N(c1ncccc1)C1OC1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 3)
    ks = res.GetLeaves().keys()
    self.assertTrue('*c1ccccc1' in ks)
    self.assertTrue('*c1ccccn1' in ks)
    self.assertTrue('*C1CO1' in ks)

    m = Chem.MolFromSmiles('c1ccccc1N(c1ncccc1)c1ccco1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 3)
    ks = res.GetLeaves().keys()
    self.assertTrue('*c1ccccc1' in ks)
    self.assertTrue('*c1ccccn1' in ks)
    self.assertTrue('*c1ccco1' in ks)

    m = Chem.MolFromSmiles('C1CCCCN1C1CC1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*N1CCCCC1' in ks)
    self.assertTrue('*C1CC1' in ks)

    m = Chem.MolFromSmiles('C1CCC2N1CC2')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 0)

  def testEtherRxn(self):
    m = Chem.MolFromSmiles('C1CC1OC1OC1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*C1CC1' in ks)
    self.assertTrue('*C1CO1' in ks)

    m = Chem.MolFromSmiles('C1CCCCO1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 0)

    m = Chem.MolFromSmiles('c1ccccc1OC1OC1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*c1ccccc1' in ks)
    self.assertTrue('*C1CO1' in ks)

    m = Chem.MolFromSmiles('c1ccccc1Oc1ncccc1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*c1ccccc1' in ks)
    self.assertTrue('*c1ccccn1' in ks)

  def testOlefinRxn(self):
    m = Chem.MolFromSmiles('ClC=CBr')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*CCl' in ks)
    self.assertTrue('*CBr' in ks)

    m = Chem.MolFromSmiles('C1CC=CC1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 0)

  def testAromNAliphCRxn(self):
    m = Chem.MolFromSmiles('c1cccn1CCCC')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*n1cccc1' in ks)
    self.assertTrue('*CCCC' in ks)

    m = Chem.MolFromSmiles('c1ccc2n1CCCC2')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 0)

  def testLactamNAliphCRxn(self):
    m = Chem.MolFromSmiles('C1CC(=O)N1CCCC')
    res = RecapDecompose(m, onlyUseReactions=[8])
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*N1CCC1=O' in ks)
    self.assertTrue('*CCCC' in ks)

    m = Chem.MolFromSmiles('O=C1CC2N1CCCC2')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 0)

  def testAromCAromCRxn(self):
    m = Chem.MolFromSmiles('c1ccccc1c1ncccc1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*c1ccccc1' in ks)
    self.assertTrue('*c1ccccn1' in ks)

    m = Chem.MolFromSmiles('c1ccccc1C1CC1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 0)

  def testAromNAromCRxn(self):
    m = Chem.MolFromSmiles('c1cccn1c1ccccc1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*n1cccc1' in ks)
    self.assertTrue('*c1ccccc1' in ks)

  def testSulfonamideRxn(self):
    m = Chem.MolFromSmiles('CCCNS(=O)(=O)CC')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*NCCC' in ks)
    self.assertTrue('*S(=O)(=O)CC' in ks)

    m = Chem.MolFromSmiles('c1cccn1S(=O)(=O)CC')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    ks = res.GetLeaves().keys()
    self.assertTrue('*n1cccc1' in ks)
    self.assertTrue('*S(=O)(=O)CC' in ks)

    m = Chem.MolFromSmiles('C1CNS(=O)(=O)CC1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 0)

  def testSFNetIssue1881803(self):
    m = Chem.MolFromSmiles('c1ccccc1n1cccc1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    m = Chem.MolFromSmiles('c1ccccc1[n+]1ccccc1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 0)

    m = Chem.MolFromSmiles('C1CC1NC(=O)CC')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    m = Chem.MolFromSmiles('C1CC1[NH+]C(=O)CC')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 0)

    m = Chem.MolFromSmiles('C1CC1NC(=O)NC1CCC1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 2)
    m = Chem.MolFromSmiles('C1CC1[NH+]C(=O)[NH+]C1CCC1')
    res = RecapDecompose(m)
    self.assertTrue(res)
    self.assertTrue(len(res.GetLeaves()) == 0)


if __name__ == '__main__':
  unittest.main()
