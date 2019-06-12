#
#  Copyright (c) 2009, Novartis Institutes for BioMedical Research Inc.
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
# Created by Greg Landrum, Nov 2008
import random
import copy
""" Implementation of the BRICS algorithm from Degen et al. ChemMedChem *3* 1503-7 (2008)

"""
import sys
import re
import random
from rdkit import Chem
from rdkit.Chem import rdChemReactions as Reactions

# These are the definitions that will be applied to fragment molecules:
environs = {
  'L1': '[C;D3]([#0,#6,#7,#8])(=O)',
  #
  # After some discussion, the L2 definitions ("N.pl3" in the original
  # paper) have been removed and incorporated into a (almost) general
  # purpose amine definition in L5 ("N.sp3" in the paper).
  #
  # The problem is one of consistency.
  #    Based on the original definitions you should get the following
  #    fragmentations:
  #      C1CCCCC1NC(=O)C -> C1CCCCC1N[2*].[1*]C(=O)C
  #      c1ccccc1NC(=O)C -> c1ccccc1[16*].[2*]N[2*].[1*]C(=O)C
  #    This difference just didn't make sense to us. By switching to
  #    the unified definition we end up with:
  #      C1CCCCC1NC(=O)C -> C1CCCCC1[15*].[5*]N[5*].[1*]C(=O)C
  #      c1ccccc1NC(=O)C -> c1ccccc1[16*].[5*]N[5*].[1*]C(=O)C
  #
  # 'L2':'[N;!R;!D1;!$(N=*)]-;!@[#0,#6]',
  # this one turned out to be too tricky to define above, so we set it off
  # in its own definition:
  # 'L2a':'[N;D3;R;$(N(@[C;!$(C=*)])@[C;!$(C=*)])]',
  'L3': '[O;D2]-;!@[#0,#6,#1]',
  'L4': '[C;!D1;!$(C=*)]-;!@[#6]',
  # 'L5':'[N;!D1;!$(N*!-*);!$(N=*);!$(N-[!C;!#0])]-[#0,C]',
  'L5': '[N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]',
  'L6': '[C;D3;!R](=O)-;!@[#0,#6,#7,#8]',
  'L7a': '[C;D2,D3]-[#6]',
  'L7b': '[C;D2,D3]-[#6]',
  '#L8': '[C;!R;!D1]-;!@[#6]',
  'L8': '[C;!R;!D1;!$(C!-*)]',
  'L9': '[n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]',
  'L10': '[N;R;$(N(@C(=O))@[C,N,O,S])]',
  'L11': '[S;D2](-;!@[#0,#6])',
  'L12': '[S;D4]([#6,#0])(=O)(=O)',
  'L13': '[C;$(C(-;@[C,N,O,S])-;@[N,O,S])]',
  'L14': '[c;$(c(:[c,n,o,s]):[n,o,s])]',
  'L14b': '[c;$(c(:[c,n,o,s]):[n,o,s])]',
  'L15': '[C;$(C(-;@C)-;@C)]',
  'L16': '[c;$(c(:c):c)]',
  'L16b': '[c;$(c(:c):c)]',
}
reactionDefs = (
  # L1
  [
    ('1', '3', '-'),
    ('1', '5', '-'),
    ('1', '10', '-'),
  ],

  # L3
  [
    ('3', '4', '-'),
    ('3', '13', '-'),
    ('3', '14', '-'),
    ('3', '15', '-'),
    ('3', '16', '-'),
  ],

  # L4
  [
    ('4', '5', '-'),
    ('4', '11', '-'),
  ],

  # L5
  [
    ('5', '12', '-'),
    ('5', '14', '-'),
    ('5', '16', '-'),
    ('5', '13', '-'),
    ('5', '15', '-'),
  ],

  # L6
  [
    ('6', '13', '-'),
    ('6', '14', '-'),
    ('6', '15', '-'),
    ('6', '16', '-'),
  ],

  # L7
  [
    ('7a', '7b', '='),
  ],

  # L8
  [
    ('8', '9', '-'),
    ('8', '10', '-'),
    ('8', '13', '-'),
    ('8', '14', '-'),
    ('8', '15', '-'),
    ('8', '16', '-'),
  ],

  # L9
  [
    ('9', '13', '-'),  # not in original paper
    ('9', '14', '-'),  # not in original paper
    ('9', '15', '-'),
    ('9', '16', '-'),
  ],

  # L10
  [
    ('10', '13', '-'),
    ('10', '14', '-'),
    ('10', '15', '-'),
    ('10', '16', '-'),
  ],

  # L11
  [
    ('11', '13', '-'),
    ('11', '14', '-'),
    ('11', '15', '-'),
    ('11', '16', '-'),
  ],

  # L12
  # none left

  # L13
  [
    ('13', '14', '-'),
    ('13', '15', '-'),
    ('13', '16', '-'),
  ],

  # L14
  [
    ('14', '14', '-'),  # not in original paper
    ('14', '15', '-'),
    ('14', '16', '-'),
  ],

  # L15
  [
    ('15', '16', '-'),
  ],

  # L16
  [
    ('16', '16', '-'),  # not in original paper
  ], )
smartsGps = copy.deepcopy(reactionDefs)
for gp in smartsGps:
    for j, defn in enumerate(gp):
        g1, g2, bnd = defn
        r1 = environs['L' + g1]
        r2 = environs['L' + g2]
        g1 = re.sub('[a-z,A-Z]', '', g1)
        g2 = re.sub('[a-z,A-Z]', '', g2)
        sma = '[$(%s):1]%s;!@[$(%s):2]>>[%s*]-[*:1].[%s*]-[*:2]' % (r1, bnd, r2, g1, g2)
        gp[j] = sma

for gp in smartsGps:
    for defn in gp:
        try:
            t = Reactions.ReactionFromSmarts(defn)
            t.Initialize()
        except Exception:
            print(defn)
            raise

environMatchers = {}
for env, sma in environs.items():
    environMatchers[env] = Chem.MolFromSmarts(sma)

bondMatchers = []
for i, compats in enumerate(reactionDefs):
    tmp = []
    for i1, i2, bType in compats:
        e1 = environs['L%s' % i1]
        e2 = environs['L%s' % i2]
        patt = '[$(%s)]%s;!@[$(%s)]' % (e1, bType, e2)
        patt = Chem.MolFromSmarts(patt)
        tmp.append((i1, i2, bType, patt))
    bondMatchers.append(tmp)

reactions = tuple([[Reactions.ReactionFromSmarts(y) for y in x] for x in smartsGps])
reverseReactions = []
for i, rxnSet in enumerate(smartsGps):
    for j, sma in enumerate(rxnSet):
        rs, ps = sma.split('>>')
        sma = '%s>>%s' % (ps, rs)
        rxn = Reactions.ReactionFromSmarts(sma)
        labels = re.findall(r'\[([0-9]+?)\*\]', ps)
        rxn._matchers = [Chem.MolFromSmiles('[%s*]' % x) for x in labels]
        reverseReactions.append(rxn)


def FindBRICSBonds(mol, randomizeOrder=False, silent=True):
    """ returns the bonds in a molecule that BRICS would cleave

    >>> from rdkit import Chem
    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> res = list(FindBRICSBonds(m))
    >>> res
    [((3, 2), ('3', '4')), ((3, 4), ('3', '4'))]

    a more complicated case:

    >>> m = Chem.MolFromSmiles('CCCOCCC(=O)c1ccccc1')
    >>> res = list(FindBRICSBonds(m))
    >>> res
    [((3, 2), ('3', '4')), ((3, 4), ('3', '4')), ((6, 8), ('6', '16'))]

    we can also randomize the order of the results:

    >>> random.seed(23)
    >>> res = list(FindBRICSBonds(m,randomizeOrder=True))
    >>> sorted(res)
    [((3, 2), ('3', '4')), ((3, 4), ('3', '4')), ((6, 8), ('6', '16'))]

    Note that this is a generator function :

    >>> res = FindBRICSBonds(m)
    >>> res
    <generator object ...>
    >>> next(res)
    ((3, 2), ('3', '4'))

    >>> m = Chem.MolFromSmiles('CC=CC')
    >>> res = list(FindBRICSBonds(m))
    >>> sorted(res)
    [((1, 2), ('7', '7'))]

    make sure we don't match ring bonds:

    >>> m = Chem.MolFromSmiles('O=C1NCCC1')
    >>> list(FindBRICSBonds(m))
    []

    another nice one, make sure environment 8 doesn't match something connected
    to a ring atom:

    >>> m = Chem.MolFromSmiles('CC1(C)CCCCC1')
    >>> list(FindBRICSBonds(m))
    []

    """
    letter = re.compile('[a-z,A-Z]')
    indices = list(range(len(bondMatchers)))
    bondsDone = set()
    if randomizeOrder:
        random.shuffle(indices, random=random.random)

    envMatches = {}
    for env, patt in environMatchers.items():
        envMatches[env] = mol.HasSubstructMatch(patt)
    for gpIdx in indices:
        if randomizeOrder:
            compats = bondMatchers[gpIdx][:]
            random.shuffle(compats, random=random.random)
        else:
            compats = bondMatchers[gpIdx]
        for i1, i2, bType, patt in compats:
            if not envMatches['L' + i1] or not envMatches['L' + i2]:
                continue
            matches = mol.GetSubstructMatches(patt)
            i1 = letter.sub('', i1)
            i2 = letter.sub('', i2)
            for match in matches:
                if match not in bondsDone and (match[1], match[0]) not in bondsDone:
                    bondsDone.add(match)
                    yield (((match[0], match[1]), (i1, i2)))


def BreakBRICSBonds(mol, bonds=None, sanitize=True, silent=True):
    """ breaks the BRICS bonds in a molecule and returns the results

    >>> from rdkit import Chem
    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> m2=BreakBRICSBonds(m)
    >>> Chem.MolToSmiles(m2,True)
    '[3*]O[3*].[4*]CC.[4*]CCC'

    a more complicated case:

    >>> m = Chem.MolFromSmiles('CCCOCCC(=O)c1ccccc1')
    >>> m2=BreakBRICSBonds(m)
    >>> Chem.MolToSmiles(m2,True)
    '[16*]c1ccccc1.[3*]O[3*].[4*]CCC.[4*]CCC([6*])=O'


    can also specify a limited set of bonds to work with:

    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> m2 = BreakBRICSBonds(m,[((3, 2), ('3', '4'))])
    >>> Chem.MolToSmiles(m2,True)
    '[3*]OCC.[4*]CCC'

    this can be used as an alternate approach for doing a BRICS decomposition by
    following BreakBRICSBonds with a call to Chem.GetMolFrags:

    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> m2=BreakBRICSBonds(m)
    >>> frags = Chem.GetMolFrags(m2,asMols=True)
    >>> [Chem.MolToSmiles(x,True) for x in frags]
    ['[4*]CCC', '[3*]O[3*]', '[4*]CC']

    """
    if not bonds:
        #bonds = FindBRICSBonds(mol)
        res = Chem.FragmentOnBRICSBonds(mol)
        if sanitize:
            Chem.SanitizeMol(res)
        return res
    eMol = Chem.EditableMol(mol)
    nAts = mol.GetNumAtoms()

    dummyPositions = []
    for indices, dummyTypes in bonds:
        ia, ib = indices
        obond = mol.GetBondBetweenAtoms(ia, ib)
        bondType = obond.GetBondType()
        eMol.RemoveBond(ia, ib)

        da, db = dummyTypes
        atoma = Chem.Atom(0)
        atoma.SetIsotope(int(da))
        atoma.SetNoImplicit(True)
        idxa = nAts
        nAts += 1
        eMol.AddAtom(atoma)
        eMol.AddBond(ia, idxa, bondType)

        atomb = Chem.Atom(0)
        atomb.SetIsotope(int(db))
        atomb.SetNoImplicit(True)
        idxb = nAts
        nAts += 1
        eMol.AddAtom(atomb)
        eMol.AddBond(ib, idxb, bondType)
        if mol.GetNumConformers():
            dummyPositions.append((idxa, ib))
            dummyPositions.append((idxb, ia))
    res = eMol.GetMol()
    if sanitize:
        Chem.SanitizeMol(res)
    if mol.GetNumConformers():
        for conf in mol.GetConformers():
            resConf = res.GetConformer(conf.GetId())
            for ia, pa in dummyPositions:
                resConf.SetAtomPosition(ia, conf.GetAtomPosition(pa))
    return res


def BRICSDecompose(mol, allNodes=None, minFragmentSize=1, onlyUseReactions=None, silent=True,
                   keepNonLeafNodes=False, singlePass=False, returnMols=False):
    """ returns the BRICS decomposition for a molecule

    >>> from rdkit import Chem
    >>> m = Chem.MolFromSmiles('CCCOCc1cc(c2ncccc2)ccc1')
    >>> res = list(BRICSDecompose(m))
    >>> sorted(res)
    ['[14*]c1ccccn1', '[16*]c1cccc([16*])c1', '[3*]O[3*]', '[4*]CCC', '[4*]C[8*]']

    >>> res = list(BRICSDecompose(m,returnMols=True))
    >>> res[0]
    <rdkit.Chem.rdchem.Mol object ...>
    >>> smis = [Chem.MolToSmiles(x,True) for x in res]
    >>> sorted(smis)
    ['[14*]c1ccccn1', '[16*]c1cccc([16*])c1', '[3*]O[3*]', '[4*]CCC', '[4*]C[8*]']

    nexavar, an example from the paper (corrected):

    >>> m = Chem.MolFromSmiles('CNC(=O)C1=NC=CC(OC2=CC=C(NC(=O)NC3=CC(=C(Cl)C=C3)C(F)(F)F)C=C2)=C1')
    >>> res = list(BRICSDecompose(m))
    >>> sorted(res)
    ['[1*]C([1*])=O', '[1*]C([6*])=O', '[14*]c1cc([16*])ccn1', '[16*]c1ccc(Cl)c([16*])c1', '[16*]c1ccc([16*])cc1', '[3*]O[3*]', '[5*]NC', '[5*]N[5*]', '[8*]C(F)(F)F']

    it's also possible to keep pieces that haven't been fully decomposed:

    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> res = list(BRICSDecompose(m,keepNonLeafNodes=True))
    >>> sorted(res)
    ['CCCOCC', '[3*]OCC', '[3*]OCCC', '[3*]O[3*]', '[4*]CC', '[4*]CCC']

    >>> m = Chem.MolFromSmiles('CCCOCc1cc(c2ncccc2)ccc1')
    >>> res = list(BRICSDecompose(m,keepNonLeafNodes=True))
    >>> sorted(res)
    ['CCCOCc1cccc(-c2ccccn2)c1', '[14*]c1ccccn1', '[16*]c1cccc(-c2ccccn2)c1', '[16*]c1cccc(COCCC)c1', '[16*]c1cccc([16*])c1', '[3*]OCCC', '[3*]OC[8*]', '[3*]OCc1cccc(-c2ccccn2)c1', '[3*]OCc1cccc([16*])c1', '[3*]O[3*]', '[4*]CCC', '[4*]C[8*]', '[4*]Cc1cccc(-c2ccccn2)c1', '[4*]Cc1cccc([16*])c1', '[8*]COCCC']

    or to only do a single pass of decomposition:

    >>> m = Chem.MolFromSmiles('CCCOCc1cc(c2ncccc2)ccc1')
    >>> res = list(BRICSDecompose(m,singlePass=True))
    >>> sorted(res)
    ['CCCOCc1cccc(-c2ccccn2)c1', '[14*]c1ccccn1', '[16*]c1cccc(-c2ccccn2)c1', '[16*]c1cccc(COCCC)c1', '[3*]OCCC', '[3*]OCc1cccc(-c2ccccn2)c1', '[4*]CCC', '[4*]Cc1cccc(-c2ccccn2)c1', '[8*]COCCC']

    setting a minimum size for the fragments:

    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> res = list(BRICSDecompose(m,keepNonLeafNodes=True,minFragmentSize=2))
    >>> sorted(res)
    ['CCCOCC', '[3*]OCC', '[3*]OCCC', '[4*]CC', '[4*]CCC']
    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> res = list(BRICSDecompose(m,keepNonLeafNodes=True,minFragmentSize=3))
    >>> sorted(res)
    ['CCCOCC', '[3*]OCC', '[4*]CCC']
    >>> res = list(BRICSDecompose(m,minFragmentSize=2))
    >>> sorted(res)
    ['[3*]OCC', '[3*]OCCC', '[4*]CC', '[4*]CCC']


    """
    global reactions
    mSmi = Chem.MolToSmiles(mol, 1)

    if allNodes is None:
        allNodes = set()

    if mSmi in allNodes:
        return set()

    activePool = {mSmi: mol}
    allNodes.add(mSmi)
    foundMols = {mSmi: mol}
    for gpIdx, reactionGp in enumerate(reactions):
        newPool = {}
        while activePool:
            matched = False
            nSmi = next(iter(activePool))
            mol = activePool.pop(nSmi)
            for rxnIdx, reaction in enumerate(reactionGp):
                if onlyUseReactions and (gpIdx, rxnIdx) not in onlyUseReactions:
                    continue
                if not silent:
                    print('--------')
                    print(smartsGps[gpIdx][rxnIdx])
                ps = reaction.RunReactants((mol, ))
                if ps:
                    if not silent:
                        print(nSmi, '->', len(ps), 'products')
                    for prodSeq in ps:
                        seqOk = True
                        # we want to disqualify small fragments, so sort the product sequence by size
                        tSeq = [(prod.GetNumAtoms(onlyExplicit=True), idx)
                                for idx, prod in enumerate(prodSeq)]
                        tSeq.sort()
                        for nats, idx in tSeq:
                            prod = prodSeq[idx]
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
                            prod.pSmi = pSmi
                        ts = [(x, prodSeq[y]) for x, y in tSeq]
                        prodSeq = ts
                        if seqOk:
                            matched = True
                            for nats, prod in prodSeq:
                                pSmi = prod.pSmi
                                # print('\t',nats,pSmi)
                                if pSmi not in allNodes:
                                    if not singlePass:
                                        activePool[pSmi] = prod
                                    allNodes.add(pSmi)
                                    foundMols[pSmi] = prod
            if singlePass or keepNonLeafNodes or not matched:
                newPool[nSmi] = mol
        activePool = newPool
    if not (singlePass or keepNonLeafNodes):
        if not returnMols:
            res = set(activePool.keys())
        else:
            res = activePool.values()
    else:
        if not returnMols:
            res = allNodes
        else:
            res = foundMols.values()
    return res


dummyPattern = Chem.MolFromSmiles('[*]')


def BRICSBuild(fragments, onlyCompleteMols=True, seeds=None, uniquify=True, scrambleReagents=True,
               maxDepth=3):
    seen = set()
    if not seeds:
        seeds = list(fragments)
    if scrambleReagents:
        seeds = list(seeds)
        random.shuffle(seeds, random=random.random)
    if scrambleReagents:
        tempReactions = list(reverseReactions)
        random.shuffle(tempReactions, random=random.random)
    else:
        tempReactions = reverseReactions
    for seed in seeds:
        seedIsR1 = False
        seedIsR2 = False
        nextSteps = []
        for rxn in tempReactions:
            if seed.HasSubstructMatch(rxn._matchers[0]):
                seedIsR1 = True
            if seed.HasSubstructMatch(rxn._matchers[1]):
                seedIsR2 = True
            for fragment in fragments:
                ps = None
                if fragment.HasSubstructMatch(rxn._matchers[0]):
                    if seedIsR2:
                        ps = rxn.RunReactants((fragment, seed))
                if fragment.HasSubstructMatch(rxn._matchers[1]):
                    if seedIsR1:
                        ps = rxn.RunReactants((seed, fragment))
                if ps:
                    for p in ps:
                        if uniquify:
                            pSmi = Chem.MolToSmiles(p[0], True)
                            if pSmi in seen:
                                continue
                            else:
                                seen.add(pSmi)
                        if p[0].HasSubstructMatch(dummyPattern):
                            nextSteps.append(p[0])
                            if not onlyCompleteMols:
                                yield p[0]
                        else:
                            yield p[0]
        if nextSteps and maxDepth > 0:
            for p in BRICSBuild(fragments, onlyCompleteMols=onlyCompleteMols, seeds=nextSteps,
                                uniquify=uniquify, maxDepth=maxDepth - 1,
                                scrambleReagents=scrambleReagents):
                if uniquify:
                    pSmi = Chem.MolToSmiles(p, True)
                    if pSmi in seen:
                        continue
                    else:
                        seen.add(pSmi)
                yield p

    # ------- ------- ------- ------- ------- ------- ------- -------
    # Begin testing code

    # ------------------------------------
    #
    #  doctest boilerplate
    #


def _test():
    import doctest
    import sys
    return doctest.testmod(sys.modules["__main__"],
                           optionflags=doctest.ELLIPSIS + doctest.NORMALIZE_WHITESPACE)


if __name__ == '__main__':
    import unittest

    class TestCase(unittest.TestCase):

        def test1(self):
            m = Chem.MolFromSmiles('CC(=O)OC')
            res = BRICSDecompose(m)
            self.assertTrue(res)
            self.assertTrue(len(res) == 2)

            m = Chem.MolFromSmiles('CC(=O)N1CCC1=O')
            res = BRICSDecompose(m)
            self.assertTrue(res)
            self.assertTrue(len(res) == 2, res)

            m = Chem.MolFromSmiles('c1ccccc1N(C)C')
            res = BRICSDecompose(m)
            self.assertTrue(res)
            self.assertTrue(len(res) == 2, res)

            m = Chem.MolFromSmiles('c1cccnc1N(C)C')
            res = BRICSDecompose(m)
            self.assertTrue(res)
            self.assertTrue(len(res) == 2, res)

            m = Chem.MolFromSmiles('o1ccnc1N(C)C')
            res = BRICSDecompose(m)
            self.assertTrue(res)
            self.assertTrue(len(res) == 2)

            m = Chem.MolFromSmiles('c1ccccc1OC')
            res = BRICSDecompose(m)
            self.assertTrue(res)
            self.assertTrue(len(res) == 2)

            m = Chem.MolFromSmiles('o1ccnc1OC')
            res = BRICSDecompose(m)
            self.assertTrue(res)
            self.assertTrue(len(res) == 2)

            m = Chem.MolFromSmiles('O1CCNC1OC')
            res = BRICSDecompose(m)
            self.assertTrue(res)
            self.assertTrue(len(res) == 2)

            m = Chem.MolFromSmiles('CCCSCC')
            res = BRICSDecompose(m)
            self.assertTrue(res)
            self.assertTrue(len(res) == 3, res)
            self.assertTrue('[11*]S[11*]' in res, res)

            m = Chem.MolFromSmiles('CCNC(=O)C1CC1')
            res = BRICSDecompose(m)
            self.assertTrue(res)
            self.assertTrue(len(res) == 4, res)
            self.assertTrue('[5*]N[5*]' in res, res)

        def test2(self):
            # example from the paper, nexavar:
            m = Chem.MolFromSmiles(
                'CNC(=O)C1=NC=CC(OC2=CC=C(NC(=O)NC3=CC(=C(Cl)C=C3)C(F)(F)F)C=C2)=C1')
            res = BRICSDecompose(m)
            self.assertTrue(res)
            self.assertTrue(len(res) == 9, res)

        def test3(self):
            m = Chem.MolFromSmiles('FC(F)(F)C1=C(Cl)C=CC(NC(=O)NC2=CC=CC=C2)=C1')
            res = BRICSDecompose(m)
            self.assertTrue(res)
            self.assertTrue(len(res) == 5, res)
            self.assertTrue('[5*]N[5*]' in res, res)
            self.assertTrue('[16*]c1ccccc1' in res, res)
            self.assertTrue('[8*]C(F)(F)F' in res, res)

        def test4(self):
            allNodes = set()
            m = Chem.MolFromSmiles('c1ccccc1OCCC')
            res = BRICSDecompose(m, allNodes=allNodes)
            self.assertTrue(res)
            leaves = res
            self.assertTrue(len(leaves) == 3, leaves)
            self.assertTrue(len(allNodes) == 6, allNodes)
            res = BRICSDecompose(m, allNodes=allNodes)
            self.assertFalse(res)
            self.assertTrue(len(allNodes) == 6, allNodes)

            m = Chem.MolFromSmiles('c1ccccc1OCCCC')
            res = BRICSDecompose(m, allNodes=allNodes)
            self.assertTrue(res)
            leaves.update(res)
            self.assertTrue(len(allNodes) == 9, allNodes)
            self.assertTrue(len(leaves) == 4, leaves)

            m = Chem.MolFromSmiles('c1cc(C(=O)NCC)ccc1OCCC')
            res = BRICSDecompose(m, allNodes=allNodes)
            self.assertTrue(res)
            leaves.update(res)
            self.assertTrue(len(leaves) == 8, leaves)
            self.assertTrue(len(allNodes) == 18, allNodes)

        def test5(self):
            allNodes = set()
            frags = [
              '[14*]c1ncncn1',
              '[16*]c1ccccc1',
              '[14*]c1ncccc1',
            ]
            frags = [Chem.MolFromSmiles(x) for x in frags]
            res = BRICSBuild(frags)
            self.assertTrue(res)
            res = list(res)
            self.assertTrue(len(res) == 6)
            smis = [Chem.MolToSmiles(x, True) for x in res]
            self.assertTrue('c1ccc(-c2ccccc2)cc1' in smis)
            self.assertTrue('c1ccc(-c2ccccn2)cc1' in smis)

        def test5a(self):
            allNodes = set()
            frags = [
              '[3*]O[3*]',
              '[16*]c1ccccc1',
            ]
            frags = [Chem.MolFromSmiles(x) for x in frags]
            res = BRICSBuild(frags)
            self.assertTrue(res)
            res = list(res)
            smis = [Chem.MolToSmiles(x, True) for x in res]
            self.assertTrue(len(smis) == 2, smis)
            self.assertTrue('c1ccc(Oc2ccccc2)cc1' in smis)
            self.assertTrue('c1ccc(-c2ccccc2)cc1' in smis)

        def test6(self):
            allNodes = set()
            frags = [
              '[16*]c1ccccc1',
              '[3*]OC',
              '[9*]n1cccc1',
            ]
            frags = [Chem.MolFromSmiles(x) for x in frags]
            res = BRICSBuild(frags)
            self.assertTrue(res)
            res = list(res)
            self.assertTrue(len(res) == 3)
            smis = [Chem.MolToSmiles(x, True) for x in res]
            self.assertTrue('c1ccc(-c2ccccc2)cc1' in smis)
            self.assertTrue('COc1ccccc1' in smis)
            self.assertTrue('c1ccc(-n2cccc2)cc1' in smis, smis)

        def test7(self):
            allNodes = set()
            frags = [
              '[16*]c1ccccc1',
              '[3*]OC',
              '[3*]OCC(=O)[6*]',
            ]
            frags = [Chem.MolFromSmiles(x) for x in frags]
            res = BRICSBuild(frags)
            self.assertTrue(res)
            res = list(res)
            smis = [Chem.MolToSmiles(x, True) for x in res]
            self.assertTrue(len(res) == 3)
            self.assertTrue('c1ccc(-c2ccccc2)cc1' in smis)
            self.assertTrue('COc1ccccc1' in smis)
            self.assertTrue('O=C(COc1ccccc1)c1ccccc1' in smis)

        def test8(self):
            random.seed(23)
            base = Chem.MolFromSmiles("n1cncnc1OCC(C1CC1)OC1CNC1")
            catalog = BRICSDecompose(base)
            self.assertTrue(len(catalog) == 5, catalog)
            catalog = [Chem.MolFromSmiles(x) for x in catalog]
            ms = list(BRICSBuild(catalog, maxDepth=4, scrambleReagents=False))
            for m in ms:
                Chem.SanitizeMol(m)
            ms = [Chem.MolToSmiles(x) for x in ms]
            self.assertEqual(len(ms), 36)

            ts = ['n1cnc(C2CNC2)nc1', 'n1cnc(-c2ncncn2)nc1', 'C(OC1CNC1)C(C1CC1)OC1CNC1',
                  'n1cnc(OC(COC2CNC2)C2CC2)nc1', 'n1cnc(OCC(OC2CNC2)C2CNC2)nc1']
            ts = [Chem.MolToSmiles(Chem.MolFromSmiles(x), True) for x in ts]
            for t in ts:
                self.assertTrue(t in ms, (t, ms))

            ms2 = list(BRICSBuild(catalog, maxDepth=4, scrambleReagents=False))
            for m in ms2:
                Chem.SanitizeMol(m)
            ms2 = [Chem.MolToSmiles(x) for x in ms2]
            self.assertEqual(ms, ms2)

            ms2 = list(BRICSBuild(catalog, maxDepth=4, scrambleReagents=True))
            for m in ms2:
                Chem.SanitizeMol(m)
            ms2 = [Chem.MolToSmiles(x) for x in ms2]
            self.assertNotEqual(ms, ms2)

        def test9(self):
            m = Chem.MolFromSmiles('CCOc1ccccc1c1ncc(c2nc(NCCCC)ncn2)cc1')
            res = BRICSDecompose(m)
            self.assertEqual(len(res), 7)
            self.assertTrue('[3*]O[3*]' in res)
            self.assertFalse('[14*]c1ncnc(NCCCC)n1' in res)
            res = BRICSDecompose(m, singlePass=True)
            self.assertEqual(len(res), 13)
            self.assertTrue('[3*]OCC' in res)
            self.assertTrue('[14*]c1ncnc(NCCCC)n1' in res)

        def test10(self):
            m = Chem.MolFromSmiles('C1CCCCN1c1ccccc1')
            res = BRICSDecompose(m)
            self.assertEqual(len(res), 2, res)

        def test11(self):
            # test coordinate preservation:
            molblock = """
     RDKit          3D

 13 14  0  0  0  0  0  0  0  0999 V2000
   -1.2004    0.5900    0.6110 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2328    1.3173    0.0343 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4299    0.6533   -0.1500 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3633   -0.7217   -0.3299 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1552   -1.3791   -0.2207 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1425   -0.7969    0.5335 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1458   -1.4244    0.4108 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2976   -0.7398   -0.1026 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4889   -0.7939    0.5501 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.4615    0.1460    0.3535 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0116    1.4034   -0.0296 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9786    1.4264   -0.9435 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1399    0.3193   -0.9885 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  1  0
  8  9  2  0
  9 10  1  0
 10 11  2  0
 11 12  1  0
 12 13  2  0
  6  1  1  0
 13  8  1  0
M  END
"""
            m = Chem.MolFromMolBlock(molblock)
            pieces = BreakBRICSBonds(m)

            frags = Chem.GetMolFrags(pieces, asMols=True)
            self.assertEqual(len(frags), 3)
            self.assertEqual(frags[0].GetNumAtoms(), 7)
            self.assertEqual(frags[1].GetNumAtoms(), 3)
            self.assertEqual(frags[2].GetNumAtoms(), 7)

            c1 = m.GetConformer()
            c2 = frags[0].GetConformer()
            for i in range(6):
                p1 = c1.GetAtomPosition(i)
                p2 = c2.GetAtomPosition(i)
                self.assertEqual((p1 - p2).Length(), 0.0)
            p1 = c1.GetAtomPosition(6)
            p2 = c2.GetAtomPosition(6)
            self.assertEqual((p1 - p2).Length(), 0.0)

            c2 = frags[2].GetConformer()
            for i in range(6):
                p1 = c1.GetAtomPosition(i + 7)
                p2 = c2.GetAtomPosition(i)
                self.assertEqual((p1 - p2).Length(), 0.0)
            p1 = c1.GetAtomPosition(6)
            p2 = c2.GetAtomPosition(6)
            self.assertEqual((p1 - p2).Length(), 0.0)

            c2 = frags[1].GetConformer()
            for i in range(1):
                p1 = c1.GetAtomPosition(i + 6)
                p2 = c2.GetAtomPosition(i)
                self.assertEqual((p1 - p2).Length(), 0.0)
            p1 = c1.GetAtomPosition(5)
            p2 = c2.GetAtomPosition(1)
            self.assertEqual((p1 - p2).Length(), 0.0)
            p1 = c1.GetAtomPosition(6)
            p2 = c2.GetAtomPosition(0)
            self.assertEqual((p1 - p2).Length(), 0.0)

            # make sure multiple conformations (include 2D) also work:
            molblock = """
     RDKit          2D

 13 14  0  0  0  0  0  0  0  0999 V2000
   -1.2990   -0.8654    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5981   -1.6154    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8971   -0.8654    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8971    0.6346    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5981    1.3846    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    0.6346    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000    1.3846    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.6346    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990   -0.8654    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -1.6154    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8971   -0.8654    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8971    0.6346    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981    1.3846    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  1  0
  8  9  2  0
  9 10  1  0
 10 11  2  0
 11 12  1  0
 12 13  2  0
  6  1  1  0
 13  8  1  0
M  END
"""
            m2 = Chem.MolFromMolBlock(molblock)
            m.AddConformer(m2.GetConformer(), assignId=True)
            self.assertEqual(m.GetNumConformers(), 2)

            pieces = BreakBRICSBonds(m)
            frags = Chem.GetMolFrags(pieces, asMols=True)
            self.assertEqual(len(frags), 3)
            self.assertEqual(frags[0].GetNumAtoms(), 7)
            self.assertEqual(frags[1].GetNumAtoms(), 3)
            self.assertEqual(frags[2].GetNumAtoms(), 7)
            self.assertEqual(frags[0].GetNumConformers(), 2)
            self.assertEqual(frags[1].GetNumConformers(), 2)
            self.assertEqual(frags[2].GetNumConformers(), 2)

            c1 = m.GetConformer(0)
            c2 = frags[0].GetConformer(0)
            for i in range(6):
                p1 = c1.GetAtomPosition(i)
                p2 = c2.GetAtomPosition(i)
                self.assertEqual((p1 - p2).Length(), 0.0)
            p1 = c1.GetAtomPosition(6)
            p2 = c2.GetAtomPosition(6)
            self.assertEqual((p1 - p2).Length(), 0.0)

            c2 = frags[2].GetConformer(0)
            for i in range(6):
                p1 = c1.GetAtomPosition(i + 7)
                p2 = c2.GetAtomPosition(i)
                self.assertEqual((p1 - p2).Length(), 0.0)
            p1 = c1.GetAtomPosition(6)
            p2 = c2.GetAtomPosition(6)
            self.assertEqual((p1 - p2).Length(), 0.0)

            c2 = frags[1].GetConformer(0)
            for i in range(1):
                p1 = c1.GetAtomPosition(i + 6)
                p2 = c2.GetAtomPosition(i)
                self.assertEqual((p1 - p2).Length(), 0.0)
            p1 = c1.GetAtomPosition(5)
            p2 = c2.GetAtomPosition(1)
            self.assertEqual((p1 - p2).Length(), 0.0)
            p1 = c1.GetAtomPosition(6)
            p2 = c2.GetAtomPosition(0)
            self.assertEqual((p1 - p2).Length(), 0.0)

            c1 = m.GetConformer(1)
            c2 = frags[0].GetConformer(1)
            for i in range(6):
                p1 = c1.GetAtomPosition(i)
                p2 = c2.GetAtomPosition(i)
                self.assertEqual((p1 - p2).Length(), 0.0)
            p1 = c1.GetAtomPosition(6)
            p2 = c2.GetAtomPosition(6)
            self.assertEqual((p1 - p2).Length(), 0.0)

            c2 = frags[2].GetConformer(1)
            for i in range(6):
                p1 = c1.GetAtomPosition(i + 7)
                p2 = c2.GetAtomPosition(i)
                self.assertEqual((p1 - p2).Length(), 0.0)
            p1 = c1.GetAtomPosition(6)
            p2 = c2.GetAtomPosition(6)
            self.assertEqual((p1 - p2).Length(), 0.0)

            c2 = frags[1].GetConformer(1)
            for i in range(1):
                p1 = c1.GetAtomPosition(i + 6)
                p2 = c2.GetAtomPosition(i)
                self.assertEqual((p1 - p2).Length(), 0.0)
            p1 = c1.GetAtomPosition(5)
            p2 = c2.GetAtomPosition(1)
            self.assertEqual((p1 - p2).Length(), 0.0)
            p1 = c1.GetAtomPosition(6)
            p2 = c2.GetAtomPosition(0)
            self.assertEqual((p1 - p2).Length(), 0.0)

        def test12(self):
            m = Chem.MolFromSmiles('CCS(=O)(=O)NCC')
            res = list(FindBRICSBonds(m))
            self.assertEqual(len(res), 2, res)
            atIds = [x[0] for x in res]
            atIds.sort()
            self.assertEqual(atIds, [(5, 2), (6, 5)])

        def testGithub1734(self):
            m = Chem.MolFromSmiles('c1ccccc1[C@H](C)NC')
            res = BRICSDecompose(m)
            self.assertEqual(len(res), 3)
            self.assertTrue('[4*][C@H]([8*])C' in res)
            res = BreakBRICSBonds(m)
            self.assertEqual(Chem.MolToSmiles(res, isomericSmiles=True),
                             '[16*]c1ccccc1.[4*][C@H]([8*])C.[5*]NC')

    failed, tried = _test()
    if failed:
        sys.exit(failed)

    unittest.main()
