# Copyright (c) 2013, GlaxoSmithKline Research & Development Ltd.
# All rights reserved.
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
#     * Neither the name of GlaxoSmithKline Research & Development Ltd.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written
#       permission.
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
# Created by Jameed Hussain, May 2013
"""
Fragmentation algorithm
-----------------------

identify acyclic bonds
enumerate all single cuts
make sure you chop off more that 1 atom
keeps bits which are >60% query mol
enumerate all double cuts
keeps bits with 1 attachment point (i.e throw middle bit away)
need to be >60% query mol

identify exocyclic bonds
enumerate all single "ring" cuts
Check if it results in more that one component
keep correct bit if >40% query mol

enumerate successful "rings" cuts with an acyclic cut
Check if it results in more that one component
keep correct if >60% query mol

"""
import sys
from itertools import combinations

from rdkit import Chem, DataStructs
from rdkit.Chem import rdqueries

# our default rdkit fingerprinter parameters:
rdkitFpParams = {'maxPath': 5, 'fpSize': 1024, 'nBitsPerHash': 2}

# Considered fragment types
FTYPE_ACYCLIC = 'acyclic'
FTYPE_CYCLIC = 'cyclic'
FTYPE_CYCLIC_ACYCLIC = 'cyclic_and_acyclic'

# Global SMARTS used by the program

# acyclic bond smarts
ACYC_SMARTS = Chem.MolFromSmarts("*!@!=!#*")
# exocyclic/fused exocyclic bond smarts
CYC_SMARTS = Chem.MolFromSmarts("[R1,R2]@[r;!R1]")

# smarts used to find appropriate fragment for
# would use SMARTS: [$([#0][r].[r][#0]),$([#0][r][#0])]
# but RDkit doesn't support component SMARTS in recursive one - $([#0][r].[r][#0])
# hence split into two
cSma1 = Chem.MolFromSmarts("[#0][r].[r][#0]")
cSma2 = Chem.MolFromSmarts("[#0][r][#0]")
dummyAtomQuery = rdqueries.AtomNumEqualsQueryAtom(0)


def delete_bonds(mol, bonds, ftype, hac):
  """ Fragment molecule on bonds and reduce to fraggle fragmentation SMILES.
  If none exists, returns None """

  # Replace the given bonds with attachment points (B1-B2 -> B1-*.*-B2)
  bondIdx = [mol.GetBondBetweenAtoms(*bond).GetIdx() for bond in bonds]
  modifiedMol = Chem.FragmentOnBonds(mol, bondIdx, dummyLabels=[(0, 0)] * len(bondIdx))

  # should be able to get away without sanitising mol as the valencies should be okay
  # do not do a full sanitization, but do find rings and calculate valences:
  Chem.SanitizeMol(modifiedMol,
                   Chem.SanitizeFlags.SANITIZE_PROPERTIES | Chem.SanitizeFlags.SANITIZE_SYMMRINGS)

  fragments = Chem.GetMolFrags(modifiedMol, asMols=True, sanitizeFrags=False)
  return select_fragments(fragments, ftype, hac)


def select_fragments(fragments, ftype, hac):
  if ftype == FTYPE_ACYCLIC:
    result = []
    result_hcount = 0
    for fMol in fragments:
      nAttachments = len(fMol.GetAtomsMatchingQuery(dummyAtomQuery))
      # check if terminal fragment
      if nAttachments == 1:
        fhac = fMol.GetNumAtoms()

        # if the fragment is 2 atoms (or less - includes attachment) it is too small
        # to be interesting. This check has the additional benefit
        # of pulling out the relevant single cuts as it discards
        # fragments where we only chop off a small part of the input cmpd
        if fhac > 3:
          result.append(Chem.MolToSmiles(fMol))
          result_hcount += fhac

    # needs to be greater than 60% of parent mol
    if result and (result_hcount > 0.6 * hac):
      return '.'.join(result)
    return None

  if ftype == FTYPE_CYCLIC:
    # make sure it is 2 components
    if len(fragments) != 2:
      return None
    result = None
    for fMol in fragments:
      f = Chem.MolToSmiles(fMol)
      # check if a valid cut
      # needs to be greater 3 heavy atoms and greater than 40% of parent mol
      if isValidRingCut(fMol):
        result_hcount = fMol.GetNumAtoms()
        if result_hcount > 3 and result_hcount > 0.4 * hac:
          result = f
    return result

  if ftype == FTYPE_CYCLIC_ACYCLIC:
    # need to find the fragments which are valid which means they must be:
    #  Terminal (one attachment point) or valid ring cut
    result = []
    result_hcount = 0
    for fMol in fragments:
      nAttachments = len(fMol.GetAtomsMatchingQuery(dummyAtomQuery))
      # We need to have a fragment that has 1 or 2 attachment points and that has more than 3 atoms
      if nAttachments >= 3:
        continue
      fhac = fMol.GetNumAtoms()
      if fhac <= 3:
        continue

      if nAttachments == 2:
        # check if a valid cut
        if isValidRingCut(fMol):
          result.append(Chem.MolToSmiles(fMol))
          result_hcount += fhac
      elif nAttachments == 1:
        result.append(Chem.MolToSmiles(fMol))
        result_hcount += fhac

    # appropriate fragmentation must have 2 components and needs to be greater than 60% of
    # parent mol
    if len(result) == 2 and result_hcount > 0.6 * hac:
      return '.'.join(result)
    return None
  raise NotImplementedError(f'Invalid fragmentation type {type}')


def isValidRingCut(mol):
  """ to check is a fragment is a valid ring cut, it needs to match the
  SMARTS: [$([#0][r].[r][#0]),$([#0][r][#0])] """
  # At this point, the molecule requires the identification of rings, so we need to sanitize
  Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_SYMMRINGS)
  return mol.HasSubstructMatch(cSma1) or mol.HasSubstructMatch(cSma2)


def generate_fraggle_fragmentation(mol, verbose=False):
  """ Create all possible fragmentations for molecule
    >>> q = Chem.MolFromSmiles('COc1cc(CN2CCC(NC(=O)c3cncc(C)c3)CC2)c(OC)c2ccccc12')
    >>> fragments = generate_fraggle_fragmentation(q)
    >>> fragments = sorted(['.'.join(sorted(s.split('.'))) for s in fragments])
    >>> fragments
     ['*C(=O)NC1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1',
      '*C(=O)c1cncc(C)c1.*C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1',
      '*C(=O)c1cncc(C)c1.*Cc1cc(OC)c2ccccc2c1OC',
      '*C(=O)c1cncc(C)c1.*c1cc(OC)c2ccccc2c1OC',
      '*C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1',
      '*C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1.*c1cncc(C)c1',
      '*Cc1cc(OC)c2ccccc2c1OC.*NC(=O)c1cncc(C)c1',
      '*Cc1cc(OC)c2ccccc2c1OC.*c1cncc(C)c1',
      '*N1CCC(NC(=O)c2cncc(C)c2)CC1.*c1cc(OC)c2ccccc2c1OC',
      '*NC(=O)c1cncc(C)c1.*c1cc(OC)c2ccccc2c1OC',
      '*NC1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1',
      '*NC1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1.*c1cncc(C)c1',
      '*c1c(CN2CCC(NC(=O)c3cncc(C)c3)CC2)cc(OC)c2ccccc12',
      '*c1c(OC)cc(CN2CCC(NC(=O)c3cncc(C)c3)CC2)c(OC)c1*',
      '*c1cc(CN2CCC(NC(=O)c3cncc(C)c3)CC2)c(OC)c2ccccc12',
      '*c1cc(OC)c2ccccc2c1OC.*c1cncc(C)c1']
  """
  # query mol heavy atom count
  hac = mol.GetNumAtoms()

  # find the relevant bonds to break
  acyclic_matching_atoms = mol.GetSubstructMatches(ACYC_SMARTS)
  cyclic_matching_atoms = mol.GetSubstructMatches(CYC_SMARTS)
  if verbose:
    print("Matching Atoms:")
    print("acyclic matching atoms: ", acyclic_matching_atoms)
    print("cyclic matching atoms: ", cyclic_matching_atoms)

  # different cuts can give the same fragments
  # to use out_fragments to remove them
  out_fragments = set()

  ######################
  # Single acyclic Cuts
  ######################
  # loop to generate every single and double cut in the molecule
  # single cuts are not required as relevant single cut fragments can be found
  # from the double cuts. For explanation see check_fragments method
  for bond1, bond2 in combinations(acyclic_matching_atoms, 2):
    fragment = delete_bonds(mol, [bond1, bond2], FTYPE_ACYCLIC, hac)
    if fragment is not None:
      out_fragments.add(fragment)

  ##################################
  # Fused/Spiro exocyclic bond Cuts
  ##################################
  for bond1, bond2 in combinations(cyclic_matching_atoms, 2):
    fragment = delete_bonds(mol, [bond1, bond2], FTYPE_CYCLIC, hac)
    if fragment is None:
      continue
    out_fragments.add(fragment)
    # now do an acyclic cut with the successful cyclic cut
    for abond in acyclic_matching_atoms:
      fragment = delete_bonds(mol, [bond1, bond2, abond], FTYPE_CYCLIC_ACYCLIC, hac)
      if fragment is not None:
        out_fragments.add(fragment)

  return sorted(out_fragments)


def atomContrib(subs, mol, tverskyThresh=0.8):
  """ atomContrib algorithm
  generate fp of query_substructs (qfp)

  loop through atoms of smiles
    For each atom
    Generate partial fp of the atom (pfp)
    Find Tversky sim of pfp in qfp
    If Tversky < 0.8, mark atom in smiles

  Loop through marked atoms
    If marked atom in ring - turn all atoms in that ring to * (aromatic) or Sc (aliphatic)
    For each marked atom
      If aromatic turn to a *
      If aliphatic turn to a Sc

  Return modified smiles
  """

  def partialSimilarity(atomID):
    """ Determine similarity for the atoms set by atomID """
    # create empty fp
    modifiedFP = DataStructs.ExplicitBitVect(1024)
    modifiedFP.SetBitsFromList(aBits[atomID])
    return DataStructs.TverskySimilarity(subsFp, modifiedFP, 0, 1)

  # generate mol object & fp for input mol (we are interested in the bits each atom sets)
  pMol = Chem.Mol(mol)
  aBits = []
  _ = Chem.RDKFingerprint(pMol, atomBits=aBits, **rdkitFpParams)

  # generate fp of query_substructs
  qsMol = Chem.MolFromSmiles(subs)
  subsFp = Chem.RDKFingerprint(qsMol, **rdkitFpParams)

  # loop through atoms of smiles get atoms that have a high similarity with substructure
  marked = set()
  for atom in pMol.GetAtoms():
    atomIdx = atom.GetIdx()
    if partialSimilarity(atomIdx) < tverskyThresh:
      marked.add(atomIdx)

  # get rings to change

  # If a marked atom is within a ring, mark the whole ring
  markRingAtoms = set()
  for ring in pMol.GetRingInfo().AtomRings():
    if any(ringAtom in marked for ringAtom in ring):
      markRingAtoms.update(ring)
  marked.update(markRingAtoms)

  if marked:
    # now mutate the marked atoms
    for idx in marked:
      if pMol.GetAtomWithIdx(idx).GetIsAromatic():
        pMol.GetAtomWithIdx(idx).SetAtomicNum(0)
        pMol.GetAtomWithIdx(idx).SetNoImplicit(True)
      else:
        # gives best sim
        pMol.GetAtomWithIdx(idx).SetAtomicNum(21)
        # works better but when replace S it fails due to valency
        # pMol.GetAtomWithIdx(idx).SetAtomicNum(6)

    try:
      Chem.SanitizeMol(
        pMol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE ^ Chem.SANITIZE_SETAROMATICITY)
    except Exception:
      sys.stderr.write(f"Can't parse smiles: {Chem.MolToSmiles(pMol)}\n")
      pMol = Chem.Mol(mol)
  return pMol


modified_query_fps = {}


def compute_fraggle_similarity_for_subs(inMol, qMol, qSmi, qSubs, tverskyThresh=0.8):
  qFP = Chem.RDKFingerprint(qMol, **rdkitFpParams)
  iFP = Chem.RDKFingerprint(inMol, **rdkitFpParams)

  rdkit_sim = DataStructs.TanimotoSimilarity(qFP, iFP)

  qm_key = f"{qSubs}_{qSmi}"
  if qm_key in modified_query_fps:
    qmMolFp = modified_query_fps[qm_key]
  else:
    qmMol = atomContrib(qSubs, qMol, tverskyThresh)
    qmMolFp = Chem.RDKFingerprint(qmMol, **rdkitFpParams)
    modified_query_fps[qm_key] = qmMolFp

  rmMol = atomContrib(qSubs, inMol, tverskyThresh)

  # wrap in a try, catch
  try:
    rmMolFp = Chem.RDKFingerprint(rmMol, **rdkitFpParams)
    fraggle_sim = max(DataStructs.FingerprintSimilarity(qmMolFp, rmMolFp), rdkit_sim)
  except Exception:
    sys.stderr.write(f"Can't generate fp for: {Chem.MolToSmiles(rmMol)}\n")
    fraggle_sim = 0.0

  return rdkit_sim, fraggle_sim


def GetFraggleSimilarity(queryMol, refMol, tverskyThresh=0.8):
  """ return the Fraggle similarity between two molecules

    >>> q = Chem.MolFromSmiles('COc1cc(CN2CCC(NC(=O)c3cncc(C)c3)CC2)c(OC)c2ccccc12')
    >>> m = Chem.MolFromSmiles('COc1cc(CN2CCC(NC(=O)c3ccccc3)CC2)c(OC)c2ccccc12')
    >>> sim,match = GetFraggleSimilarity(q,m)
    >>> sim
    0.980...
    >>> match
    '*C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1'

    >>> m = Chem.MolFromSmiles('COc1cc(CN2CCC(Nc3nc4ccccc4s3)CC2)c(OC)c2ccccc12')
    >>> sim,match = GetFraggleSimilarity(q,m)
    >>> sim
    0.794...
    >>> match
    '*C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1'

    >>> q = Chem.MolFromSmiles('COc1ccccc1')
    >>> sim,match = GetFraggleSimilarity(q,m)
    >>> sim
    0.347...
    >>> match
    '*c1ccccc1'

    """
  if hasattr(queryMol, '_fraggleDecomp'):
    frags = queryMol._fraggleDecomp
  else:
    frags = generate_fraggle_fragmentation(queryMol)
    queryMol._fraggleDecomp = frags
  qSmi = Chem.MolToSmiles(queryMol, True)
  result = 0.0
  bestMatch = None
  for frag in frags:
    _, fragsim = compute_fraggle_similarity_for_subs(refMol, queryMol, qSmi, frag, tverskyThresh)
    if fragsim > result:
      result = fragsim
      bestMatch = frag
  return result, bestMatch


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import doctest
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS + doctest.NORMALIZE_WHITESPACE,
                              verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
