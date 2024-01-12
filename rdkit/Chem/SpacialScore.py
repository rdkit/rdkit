#
#  Copyright (c) 2023, RDKit Hackathon 2023, implemented by:
# - Arsenio Cruz
# - José-Manuel Gally
# - Axel Pahl
# - Vincenzo Palmacci
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above
#     copyright notice, this list of conditions and the following
#     disclaimer in the documentation and/or other materials provided
#     with the distribution.
#   * Neither the name of Novartis Institutes for BioMedical Research Inc.
#     nor the names of its contributors may be used to endorse or promote
#     products derived from this software without specific prior written permission.
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
"""
Spacial score (SPS) is an empirical scoring system to express the spacial complexity of a compound
in an uniform manner and on a highly granular scale for ranking and comparison between molecules. [1]
SPS takes into account the fraction of sp3 hybridized carbons and the fraction of stereogenic carbons.

By default, this module generates the normalized spacial score (nSPS), which is a variation of the SPS score
that considers the size of the molecule.
To obtain the nSPS score the SPS score is divided by the total number of heavy atoms in the analyzed molecule.

SPS = sum(h*s*r*n*n)
nSPS = SPS/a

Where:
h = Atom hybridisation term
s = Stereoisomeric term
r = Non-aromatic ring term
n = Number of heavy atom neighbors
a = Total number of heavy atoms in the molecule

The SPS function in this module takes a mol object and returns either the absolute score (normalize=False) or the score normalized by the number of heavy atoms (normalize=True (default)).

The original code implementation can be found at: https://github.com/frog2000/Spacial-Score/blob/main/spacial_score.py

[1] Krzyzanowski, A.; Pahl, A.; Grigalunas, M.; Waldmann, H. Spacial Score─A Comprehensive Topological Indicator for Small-Molecule Complexity. J. Med. Chem. 2023. https://doi.org/10.1021/acs.jmedchem.3c00689.
"""
from rdkit import Chem

# import rdkit.Chem.Descriptors as Desc
from rdkit.Chem import rdmolops
from rdkit.Chem.ChemUtils.DescriptorUtilities import setDescriptorVersion
from collections import defaultdict


@setDescriptorVersion(version="1.0.0")
def SPS(mol, normalize=True):
  """Calculates the SpacialScore descriptor. By default, the score is normalized by the number of heavy atoms (nSPS) resulting in a float value,
    otherwise (normalize=False) the absolute score is returned as an integer.
    """
  return _SpacialScore(mol, normalize=normalize).score


class _SpacialScore:
  """Class intended for calculating spacial score (SPS) and size-normalised SPS (nSPS) for small organic molecules"""

  def __init__(self, mol, normalize=True):
    if mol is None:
      raise ValueError("No valid molecule object found.")
    molCp = Chem.Mol(mol)
    rdmolops.FindPotentialStereoBonds(molCp)
    self.mol = molCp  # mol is supposed to be a valid RDKit Mol object
    self.normalize = normalize  # if true nSPS, otherwise SPS
    self.hyb_score = {}
    self.stereo_score = {}
    self.ring_score = {}
    self.bond_score = {}
    self.chiral_idxs = self.findStereoAtomIdxs()
    self.doublebonds_stereo = self.findDoubleBondsStereo()
    # calculate SPS
    self.score = self.calculateSpacialScore()
    # return nSPS
    if normalize:
      self.score /= self.mol.GetNumHeavyAtoms()

  def findStereoAtomIdxs(self, includeUnassigned=True):
    """Finds indices of atoms that are (pseudo)stereo/chiralcentres, in respect to the attached groups (does not account for double bond isomers)"""
    stereo_centers = Chem.FindMolChiralCenters(
      self.mol,
      includeUnassigned=includeUnassigned,
      includeCIP=False,
      useLegacyImplementation=False,
    )
    stereo_idxs = [atom_idx for atom_idx, _ in stereo_centers]
    return stereo_idxs

  def findDoubleBondsStereo(self):
    """Finds indeces of stereo double bond atoms (E/Z)"""
    db_stereo = {}
    for bond in self.mol.GetBonds():
      if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
        db_stereo[(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())] = bond.GetStereo()
    return db_stereo

  def calculateSpacialScore(self):
    """Calculates the total spacial score for a molecule"""
    score = 0
    for atom in self.mol.GetAtoms():
      atom_idx = atom.GetIdx()
      self.hyb_score[atom_idx] = self._accountForHybridisation(atom)
      self.stereo_score[atom_idx] = self._accountForStereo(atom_idx)
      self.ring_score[atom_idx] = self._accountForRing(atom)
      self.bond_score[atom_idx] = self._accountForNeighbors(atom)
      score += self._calculateScoreForAtom(atom_idx)  # absolute score
    return score

  def _calculateScoreForAtom(self, atom_idx):
    """Calculates the total score for a single atom in a molecule"""
    atom_score = (self.hyb_score[atom_idx] * self.stereo_score[atom_idx] *
                  self.ring_score[atom_idx] * self.bond_score[atom_idx])
    return atom_score

  _hybridisations = defaultdict(lambda: 4)
  _hybridisations.update({
    Chem.HybridizationType.SP: 1,
    Chem.HybridizationType.SP2: 2,
    Chem.HybridizationType.SP3: 3
  })

  def _accountForHybridisation(self, atom):
    """Calculates the hybridisation score for a single atom in a molecule"""
    return self._hybridisations[atom.GetHybridization()]

  def _accountForStereo(self, atom_idx):
    """Calculates the stereo score for a single atom in a molecule"""
    if atom_idx in self.chiral_idxs:
      return 2
    for bond_atom_idxs, stereo in self.doublebonds_stereo.items():
      if stereo != Chem.BondStereo.STEREONONE and atom_idx in bond_atom_idxs:
        return 2
    return 1

  def _accountForRing(self, atom):
    """Calculates the ring score for a single atom in a molecule"""
    if atom.GetIsAromatic():  # aromatic rings are not promoted
      return 1
    if atom.IsInRing():
      return 2
    return 1

  def _accountForNeighbors(self, atom):
    """Calculates the neighbour score for a single atom in a molecule
        The second power allows to account for branching in the molecular structure"""
    return atom.GetDegree()**2
