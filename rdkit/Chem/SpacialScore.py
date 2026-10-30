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
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.ChemUtils.DescriptorUtilities import setDescriptorVersion


@setDescriptorVersion(version="1.0.0")
def SPS(mol, normalize=True):
  """Calculates the SpacialScore descriptor. By default, the score is normalized by the number of heavy atoms (nSPS) resulting in a float value,
    otherwise (normalize=False) the absolute score is returned as an integer.
    """
  if mol is None:
    raise ValueError("No valid molecule object found.")
  res = rdMolDescriptors.CalcSPS(mol, normalize=bool(normalize))
  return res if normalize else int(res)

