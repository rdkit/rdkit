"""
unit tests for the MolStandardize module
tests include:
reorder_tautomers
"""

import unittest

from rdkit import Chem
from rdkit.Chem import MolStandardize
from rdkit.Chem.MolStandardize import rdMolStandardize


class TestCase(unittest.TestCase):

  def testBasic(self):
    m = Chem.MolFromSmiles('Oc1c(cccc3)c3nc2ccncc12')
    enumerator = rdMolStandardize.TautomerEnumerator()
    canon = enumerator.Canonicalize(m)
    reord = MolStandardize.ReorderTautomers(m)[0]
    canonSmile = Chem.MolToSmiles(canon)
    reordSmile = Chem.MolToSmiles(reord)
    self.assertEqual(canonSmile, reordSmile)

  def testLength(self):
    m = Chem.MolFromSmiles('Oc1c(cccc3)c3nc2ccncc12')
    enumerator = rdMolStandardize.TautomerEnumerator()
    tauts = enumerator.Enumerate(m)
    reordtauts = MolStandardize.ReorderTautomers(m)
    self.assertEqual(len(reordtauts), len(tauts))

  def testNoUnassignedStereoAfterCanonicalize(self):
    # tautomer canonicalization must not produce an unassigned stereocenter ('?')
    # due to atom/bond reordering + stereochem reassignment.
    # Enhanced StereoGroup over the stereocenters to indicate the `syn` nature of mol
    molblock = r"""
  ChemDraw01232415492D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 12 13 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O 1.690391 0.714714 0.000000 0
M  V30 2 C 1.278465 0.000859 0.000000 0
M  V30 3 O 1.690391 -0.712995 0.000000 0
M  V30 4 C 0.454037 0.000286 0.000000 0
M  V30 5 C -0.260390 0.412786 0.000000 0
M  V30 6 H 0.153828 1.125495 0.000000 0
M  V30 7 C -0.975391 0.831588 0.000000 0
M  V30 8 C -1.690391 0.412786 0.000000 0
M  V30 9 C -1.690391 -0.415651 0.000000 0
M  V30 10 C -0.975391 -0.825860 0.000000 0
M  V30 11 C -0.260390 -0.412214 0.000000 0
M  V30 12 H 0.152110 -1.125495 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 4 2 CFG=3
M  V30 4 1 4 11
M  V30 5 1 4 5
M  V30 6 1 5 6 CFG=3
M  V30 7 1 5 11
M  V30 8 1 5 7
M  V30 9 1 7 8
M  V30 10 1 8 9
M  V30 11 1 9 10
M  V30 12 1 10 11
M  V30 13 1 11 12 CFG=3
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    base = Chem.MolFromMolBlock(molblock, sanitize=True, removeHs=False)
    enumerator = rdMolStandardize.TautomerEnumerator()
    enumerator.SetRemoveSp3Stereo(False)
    canon = enumerator.Canonicalize(base)
    centers = Chem.FindMolChiralCenters(
      canon,
      includeUnassigned=True,
      useLegacyImplementation=False,
    )
    unassigned = [c for c in centers if c[1] == '?']
    self.assertFalse(
      unassigned,
      msg=f"produced unassigned stereocenters: {unassigned} (centers={centers})",
    )

