# $Id$
#
#  Copyright (C) 2002-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for the EState indices

validation values are from the paper (JCICS _31_ 76-81 (1991))

"""

import unittest
from io import StringIO

import numpy as np
from rdkit import Chem
from rdkit.Chem import EState


class TestCase(unittest.TestCase):

    def _compareEstates(self, val1, val2, msg, tol=1e-2):
        maxV = max(abs(val1 - val2))
        self.assertLess(maxV, tol, msg)

    def _validate(self, vals, places=2, tol=1e-2, debug=False):
        for smi, ans in vals:
            ans = np.array(ans)
            mol = Chem.MolFromSmiles(smi)
            inds = EState.EStateIndices(mol)
            if debug:  # pragma: nocover
                print(inds)
            self._compareEstates(ans, inds, 'bad EStates for smiles: {0}'.format(smi), tol=tol)

            self.assertLess(abs(EState.MaxEStateIndex(mol) - max(ans)), tol)
            self.assertLess(abs(EState.MinEStateIndex(mol) - min(ans)), tol)
            self.assertLess(abs(EState.MaxAbsEStateIndex(mol) - max(abs(ans))), tol)
            self.assertLess(abs(EState.MinAbsEStateIndex(mol) - min(abs(ans))), tol)

    def test_simpleMolecules(self):
        data = [
          ('CCCC', [2.18, 1.32, 1.32, 2.18]),
          ('CCCCC', [2.21, 1.34, 1.39, 1.34, 2.21]),
          ('CCCCCCC', [2.24, 1.36, 1.42, 1.44, 1.42, 1.36, 2.24]),
          ('CCCCCCCCCC', [2.27, 1.37, 1.44, 1.46, 1.47, 1.47, 1.46, 1.44, 1.37, 2.27]),
        ]
        self._validate(data)

    def test_isomers(self):
        data = [
          ('CCCCCC', [2.23, 1.36, 1.41, 1.41, 1.36, 2.23]),
          ('CCC(C)CC', [2.23, 1.33, 0.94, 2.28, 1.33, 2.23]),
          ('CC(C)CCC', [2.25, 0.90, 2.25, 1.38, 1.33, 2.22]),
          ('CC(C)(C)CC', [2.24, 0.54, 2.24, 2.24, 1.27, 2.20]),
        ]
        self._validate(data)

    def test_heteroatoms1(self):
        data = [
          ('CCCCOCCCC', [2.18, 1.24, 1.21, 0.95, 5.31, 0.95, 1.21, 1.24, 2.18]),
          ('CCC(C)OC(C)CC', [2.15, 1.12, 0.43, 2.12, 5.54, 0.43, 2.12, 1.12, 2.15]),
          ('CC(C)(C)OC(C)(C)C', [2.07, -0.02, 2.07, 2.07, 5.63, -0.02, 2.07, 2.07, 2.07]),
          ('CC(C)CC', [2.22, 0.88, 2.22, 1.31, 2.20]),
          ('CC(C)CN', [2.10, 0.66, 2.10, 0.81, 5.17]),
          ('CC(C)CO', [1.97, 0.44, 1.97, 0.31, 8.14]),
          ('CC(C)CF', [1.85, 0.22, 1.85, -0.19, 11.11]),
          ('CC(C)CCl', [2.09, 0.65, 2.09, 0.78, 5.34]),
          ('CC(C)CBr', [2.17, 0.80, 2.17, 1.11, 3.31]),
          ('CC(C)CI', [2.21, 0.87, 2.21, 1.28, 2.38]),
        ]
        self._validate(data, debug=False)

    def test_heteroatoms2(self):
        data = [
          ('CC(N)C(=O)O', [1.42, -0.73, 4.84, -0.96, 9.57, 7.86]),
          ('CCOCC', [1.99, 0.84, 4.83, 0.84, 1.99]),
          # NOTE: this doesn't match the values in the paper
          ('CCSCC', [2.17, 1.26, 1.96, 1.26, 2.17]),
          ('CC(=O)OC', [1.36, -0.24, 9.59, 4.11, 1.35]),
          ('CC(=S)OC', [1.73, 0.59, 4.47, 4.48, 1.56]),
        ]
        self._validate(data, debug=False)

    def test_aromatics(self):
        # aromatics with heteroatoms
        data = [
          ('Fc1ccc(C)cc1', [12.09, -0.17, 1.45, 1.75, 1.09, 1.93, 1.75, 1.45]),
          ('Clc1ccc(C)cc1', [5.61, 0.80, 1.89, 1.99, 1.24, 2.04, 1.99, 1.89]),
          ('Brc1ccc(C)cc1', [3.35, 1.14, 2.04, 2.07, 1.30, 2.08, 2.07, 2.04]),
          ('Ic1ccc(C)cc1', [2.30, 1.30, 2.10, 2.11, 1.32, 2.09, 2.11, 2.10]),
        ]
        self._validate(data, debug=False)

    def test_GetPrincipleQuantumNumber(self):
        for principalQN, (nmin, nmax) in enumerate(
          [(1, 2), (3, 10), (11, 18), (19, 36), (37, 54), (55, 86), (87, 120)], 1):
            for n in range(nmin, nmax + 1):
                self.assertEqual(EState.GetPrincipleQuantumNumber(n), principalQN)

    def test_cacheEstate(self):
        mol = Chem.MolFromSmiles('CCCC')
        expected = [2.18, 1.32, 1.32, 2.18]

        # The mol object has no information about E-states
        self.assertFalse(hasattr(mol, '_eStateIndices'))
        inds = EState.EStateIndices(mol)
        self._compareEstates(inds, expected, 'cacheTest')

        # We now have E-states stored with the molecule
        self.assertTrue(hasattr(mol, '_eStateIndices'))

        # Let's make sure that we skip the calculation next time if force is False
        mol._eStateIndices = 'cached'
        self.assertTrue(hasattr(mol, '_eStateIndices'))

        inds = EState.EStateIndices(mol, force=False)
        self.assertEqual(inds, 'cached')

        # But with force (default) we calculate again
        inds = EState.EStateIndices(mol)
        self._compareEstates(inds, expected, 'cacheTest')
        self._compareEstates(mol._eStateIndices, expected, 'cacheTest')

    def test_exampleCode(self):
        # We make sure that the example code runs
        from rdkit.TestRunner import redirect_stdout
        f = StringIO()
        with redirect_stdout(f):
            EState.EState._exampleCode()
        s = f.getvalue()
        self.assertIn('CC(N)C(=O)O', s)


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
