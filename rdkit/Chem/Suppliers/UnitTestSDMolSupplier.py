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
"""
unit testing code for the SD file handling stuff
"""
import os
import tempfile
import unittest

from rdkit import Chem
from rdkit import RDConfig


class TestCase(unittest.TestCase):

    def setUp(self):
        self.fName = os.path.join(RDConfig.RDDataDir, 'NCI', 'first_200.props.sdf')
        with open(self.fName, 'r') as inf:
            inD = inf.read()
        self.nMolecules = inD.count('$$$$')

    def assertMolecule(self, mol, i):
        """ Assert that we have a valid molecule """
        self.assertIsNotNone(mol, 'read %d failed' % i)
        self.assertGreater(mol.GetNumAtoms(), 0, 'no atoms in mol %d' % i)

    def test_SDMolSupplier(self):
        # tests reads using a file name (file contains 200 molecules)
        supp = Chem.SDMolSupplier(self.fName)

        # Can use as an iterator
        for i in range(10):
            mol = next(supp)
            self.assertMolecule(mol, i)

        # Can access directly
        i = 100
        mol = supp[i - 1]
        self.assertMolecule(mol, i)

        # We can access the number of molecules
        self.assertEqual(len(supp), self.nMolecules, 'bad supplier length')

        # We know the number and can still access directly
        i = 12
        mol = supp[i - 1]
        self.assertMolecule(mol, i)

        # Get an exception if we access an invalid number
        with self.assertRaises(IndexError):
            _ = supp[self.nMolecules]  # out of bound read must fail

        # and we can access with negative numbers
        mol1 = supp[len(supp) - 1]
        mol2 = supp[-1]
        self.assertEqual(Chem.MolToSmiles(mol1), Chem.MolToSmiles(mol2))

    def test_SDWriter(self):
        # tests writes using a file name
        supp = Chem.SDMolSupplier(self.fName)
        outName = tempfile.NamedTemporaryFile(suffix='.sdf', delete=False).name
        writer = Chem.SDWriter(outName)
        m1 = next(supp)
        writer.SetProps(m1.GetPropNames())
        for m in supp:
            writer.write(m)
        writer.flush()
        writer.close()

        # The writer does not have an explicit "close()" so need to
        # let the garbage collector kick in to close the file.
        writer = None
        with open(outName, 'r') as inf:
            outD = inf.read()
        # The file should be closed, but if it isn't, and this
        # is Windows, then the unlink() can fail. Wait and try again.
        try:
            os.unlink(outName)
        except Exception:
            import time
            time.sleep(1)
            try:
                os.unlink(outName)
            except Exception:
                pass
        self.assertEqual(self.nMolecules, outD.count('$$$$'), 'bad nMols in output')

#   def _testStreamRoundtrip(self):
#     inD = open(self.fName).read()
#     supp = Chem.SDMolSupplier(self.fName)
#     outName = tempfile.mktemp('.sdf')
#     writer = Chem.SDWriter(outName)
#     _ = next(supp)
#     for m in supp:
#       writer.write(m)
#     writer.flush()
#     writer = None
#     outD = open(outName, 'r').read()
#     try:
#       os.unlink(outName)
#     except Exception:
#       import time
#       time.sleep(1)
#       try:
#         os.unlink(outName)
#       except Exception:
#         pass
#     assert inD.count('$$$$') == outD.count('$$$$'), 'bad nMols in output'
#     io = StringIO(outD)
#     supp = Chem.SDMolSupplier(stream=io)
#     outD2 = supp.Dump()
#     assert outD2.count('$$$$') == len(supp), 'bad nMols in output'
#     assert outD2.count('$$$$') == outD.count('$$$$'), 'bad nMols in output'
#     assert outD2 == outD, 'bad outd'

#   def _testLazyDataRoundtrip(self):
#     inD = open(self.fName).read()
#     supp = Chem.SDMolSupplier(self.fName)
#     outName = tempfile.mktemp('.sdf')
#     writer = Chem.SDWriter(outName)
#     _ = next(supp)
#     for m in supp:
#       writer.write(m)
#     writer.flush()
#     writer = None
#     outD = open(outName, 'r').read()
#     try:
#       os.unlink(outName)
#     except Exception:
#       import time
#       time.sleep(1)
#       try:
#         os.unlink(outName)
#       except Exception:
#         pass
#     assert inD.count('$$$$') == outD.count('$$$$'), 'bad nMols in output'
#     supp = Chem.SDMolSupplier.LazySDMolSupplier(inD=outD)
#     outD2 = supp.Dump()
#     assert outD2.count('$$$$') == len(supp), 'bad nMols in output'
#     assert outD2.count('$$$$') == outD.count('$$$$'), 'bad nMols in output'
#     assert outD2 == outD, 'bad outd'

#   def _testLazyIter(self):
#     " tests lazy reads using the iterator interface "
#     supp = Chem.SDMolSupplier.LazySDMolSupplier(fileN=self.fName)
#
#     nDone = 0
#     for mol in supp:
#       assert mol, 'read %d failed' % nDone
#       assert mol.GetNumAtoms(), 'no atoms in mol %d' % nDone
#       nDone += 1
#     assert nDone == 200, 'bad number of molecules: %d' % (nDone)
#
#     l = len(supp)
#     assert l == 200, 'bad supplier length: %d' % (l)
#
#     i = 12
#     m = supp[i - 1]
#     assert m, 'back index %d failed' % i
#     assert m.GetNumAtoms(), 'no atoms in mol %d' % i
#
#     try:
#       m = supp[201]
#     except IndexError:
#       fail = 1
#     else:
#       fail = 0
#     assert fail, 'out of bound read did not fail'


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
