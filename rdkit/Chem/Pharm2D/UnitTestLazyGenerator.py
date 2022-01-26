#
#  Copyright (C) 2003-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for the lazy signature generator

"""
import unittest

from rdkit import Chem
from rdkit.Chem.Pharm2D import SigFactory

try:
  from rdkit.Chem.Pharm2D import LazyGenerator
except NotImplementedError:
  LazyGenerator = None


class TestCase(unittest.TestCase):  # pragma: nocover

  def getFactory(self):
    factory = SigFactory.SigFactory()
    factory.SetPatternsFromSmarts(['O', 'N'])
    factory.SetBins([(0, 2), (2, 5), (5, 8)])
    factory.SetMinCount(2)
    factory.SetMaxCount(3)
    return factory

  def test_NotImplemented(self):
    self.assertIsNone(LazyGenerator, 'Review LazyGenerator unit tests')

  @unittest.skipIf(LazyGenerator is None, 'LazyGenerator implementation incomplete')
  def test1_simple(self):
    mol = Chem.MolFromSmiles('OCC(=O)CCCN')
    factory = self.getFactory()
    sig = factory.GetSignature()
    assert sig.GetSize() == 105, f'bad signature size: {sig.GetSize()}'
    sig.SetIncludeBondOrder(0)
    gen = LazyGenerator.Generator(sig, mol)
    assert len(gen) == sig.GetSize(), f'length mismatch {len(gen)}!={sig.GetSize()}'

    tgt = (1, 5, 48)
    for bit in tgt:
      assert gen[bit], f'bit {bit} not properly set'
      assert gen.GetBit(bit), f'bit {bit} not properly set'
      assert not gen[bit + 50], f'bit {bit + 100} improperly set'

    sig = factory.GetSignature()
    assert sig.GetSize() == 105, f'bad signature size: {sig.GetSize()}'
    sig.SetIncludeBondOrder(1)
    gen = LazyGenerator.Generator(sig, mol)
    assert len(gen) == sig.GetSize(), f'length mismatch {len(gen)}!={sig.GetSize()}'

    tgt = (1, 4, 5, 45)
    for bit in tgt:
      assert gen[bit], f'bit {bit} not properly set'
      assert gen.GetBit(bit), f'bit {bit} not properly set'
      assert not gen[bit + 50], f'bit {bit + 100} improperly set'

    try:
      gen[sig.GetSize() + 1]
    except IndexError:
      ok = 1
    else:
      ok = 0
    assert ok, 'accessing bogus bit did not fail'
    try:
      gen[-1]
    except IndexError:
      ok = 1
    else:
      ok = 0
    assert ok, 'accessing bogus bit did not fail'


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
