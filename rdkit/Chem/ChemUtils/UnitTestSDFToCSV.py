import unittest

from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem.ChemUtils.SDFToCSV import Convert
from rdkit.six.moves import cStringIO as StringIO


class TestCase(unittest.TestCase):
  def test1(self):
    import os
    fName = os.path.join(RDConfig.RDDataDir, 'NCI', 'first_200.props.sdf')
    suppl = Chem.SDMolSupplier(fName)
    io = StringIO()
    try:
      Convert(suppl, io)
    except Exception:
      import traceback
      traceback.print_exc()
      self.fail('conversion failed')
    txt = io.getvalue()
    lines = txt.split('\n')
    if not lines[-1]:
      del lines[-1]
    self.assertTrue(len(lines) == 201, 'bad num lines: %d' % len(lines))
    line0 = lines[0].split(',')
    self.assertEqual(len(line0), 20)
    self.assertTrue(line0[0] == 'SMILES')

  def test2(self):
    import os
    from rdkit.six.moves import cStringIO as StringIO  #@UnresolvedImport #pylint: disable=F0401
    fName = os.path.join(RDConfig.RDDataDir, 'NCI', 'first_200.props.sdf')
    suppl = Chem.SDMolSupplier(fName)
    io = StringIO()
    try:
      Convert(suppl, io, keyCol='AMW', stopAfter=5)
    except Exception:
      import traceback
      traceback.print_exc()
      self.fail('conversion failed')
    txt = io.getvalue()
    lines = txt.split('\n')
    if not lines[-1]:
      del lines[-1]
    self.assertTrue(len(lines) == 6, 'bad num lines: %d' % len(lines))
    line0 = lines[0].split(',')
    self.assertEqual(len(line0), 20)
    self.assertTrue(line0[0] == 'AMW')
    self.assertTrue(line0[1] == 'SMILES')

if __name__ == '__main__':  # pragma: nocover
    unittest.main()
