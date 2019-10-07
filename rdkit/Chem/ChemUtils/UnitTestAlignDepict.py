from contextlib import contextmanager
import sys
import unittest

from rdkit import Chem
from rdkit.Chem.ChemUtils.AlignDepict import initParser, processArgs, AlignDepict
from io import StringIO


class TestCase(unittest.TestCase):

    def test1(self):
        parser = initParser()
        with outputRedirect() as (out, err):
            args = parser.parse_args('--smiles CC CCC'.split())
            args.outF = out
            processArgs(args)
        self.assertIn('RDKit', out.getvalue())
        self.assertIn('2D', out.getvalue())
        self.assertEqual(err.getvalue(), '')

    def test_AlignDepict(self):
        mol = Chem.MolFromSmiles('CNC')
        core = Chem.MolFromSmiles('CC')
        pattern = Chem.MolFromSmarts('CCC')
        self.assertRaises(ValueError, AlignDepict, mol, core, pattern)

        pattern = Chem.MolFromSmarts('CN')
        self.assertRaises(ValueError, AlignDepict, mol, core, pattern)

        pattern = Chem.MolFromSmarts('CC')
        self.assertRaises(ValueError, AlignDepict, mol, core, pattern)

        pattern = Chem.MolFromSmarts('CC')
        self.assertRaises(ValueError, AlignDepict, mol, core, pattern)

        mol = Chem.MolFromSmiles('CCC')
        Chem.rdDepictor.Compute2DCoords(core)
        AlignDepict(mol, core, pattern)

        mol = Chem.MolFromSmiles('CNC')
        AlignDepict(mol, core, pattern, acceptFailure=True)


@contextmanager
def outputRedirect():
    """ Redirect standard output and error to String IO and return """
    try:
        _stdout, _stderr = sys.stdout, sys.stderr
        sys.stdout = sStdout = StringIO()
        sys.stderr = sStderr = StringIO()
        yield (sStdout, sStderr)
    finally:
        sys.stdout, sys.stderr = _stdout, _stderr


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
