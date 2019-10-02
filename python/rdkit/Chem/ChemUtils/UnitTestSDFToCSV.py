from contextlib import contextmanager
import os
import sys
import unittest

from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem.ChemUtils.SDFToCSV import Convert, initParser
from io import StringIO


class TestCase(unittest.TestCase):

    def test1(self):
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
        lines = [line for line in txt.split('\n') if line.strip() != '']
        self.assertTrue(len(lines) == 6, 'bad num lines: %d' % len(lines))
        line0 = lines[0].split(',')
        self.assertEqual(len(line0), 20)
        self.assertTrue(line0[0] == 'AMW')
        self.assertTrue(line0[1] == 'SMILES')

    def test_parser(self):
        parser = initParser()
        # User want's help
        with self.assertRaises(SystemExit), outputRedirect() as (out, err):
            parser.parse_args(['-h'])
        self.assertNotEqual(out.getvalue(), '')
        self.assertEqual(err.getvalue(), '')

        # Missing input file
        with self.assertRaises(SystemExit), outputRedirect() as (out, err):
            parser.parse_args([])
        self.assertEqual(out.getvalue(), '')
        self.assertNotEqual(err.getvalue(), '')

        # Input file doesn't exist
        with self.assertRaises(SystemExit), outputRedirect() as (out, err):
            parser.parse_args(['incorrectFilename'])
        self.assertEqual(out.getvalue(), '')
        self.assertNotEqual(err.getvalue(), '')


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
