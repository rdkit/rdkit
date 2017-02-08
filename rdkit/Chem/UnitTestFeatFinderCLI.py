from contextlib import contextmanager
import os
import sys
import unittest

from rdkit import RDConfig
from rdkit.Chem import FeatFinderCLI
from rdkit.six.moves import cStringIO as StringIO


class TestCase(unittest.TestCase):

  def test_FeatFinderCLI(self):
    smilesFile = os.path.join(RDConfig.RDDataDir, 'NCI', 'first_5K.smi')
    featureFile = os.path.join(RDConfig.RDCodeDir, 'Chem', 'Pharm2D', 'test_data',
                               'BaseFeatures.fdef')
    parser = FeatFinderCLI.initParser()
    cmd = '-n 10 {0} {1}'.format(featureFile, smilesFile)
    with outputRedirect() as (out, err):
      args = parser.parse_args(cmd.split())
      FeatFinderCLI.processArgs(args, parser)
    self.assertIn('Mol-1', out.getvalue())
    self.assertIn('Acceptor-SingleAtomAcceptor', out.getvalue())
    self.assertIn('C(1)', out.getvalue())
    self.assertNotIn('Mol-11', out.getvalue())
    self.assertEqual(err.getvalue(), '')

    cmd = '-n 2 -r {0} {1}'.format(featureFile, smilesFile)
    with outputRedirect() as (out, err):
      args = parser.parse_args(cmd.split())
      FeatFinderCLI.processArgs(args, parser)
    self.assertIn('Mol-1', out.getvalue())
    self.assertIn('Acceptor-SingleAtomAcceptor:', out.getvalue())
    self.assertIn('2, 3, 4', out.getvalue())
    self.assertNotIn('Mol-3', out.getvalue())
    self.assertEqual(err.getvalue(), '')

  def test_FeatFinderCLIexceptions(self):
    smilesFile = os.path.join(RDConfig.RDDataDir, 'NCI', 'first_5K.smi')
    featureFile = os.path.join(RDConfig.RDCodeDir, 'Chem', 'Pharm2D', 'test_data',
                               'BaseFeatures.fdef')
    parser = FeatFinderCLI.initParser()
    cmd = '-n 10 {0} {1}'.format(smilesFile, smilesFile)
    with self.assertRaises(SystemExit), outputRedirect() as (_, err):
      args = parser.parse_args(cmd.split())
      FeatFinderCLI.processArgs(args, parser)
    self.assertIn('error', err.getvalue())

    cmd = '-n 10 {0} {1}'.format(featureFile, 'incorrectFilename')
    with self.assertRaises(SystemExit), outputRedirect() as (_, err):
      args = parser.parse_args(cmd.split())
      FeatFinderCLI.processArgs(args, parser)
    self.assertIn('error', err.getvalue())


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
