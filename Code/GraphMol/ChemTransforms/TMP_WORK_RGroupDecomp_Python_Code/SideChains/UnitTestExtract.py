# $Id: UnitTestExtract.py 15142 2015-01-08 04:57:48Z gedecpe1 $
#
# Created by Greg Landrum, Dec 2008
#
#pylint: disable=W0403,C0111
import unittest, subprocess, os, sys
from collections import namedtuple
from rdkit import Chem
from io import BytesIO
try:
  from . import Extract  #@UnusedImport
except (ValueError, SystemError):
  import Extract  #@UnresolvedImport @Reimport
__file__ = Extract.__file__
ScriptResult = namedtuple('ScriptResult', 'retCode,stdout,stderr')

if sys.version_info[:2] < (2,7):
  def assertIn_(self,a,b):
    self.assertTrue(a in b)
  unittest.TestCase.assertIn = assertIn_

class TestCase(unittest.TestCase):
  def setUp(self):
    self.old_wd = os.getcwd()
    os.chdir(os.path.abspath(os.path.dirname(__file__)))
    if os.path.exists('testData/av.out.txt'):
      os.unlink('testData/av.out.txt')
    
  def tearDown(self):
    if os.path.exists('testData/av.out.txt'):
      os.unlink('testData/av.out.txt')
    os.chdir(self.old_wd)
  
  @staticmethod
  def runScript(cmd, debug=False):
    """ Run the script and return the exit code, stdin and stdout """
    cmd = list(cmd)
    if sys.executable not in cmd:
      cmd.insert(0, sys.executable)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         bufsize=0)
    stdout, stderr = p.stdout.read(), p.stderr.read()
    res = p.wait()
    p.stdout.close()
    p.stderr.close()
    #if debug:
    print(stdout)
    print(stderr)
    return ScriptResult(res, stdout, stderr)


  def test1(self):
    # Test SMILES core
    result = self.runScript(('Extract.py',
                          '--core=c1ccc2c(cncn2)c1',
                          '--coreFormat=smiles',
                          '--dataFormat=other',
                          '--outFormat=csv',
                          '--outF=testData/av.out.txt',
                          'testData/av.txt'
                          ), debug=False)
    self.assertEqual(result.retCode, 0)
    
    l1 = self.assertExtractResult(241, 37)
    self.assertEqual(l1[:4], ['NVP-BXN061-NX-1', 'COc1ccc(-c2ccc3ncnc(-c4cccc(C(=O)N5CCN(C(C)=O)CC5)c4)c3c2)cn1',
                              '0', '*c1ccc(nc1)OC'])

  def assertExtractResult(self, nrows, ncols, outfile='testData/av.out.txt'):
    self.assertTrue(os.path.exists(outfile))
    with open(outfile, 'r') as f:
      d = f.readlines()
    os.unlink(outfile)
    
    self.assertEqual(len(d), nrows)
    l1 = d[1].strip().split(',')
    self.assertEqual(len(l1), ncols)
    return l1
    
    
  def test2(self):
    # Test SMARTS core
    result = self.runScript(('Extract.py',
                          '--core=c1ccc2c(c[n,c]cn2)c1',
                          '--coreFormat=smarts',
                          '--dataFormat=other',
                          '--outFormat=csv',
                          '--outF=testData/av.out.txt',
                          'testData/av.txt'
                          ), debug=False)
    self.assertEqual(result.retCode, 0)

    l1 = self.assertExtractResult(241, 37)
    self.assertEqual(l1[:4], ['NVP-BXN061-NX-1', 'COc1ccc(-c2ccc3ncnc(-c4cccc(C(=O)N5CCN(C(C)=O)CC5)c4)c3c2)cn1',
                              '0', '*c1ccc(nc1)OC'])


  def test2b(self):
    # Test SMARTS core
    result = self.runScript(('Extract.py',
                          '--core=c1ccc2c(c[n,c]cn2)c1',
                          '--coreFormat=smarts',
                          '--dataFormat=other',
                          '--outFormat=csv',
                          '--outF=testData/av.out.txt',
                          'testData/av.txt'
                          ), debug=False)
    self.assertEqual(result.retCode, 0)

    l1 = self.assertExtractResult(241, 37)
    self.assertEqual(l1[:4], ['NVP-BXN061-NX-1', 'COc1ccc(-c2ccc3ncnc(-c4cccc(C(=O)N5CCN(C(C)=O)CC5)c4)c3c2)cn1',
                              '0', '*c1ccc(nc1)OC'])

    
  def test3(self):
    # Test core from SDF
    result = self.runScript(('Extract.py',
                          '--core=testData/core1.mol',
                          '--coreFormat=sdf',
                          '--dataFormat=other',
                          '--outFormat=csv',
                          '--outF=testData/av.out.txt',
                          'testData/av.txt'
                          ), debug=False)
    self.assertEqual(result.retCode, 0)

    l1 = self.assertExtractResult(241, 37)
    self.assertEqual(l1[:4], ['NVP-BXN061-NX-1', 'COc1ccc(-c2ccc3ncnc(-c4cccc(C(=O)N5CCN(C(C)=O)CC5)c4)c3c2)cn1',
                              '0', '*c1ccc(nc1)OC'])

    
  def test4(self):
    # Test core from avalon
    result = self.runScript(('Extract.py',
                          '--core=testData/combined.avalon',
                          '--coreFormat=avalon',
                          '--dataFormat=other',
                          '--outFormat=csv',
                          '--outF=testData/av.out.txt',
                          '--silent',
                          'testData/av.txt'
                          ), debug=False)
    self.assertEqual(result.retCode, 0)

    l1 = self.assertExtractResult(481, 37)
    self.assertEqual(l1[:4], ['NVP-BXN061-NX-1', 'COc1ccc(-c2ccc3ncnc(-c4cccc(C(=O)N5CCN(C(C)=O)CC5)c4)c3c2)cn1',
                              '0', '*c1ccc(nc1)OC'])

    
  def test5(self):
    # Test core labelled
    result = self.runScript(('Extract.py',
                          '--core=testData/dummy_in_scaffold.mol',
                          '--coreFormat=sdf',
                          '--dataFormat=other',
                          '--outFormat=csv',
                          '--outF=testData/av.out.txt',
                          '--silent',
                          '--labelledCores',
                          'testData/dummy_in_scaffold.txt'
                          ), debug=False)
    self.assertEqual(result.retCode, 0)

    l1 = self.assertExtractResult(4, 6)
    self.assertEqual(l1, ['mol1', 'Cc1c2c(c3n1CCCC3)n(C)c(=O)n(C)c2=O', '0', '*C', '*C', '*C', ])


  def test_GetMolDegenPoints(self):
    self.assertEqual(Extract.GetMolDegenPoints(Chem.MolFromSmiles('Nc1ncccc1')), set())
    self.assertEqual(Extract.GetMolDegenPoints(Chem.MolFromSmiles('Nc1ccccc1')), set([2, 3, 5, 6]))
    
    # cubane is a nice case:
    self.assertEqual(Extract.GetMolDegenPoints(Chem.MolFromSmiles('NC12C3C4C5C3C1C5C24')), set([2, 3, 5, 6, 7, 8]))

  def test_SymmetrizeSidechains(self):
    # indexing of arguments and results:  [mol][group] -> (attachIdx,smiles)
    smis = ['Nc1c(C)cccc1', 'Nc1ccccc1C']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    cores = [Chem.MolFromSmiles('Nc1ccccc1')]

    # Preparing the results from GetSidechains:
    d = Extract.GetSidechains(ms, cores)
    chains = [x[0] for x in d]

    # Starting point: non-symmetric chains:
    self.assertEqual(chains, [((2, '*C'),), ((6, '*C'),)])

    schains = Extract.SymmetrizeSidechains(ms, cores[0], chains)
    self.assertEqual(len(schains), len(chains))
    self.assertEqual(schains, [((2, '*C'),), ((2, '*C'),)])

    # When the core is nonsymmetric, this doesn't do anything:
    smis = ['Nc1c(C)nccc1', 'Nc1cnccc1C']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    cores = [Chem.MolFromSmiles('Nc1cnccc1')]
    chains = Extract.GetSidechains(ms, cores)
    schains = Extract.SymmetrizeSidechains(ms, cores[0], chains)
    self.assertTrue(schains is chains)

    # A more interesting variation, with two substituents per mol:
    smis = ['Nc1c(C)c(O)ccc1', 'Nc1cccc(O)c1C']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    cores = [Chem.MolFromSmiles('Nc1ccccc1')]

    d = Extract.GetSidechains(ms, cores)
    chains = [x[0] for x in d]
    schains = Extract.SymmetrizeSidechains(ms, cores[0], chains)
    self.assertEqual(len(schains), len(chains))
    self.assertEqual(chains, [((2, '*C'), (3, '*O')), ((5, '*O'), (6, '*C'))])
    self.assertEqual(schains, [((2, '*C'), (3, '*O')), ((2, '*C'), (3, '*O'))])

    # First "real" test case:
    smis = ['Nc1c(C)cccc1', 'Nc1ccccc1C', 'Nc1c(Cl)cccc1C']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    cores = [Chem.MolFromSmiles('Nc1ccccc1')]

    d = Extract.GetSidechains(ms, cores)
    chains = [x[0] for x in d]
    schains = Extract.SymmetrizeSidechains(ms, cores[0], chains)
    self.assertEqual(len(schains), len(chains))
    self.assertEqual(chains, [((2, '*C'),), ((6, '*C'),), ((2, '*Cl'), (6, '*C'))])
    self.assertEqual(schains, [((2, '*C'),), ((2, '*C'),), ((2, '*C'), (6, '*Cl'))])

    smis = ['Nc1c(C)cccc1', 'Nc1c(Cl)cccc1C', 'Nc1ccccc1C']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    cores = [Chem.MolFromSmiles('Nc1ccccc1')]
    d = Extract.GetSidechains(ms, cores)
    chains = [x[0] for x in d]
    schains = Extract.SymmetrizeSidechains(ms, cores[0], chains)
    self.assertEqual(len(schains), len(chains))
    self.assertEqual(chains, [((2, '*C'),), ((2, '*Cl'), (6, '*C')), ((6, '*C'),)])
    self.assertEqual(schains, [((2, '*C'),), ((2, '*C'), (6, '*Cl')), ((2, '*C'),)])

    # This case demonstrates a wart in the current procedure. This is correct
    # for the current implementation, but should ideally behave differently:
    smis = ['Nc1c(C)cccc1', 'Nc1c(C)ccc(Cl)c1', 'Nc1cccc(Cl)c1']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    cores = [Chem.MolFromSmiles('Nc1ccccc1')]
    d = Extract.GetSidechains(ms, cores)
    chains = [x[0] for x in d]
    schains = Extract.SymmetrizeSidechains(ms, cores[0], chains)
    self.assertEqual(len(schains), len(chains))
    self.assertEqual(chains, [((2, '*C'),), ((2, '*C'), (5, '*Cl')), ((5, '*Cl'),)])
    self.assertEqual(schains, [((2, '*C'),), ((2, '*C'), (5, '*Cl')), ((3, '*Cl'),)])

    # really should be: [((2, '*C'),), ((2, '*C'), (5, '*Cl')), ((5, '*Cl'),)]

    # Robust w.r.t to multiple attachment points in sidechain:
    options, _ = Extract.parser.parse_args([])
    options.silent = True
    smis = ['C2C[C]c1ccccc1C2', 'C2C[N]c1ccccc1C2', 'NC(Cc1ccc(O)cc1)C(O)=O']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    cores = [Chem.MolFromSmiles('c1ccccc1C')]
    d = Extract.GetSidechains(ms, cores, options)
    self.assertEqual(d, [(((None, 'core matches multiple times'),),),
                         (((None, 'multiple attachment points in sidechain'),),),
                         (((2, '*O'), (6, '*C(N)C(O)=O')),)])
    chains = [x[0] for x in d]
    schains = Extract.SymmetrizeSidechains(ms, cores[0], chains)
    self.assertEqual(len(schains), len(chains))
    self.assertEqual(schains, [((None, 'core matches multiple times'),),
                               ((None, 'multiple attachment points in sidechain'),),
                               ((2, '*O'), (6, '*C(N)C(O)=O'))])

  def test_GetSidechains(self):
    class Silent(object): 
      """ Suppress error messages """
      silent = True
      rejectDoubleAttachments = False
    silent = Silent()
    
    # indexing of results:   [mol][core][group] -> (attachIdx,smiles)
    smis = ['Nc1c(C)cccc1', 'Nc1c(O)cccc1']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    cores = [Chem.MolFromSmiles('Nc1ccccc1')]
    d = Extract.GetSidechains(ms, cores)
    self.assertEqual(len(d), len(ms))
    self.assertEqual(len(d[0]), len(cores))
    self.assertEqual(len(d[0][0]), len(d[1][0]))
    idx0, smi0 = d[0][0][0]
    idx1, smi1 = d[1][0][0]
    self.assertEqual(idx0, idx1)
    self.assertEqual(smi0, '*C')
    self.assertEqual(smi1, '*O')
  
    # This doesn't symmetrize things, so for symmetric scaffolds we aren't
    # guaranteed to get the same R labels:
    smis = ['Nc1c(C)cccc1', 'Nc1ccccc1C']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    cores = [Chem.MolFromSmiles('Nc1ccccc1')]
    d = Extract.GetSidechains(ms, cores)
    idx0, smi0 = d[0][0][0]
    idx1, smi1 = d[1][0][0]
    self.assertNotEqual(idx0, idx1)
    self.assertEqual(smi0, '*C')
    self.assertEqual(smi1, '*C')
    
    # If the core doesn't match, empty lists are returned:
    smis = ['Nc1c(C)cccc1']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    cores = [Chem.MolFromSmiles('Nc1cnccc1'), Chem.MolFromSmiles('Nc1cncnc1')]
    d = Extract.GetSidechains(ms, cores, options=silent)
    self.assertEqual(len(d) , len(ms))
    self.assertEqual(len(d[0]) , len(cores))
    self.assertEqual(d[0][0], ((None, 'no core matches'),))
    self.assertEqual(d[0][1], ((None, 'no core matches'),))
    
    # And we're robust w.r.t. bogus cores:
    smis = ['Nc1c(C)cccc1']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    cores = [Chem.MolFromSmiles('Nc1ccccc1'), None]
    d = Extract.GetSidechains(ms, cores)
    self.assertEqual(len(d) , len(ms))
    self.assertEqual(len(d[0]) , len(cores))
    self.assertEqual(len(d[0][0]), 1)
    self.assertEqual(len(d[0][1]), 2)
    self.assertEqual(d[0][0], ((2, '*C'),))
    self.assertEqual(d[0][1], (None, 'bad core'))
    
    # and bogus molecules:
    ms = [None, Chem.MolFromSmiles('Nc1c(C)cccc1')]
    cores = [Chem.MolFromSmiles('Nc1cnccc1')]
    d = Extract.GetSidechains(ms, cores, options=silent)
    self.assertEqual(len(d) , len(ms))
    self.assertEqual(len(d[0]) , len(cores))
    self.assertEqual(len(d[0][0]), 2)
    self.assertEqual(d[0][0], (None, 'bad molecule'))
    
    # by default we accept doubly attached sidechains:
    options, _ = Extract.parser.parse_args([])
    options.silent = True
    smis = ['C2Nc1c(C2)cccc1']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    cores = [Chem.MolFromSmiles('Cc1c(N)cccc1')]
    d = Extract.GetSidechains(ms, cores, options=silent)
    self.assertEqual(d[0], (((0, '*C*'),),))

    # but this can be toggled off with the options construct:
    options, _ = Extract.parser.parse_args([])
    options.silent = True
    options.rejectDoubleAttachments = True
    d = Extract.GetSidechains(ms, cores, options=options)
    self.assertEqual(d[0], (((None, 'multiple attachment points in sidechain'),),))
    
    
  def test_ProcessCoreLabels(self):
    options, _ = Extract.parser.parse_args([])
    options.silent = True
    options.labelledCores = 1
    core = Chem.MolFromSmiles('c1ccccc1[2*]')
    cores = [core]
    Extract.ProcessCoreLabels(cores, options)
    self.assertEqual(Chem.MolToSmiles(cores[0]), 'c1ccccc1')
    self.assertEqual(cores[0].GetAtomWithIdx(5).HasProp('_RLabel'), 1)
    self.assertEqual(cores[0].GetAtomWithIdx(5).GetProp('_RLabel'), '2')
    core = Chem.MolFromSmarts('c1ccccc1[*:2]')
    cores = [core]
    Extract.ProcessCoreLabels(cores, options)
    self.assertEqual(Chem.MolToSmiles(cores[0]), 'c1ccccc1')
    self.assertEqual(cores[0].GetAtomWithIdx(5).HasProp('_RLabel'), 1)
    self.assertEqual(cores[0].GetAtomWithIdx(5).GetProp('_RLabel'), '2')

    core = Chem.MolFromSmarts('c1ccc[c,n]c1[*:2]')
    cores = [core]
    Extract.ProcessCoreLabels(cores, options)
    self.assertEqual(Chem.MolToSmiles(cores[0]), 'c1ccccc1')
  
  def test_RunDecomposition(self):
    class Silent(object): 
      """ Suppress error messages """
      silent = True
      rejectDoubleAttachments = False
    silent = Silent()

    options, _ = Extract.parser.parse_args([])
    options.silent = True
    options.labelledCores = 1
    options.requireLabels = 1
    smis = ['Nc1ccccc1', 'c1cc(F)ccc1', 'c1cc(F)c(C)cc1']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    options.cores = [Chem.MolFromSmiles('[1*]c1ccccc1')]
    _, d = Extract.RunDecomposition(ms, options)

    self.assertEqual(len(d), len(ms))
    self.assertEqual(len(d[0]), len(options.cores))
    self.assertEqual(len(d[0][0]), 1)
    self.assertEqual(d[0][0][0], (1, '*N'))
    self.assertEqual(d[1][0][0], (1, '*F'))
    self.assertEqual(d[2], (None, 'substitution at unmarked position'))

    options, _ = Extract.parser.parse_args([])
    options.silent = True
    options.labelledCores = 1
    options.requireLabels = 1
    smis = ['Nc1cc(F)ccc1', 'c1cc(F)cc(F)c1', 'c1c(C)cc(C)cc1', 'c1c(C)c(C)ccc1']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    options.cores = [Chem.MolFromSmiles('[1*]c1cc([2*])ccc1')]
    _, d = Extract.RunDecomposition(ms, options)

    self.assertEqual(len(d), len(ms))
    self.assertEqual(len(d[0]), len(options.cores))
    self.assertEqual(len(d[0][0]), 2)
    self.assertEqual(d[0][0], ((2, '*N'), (1, '*F')))
    self.assertEqual(d[1][0], ((1, '*F'), (2, '*F')))
    self.assertEqual(d[2][0], ((1, '*C'), (2, '*C')))
    self.assertEqual(d[3], (None, 'substitution at unmarked position'))

    options, _ = Extract.parser.parse_args([])
    options.silent = True
    options.labelledCores = 1
    options.requireLabels = 1
    smis = ['CCN3C=Nc2c(Nc1cccc(Cl)c1)nc(nc23)N4CCC(O)CC4', 'COc4ccc(CNC2=NNc3ncnc(Nc1cccc(Cl)c1)c23)cc4Cl']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    options.cores = [Chem.MolFromSmiles('[1*]c1cccc([2*])c1')]
    _, d = Extract.RunDecomposition(ms, options)

    self.assertEqual(len(d), len(ms))
    self.assertEqual(len(d[0]), len(options.cores))
    self.assertEqual(d[0][0], ((2, '*Nc1nc(nc2c1N=CN2CC)N3CCC(O)CC3'), (1, '*Cl')))

    self.assertEqual(d[1], (None, 'core matches multiple times'))

    options, _ = Extract.parser.parse_args([])
    options.silent = True
    options.labelledCores = 1
    options.requireLabels = 1
    smis = ['CC1CC1', 'C1CCC1']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    options.cores = [Chem.MolFromSmiles('C1CC1[1*]')]
    _, d = Extract.RunDecomposition(ms, options)

    self.assertEqual(len(d), len(ms))
    self.assertEqual(len(d[0]), len(options.cores))
    self.assertEqual(d[0][0], ((1, '*C'),))
    self.assertEqual(d[1], (None, 'no core matches'))

    options, _ = Extract.parser.parse_args([])
    options.silent = True
    options.labelledCores = 1
    options.nonLabelledSubstituentHandling = 'FAIL'
    smis = ['CCC1CC1', 'C1C(Cl)C1CC']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    options.cores = [Chem.MolFromSmiles('C1CC1C[1*]')]
    _, d = Extract.RunDecomposition(ms, options)
    self.assertEqual(len(d), len(ms))
    self.assertEqual(len(d[0]), len(options.cores))
    self.assertEqual(d[0][0], ((1, '*C'),))
    self.assertEqual(d[1], (None, 'substitution at unmarked position'))

    options, _ = Extract.parser.parse_args([])
    options.silent = True
    options.labelledCores = 1
    options.nonLabelledSubstituentHandling = 'PASS'
    smis = ['CCC1CC1', 'C1C(Cl)C1CC']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    options.cores = [Chem.MolFromSmiles('C1CC1C[1*]')]
    _, d = Extract.RunDecomposition(ms, options)
    self.assertEqual(len(d), len(ms))
    self.assertEqual(len(d[0]), len(options.cores))
    self.assertEqual(d[0][0], ((1, '*C'),))
    self.assertEqual(d[1][0], ((101, '*Cl'), (1, '*C')))


    options, _ = Extract.parser.parse_args([])
    options.silent = True
    options.labelledCores = 1
    options.nonLabelledSubstituentHandling = 'IGNORE'
    smis = ['CCC1CC1', 'C1C(Cl)C1CC']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    options.cores = [Chem.MolFromSmiles('C1CC1C[1*]')]
    _, d = Extract.RunDecomposition(ms, options)
    self.assertEqual(len(d), len(ms))
    self.assertEqual(len(d[0]), len(options.cores))
    self.assertEqual(d[0][0], ((1, '*C'),))
    self.assertEqual(d[1][0], ((1, '*C'),))



  def test_CreateOutput(self):
    # repeat a test case we used above:
    options, _ = Extract.parser.parse_args([])
    options.silent = True
    
    options.labelledCores = 1
    options.requireLabels = 1
    smis = ['CCN3C=Nc2c(Nc1cccc(Cl)c1)nc(nc23)N4CCC(O)CC4', 'COc4ccc(CNC2=NNc3ncnc(Nc1cccc(Cl)c1)c23)cc4Cl']
    ms = [Chem.MolFromSmiles(x) for x in smis]
    options.cores = [Chem.MolFromSmiles('[1*]c1cccc([2*])c1')]
    mols, d = Extract.RunDecomposition(ms, options)
    sio = BytesIO()
    options.outFile = sio
    options.outputFormat = 'csv'
    Extract.CreateOutput(mols, d, options)
    self.assertIn('Mol_1', sio.getvalue())
    self.assertIn('compound_id,smiles,Core,R1,R2', sio.getvalue())

    # without labels:
    ms.append(Chem.MolFromSmiles('[C]1CCCc2c1cccc2'))
    sio = BytesIO()
    options.outFile = sio
    options.labelledCores = 0
    options.requireLabels = 0
    options.cores = [Chem.MolFromSmiles('c1ccccc1')]
    mols, d = Extract.RunDecomposition(ms, options)
    Extract.CreateOutput(mols, d, options)
    self.assertIn('compound_id,smiles,Core,R1,R5', sio.getvalue())
    self.assertIn('Mol_3,[C]1CCCc2ccccc21,ERROR: multiple attachment points in sidechain', sio.getvalue())
  
    # robust r.s.t bogus core when includeCoreSmiles is ON
    sio = BytesIO()
    options.outFile = sio
    options.includeCoreSmiles = '1'
    ms.append(Chem.MolFromSmiles('CCC'))
    mols, d = Extract.RunDecomposition(ms, options)
    Extract.CreateOutput(mols, d, options)
    self.assertIn('Mol_4,CCC,ERROR: no core matches', sio.getvalue())


#if __name__ == '__main__':
unittest.main()

