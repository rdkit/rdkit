import os
import sys
import unittest
from rdkit import RDConfig
from rdkit import Chem
from rdkit.Chem import ChemicalForceFields

def get_test_mol():
    fName = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol',
            'ForceFieldHelpers', 'UFF', 'test_data', 'CH4.mol')  # CH4.mol, toluene.mol
    return Chem.MolFromMolFile(fName, True, False)

class TestCase(unittest.TestCase):
  def testANIForceField(self):
    mol = get_test_mol()
    ff = ChemicalForceFields.ANIGetMoleculeForceField(mol, 'ANI-1ccx')
    import time
    start = time.time()
    E_start = ff.CalcEnergy()
    end = time.time()
    print('E_start =', E_start, end - start, 's')
    start = time.time()
    ff.Minimize(maxIts=20)
    E_min = ff.CalcEnergy()
    end = time.time()
    print('E_min =', E_min, end - start, 's')
    self.assertAlmostEqual(E_start, -40.0552925, delta=1e-4)  # CH4
    self.assertAlmostEqual(E_min, -40.4576525, delta=1e-4)  # CH4
    #self.assertAlmostEqual(E_start, -271.2003578, delta=1e-4)  # toluene
    #self.assertAlmostEqual(E_min, -271.2028856, delta=1e-4)  # toluene

  def _testANIForceFieldSpeed(self):
    from rdkit.Chem import AllChem
    for nC in range(20, 40):
        mol = Chem.MolFromSmiles('C' * nC)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        ff = ChemicalForceFields.ANIGetMoleculeForceField(mol, "ANI-1ccx")
        import time
        start = time.time()
        E_start = ff.CalcEnergy()
        end = time.time()
        sp_time = end - start
        start = time.time()
        ff.Minimize(maxIts=20)
        E_min = ff.CalcEnergy()
        end = time.time()
        opt_time = end - start
        print(nC, mol.GetNumAtoms(), sp_time, opt_time)

if __name__ == '__main__':
  unittest.main()