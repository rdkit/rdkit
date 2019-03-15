#
# Copyright (C) 2018 Susan H. Leung
#         All Rights Reserved
#
from rdkit import RDConfig
import os
import sys
import math
import unittest
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Geometry import rdGeometry as geom
from rdkit.Chem.rdchem import Atom
from rdkit.Chem.MolStandardize import rdMolStandardize


class TestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test1Cleanup(self):
        mol = Chem.MolFromSmiles("CCC(=O)O[Na]")
        nmol = rdMolStandardize.Cleanup(mol)
        self.assertEqual(Chem.MolToSmiles(nmol), "CCC(=O)[O-].[Na+]")

    def test2StandardizeSmiles(self):
        self.assertEqual(rdMolStandardize.StandardizeSmiles("CCC(=O)O[Na]"), "CCC(=O)[O-].[Na+]")
        # should get ValueError


#    rdMolStandardize.StandardizeSmiles("C1CCC1C(=O)O.Na")


    def test3FragmentParent(self):
        mol = Chem.MolFromSmiles("[Na]OC(=O)c1ccccc1")
        nmol = rdMolStandardize.FragmentParent(mol)
        self.assertEqual(Chem.MolToSmiles(nmol), "O=C([O-])c1ccccc1")

    def test4ChargeParent(self):
        mol = Chem.MolFromSmiles("C[NH+](C)(C).[Cl-]")
        nmol = rdMolStandardize.ChargeParent(mol)
        self.assertEqual(Chem.MolToSmiles(nmol), "CN(C)C")

    def test4Normalize(self):
        mol = Chem.MolFromSmiles("C[N+](C)=C\C=C\[O-]")
        nmol = rdMolStandardize.Normalize(mol)
        self.assertEqual(Chem.MolToSmiles(nmol), "CN(C)C=CC=O")

    def test4Reionize(self):
        mol = Chem.MolFromSmiles("C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O")
        nmol = rdMolStandardize.Reionize(mol)
        self.assertEqual(Chem.MolToSmiles(nmol), "O=S(O)c1ccc(S(=O)(=O)[O-])cc1")

    def test5Metal(self):
        mol = Chem.MolFromSmiles("C1(CCCCC1)[Zn]Br")
        md = rdMolStandardize.MetalDisconnector()
        nm = md.Disconnect(mol)
        #    Metal.MetalDisconnector.Disconnect(mol)
        self.assertEqual(Chem.MolToSmiles(nm), "[Br-].[CH-]1CCCCC1.[Zn+2]")

        # test user defined metal_nof
        md.SetMetalNof(
          Chem.MolFromSmarts(
            "[Li,K,Rb,Cs,Fr,Be,Mg,Ca,Sr,Ba,Ra,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Al,Ga,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi]~[N,O,F]"
          ))
        mol2 = Chem.MolFromSmiles("CCC(=O)O[Na]")
        nm2 = md.Disconnect(mol2)
        self.assertEqual(Chem.MolToSmiles(nm2), "CCC(=O)O[Na]")

    def test6Charge(self):
        mol = Chem.MolFromSmiles("C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O")
        # instantiate with default acid base pair library
        reionizer = rdMolStandardize.Reionizer()
        nm = reionizer.reionize(mol)
        self.assertEqual(Chem.MolToSmiles(nm), "O=S(O)c1ccc(S(=O)(=O)[O-])cc1")

        # try reionize with another acid base pair library without the right
        # pairs
        abfile = os.path.join(RDConfig.RDDataDir, 'MolStandardize', 'acid_base_pairs2.txt')
        reionizer2 = rdMolStandardize.Reionizer(abfile)
        nm2 = reionizer2.reionize(mol)
        self.assertEqual(Chem.MolToSmiles(nm2), "O=S([O-])c1ccc(S(=O)(=O)O)cc1")

        # test Uncharger
        uncharger = rdMolStandardize.Uncharger()
        mol3 = Chem.MolFromSmiles("O=C([O-])c1ccccc1")
        nm3 = uncharger.uncharge(mol3)
        self.assertEqual(Chem.MolToSmiles(nm3), "O=C(O)c1ccccc1")

    def test7Fragment(self):
        fragremover = rdMolStandardize.FragmentRemover()
        mol = Chem.MolFromSmiles("CN(C)C.Cl.Cl.Br")
        nm = fragremover.remove(mol)
        self.assertEqual(Chem.MolToSmiles(nm), "CN(C)C")

        lfragchooser = rdMolStandardize.LargestFragmentChooser()
        mol2 = Chem.MolFromSmiles("[N+](=O)([O-])[O-].[CH3+]")
        nm2 = lfragchooser.choose(mol2)
        self.assertEqual(Chem.MolToSmiles(nm2), "O=[N+]([O-])[O-]")

        lfragchooser2 = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
        nm3 = lfragchooser2.choose(mol2)
        self.assertEqual(Chem.MolToSmiles(nm3), "[CH3+]")

        fragremover = rdMolStandardize.FragmentRemover(skip_if_all_match=True)
        mol = Chem.MolFromSmiles("[Na+].Cl.Cl.Br")
        nm = fragremover.remove(mol)
        self.assertEqual(nm.GetNumAtoms(), mol.GetNumAtoms())

    def test8Normalize(self):
        normalizer = rdMolStandardize.Normalizer()
        mol = Chem.MolFromSmiles("C[n+]1ccccc1[O-]")
        nm = normalizer.normalize(mol)
        self.assertEqual(Chem.MolToSmiles(nm), "Cn1ccccc1=O")

    def test9Validate(self):
        vm = rdMolStandardize.RDKitValidation()
        mol = Chem.MolFromSmiles("CO(C)C", sanitize=False)
        msg = vm.validate(mol)
        self.assertEqual(len(msg), 1)
        self.assertEqual
        ("""INFO: [ValenceValidation] Explicit valence for atom # 1 O, 3, is greater than permitted""",
         msg[0])

        vm2 = rdMolStandardize.MolVSValidation([rdMolStandardize.FragmentValidation()])
        # with no argument it also works
        #    vm2 = rdMolStandardize.MolVSValidation()
        mol2 = Chem.MolFromSmiles("COc1cccc(C=N[N-]C(N)=O)c1[O-].O.O.O.O=[U+2]=O")
        msg2 = vm2.validate(mol2)
        self.assertEqual(len(msg2), 1)
        self.assertEqual
        ("""INFO: [FragmentValidation] water/hydroxide is present""", msg2[0])

        vm3 = rdMolStandardize.MolVSValidation()
        mol3 = Chem.MolFromSmiles("C1COCCO1.O=C(NO)NO")
        msg3 = vm3.validate(mol3)
        self.assertEqual(len(msg3), 2)
        self.assertEqual
        ("""INFO: [FragmentValidation] 1,2-dimethoxyethane is present""", msg3[0])
        self.assertEqual
        ("""INFO: [FragmentValidation] 1,4-dioxane is present""", msg3[1])

        atomic_no = [6, 7, 8]
        allowed_atoms = [Atom(i) for i in atomic_no]
        vm4 = rdMolStandardize.AllowedAtomsValidation(allowed_atoms)
        mol4 = Chem.MolFromSmiles("CC(=O)CF")
        msg4 = vm4.validate(mol4)
        self.assertEqual(len(msg4), 1)
        self.assertEqual
        ("""INFO: [AllowedAtomsValidation] Atom F is not in allowedAtoms list""", msg4[0])

        atomic_no = [9, 17, 35]
        disallowed_atoms = [Atom(i) for i in atomic_no]
        vm5 = rdMolStandardize.DisallowedAtomsValidation(disallowed_atoms)
        mol5 = Chem.MolFromSmiles("CC(=O)CF")
        msg5 = vm4.validate(mol5)
        self.assertEqual(len(msg5), 1)
        self.assertEqual
        ("""INFO: [DisallowedAtomsValidation] Atom F is in disallowedAtoms list""", msg5[0])

        msg6 = rdMolStandardize.ValidateSmiles("ClCCCl.c1ccccc1O")
        self.assertEqual(len(msg6), 1)
        self.assertEqual
        ("""INFO: [FragmentValidation] 1,2-dichloroethane is present""", msg6[0])


if __name__ == "__main__":
    unittest.main()
