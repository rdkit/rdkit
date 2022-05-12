#
# Copyright (C) 2018 Susan H. Leung
#         All Rights Reserved
#
from rdkit import RDConfig
import os
import sys
import math
from datetime import datetime, timedelta
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

  def test3Parents(self):
    mol = Chem.MolFromSmiles("[Na]OC(=O)c1ccccc1")
    nmol = rdMolStandardize.FragmentParent(mol)
    self.assertEqual(Chem.MolToSmiles(nmol), "O=C([O-])c1ccccc1")

    mol = Chem.MolFromSmiles("C[NH+](C)(C).[Cl-]")
    nmol = rdMolStandardize.ChargeParent(mol)
    self.assertEqual(Chem.MolToSmiles(nmol), "CN(C)C")

    mol = Chem.MolFromSmiles("[O-]CCCC=CO.[Na+]")
    nmol = rdMolStandardize.TautomerParent(mol)
    self.assertEqual(Chem.MolToSmiles(nmol), "O=CCCCC[O-].[Na+]")
    nmol = rdMolStandardize.TautomerParent(mol, skipStandardize=True)
    # same answer because of the standardization at the end
    self.assertEqual(Chem.MolToSmiles(nmol), "O=CCCCC[O-].[Na+]")

    mol = Chem.MolFromSmiles("C[C@](F)(Cl)C/C=C/[C@H](F)Cl")
    nmol = rdMolStandardize.StereoParent(mol)
    self.assertEqual(Chem.MolToSmiles(nmol), "CC(F)(Cl)CC=CC(F)Cl")

    mol = Chem.MolFromSmiles("[12CH3][13CH3]")
    nmol = rdMolStandardize.IsotopeParent(mol)
    self.assertEqual(Chem.MolToSmiles(nmol), "CC")

    mol = Chem.MolFromSmiles("[Na]Oc1c([12C@H](F)Cl)c(O[2H])c(C(=O)O)cc1CC=CO")
    nmol = rdMolStandardize.SuperParent(mol)
    self.assertEqual(Chem.MolToSmiles(nmol), "O=CCCc1cc(C(=O)O)c(O)c(C(F)Cl)c1O")
    mol = Chem.MolFromSmiles("[Na]Oc1c([12C@H](F)Cl)c(O[2H])c(C(=O)O)cc1CC=CO")
    nmol = rdMolStandardize.SuperParent(mol, skipStandardize=True)
    self.assertEqual(Chem.MolToSmiles(nmol), "O=CCCc1cc(C(=O)[O-])c(O)c(C(F)Cl)c1O.[Na+]")

  def test4Normalize(self):
    mol = Chem.MolFromSmiles(r"C[N+](C)=C\C=C\[O-]")
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
    abfile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MolStandardize', 'test_data',
                          'acid_base_pairs2.txt')
    reionizer2 = rdMolStandardize.Reionizer(abfile)
    nm2 = reionizer2.reionize(mol)
    self.assertEqual(Chem.MolToSmiles(nm2), "O=S([O-])c1ccc(S(=O)(=O)O)cc1")

    # test Uncharger
    uncharger = rdMolStandardize.Uncharger()
    mol3 = Chem.MolFromSmiles("O=C([O-])c1ccccc1")
    nm3 = uncharger.uncharge(mol3)
    self.assertEqual(Chem.MolToSmiles(nm3), "O=C(O)c1ccccc1")

    # test canonical Uncharger
    uncharger = rdMolStandardize.Uncharger(canonicalOrder=False)
    mol3 = Chem.MolFromSmiles("C[N+](C)(C)CC(C(=O)[O-])CC(=O)[O-]")
    nm3 = uncharger.uncharge(mol3)
    self.assertEqual(Chem.MolToSmiles(nm3), "C[N+](C)(C)CC(CC(=O)[O-])C(=O)O")
    uncharger = rdMolStandardize.Uncharger(canonicalOrder=True)
    nm3 = uncharger.uncharge(mol3)
    self.assertEqual(Chem.MolToSmiles(nm3), "C[N+](C)(C)CC(CC(=O)O)C(=O)[O-]")

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

    smi3 = "CNC[C@@H]([C@H]([C@@H]([C@@H](CO)O)O)O)O.c1cc2c(cc1C(=O)O)oc(n2)c3cc(cc(c3)Cl)Cl"

    lfParams = rdMolStandardize.CleanupParameters()
    lfrag_params = rdMolStandardize.LargestFragmentChooser(lfParams)
    mol3 = Chem.MolFromSmiles(smi3)
    lfrag3 = lfrag_params.choose(mol3)
    self.assertEqual(Chem.MolToSmiles(lfrag3), "CNC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO")

    lfParams = rdMolStandardize.CleanupParameters()
    lfParams.largestFragmentChooserCountHeavyAtomsOnly = True
    lfrag_params = rdMolStandardize.LargestFragmentChooser(lfParams)
    mol3 = Chem.MolFromSmiles(smi3)
    lfrag3 = lfrag_params.choose(mol3)
    self.assertEqual(Chem.MolToSmiles(lfrag3), "O=C(O)c1ccc2nc(-c3cc(Cl)cc(Cl)c3)oc2c1")

    lfParams = rdMolStandardize.CleanupParameters()
    lfParams.largestFragmentChooserUseAtomCount = False
    lfrag_params = rdMolStandardize.LargestFragmentChooser(lfParams)
    mol3 = Chem.MolFromSmiles(smi3)
    lfrag3 = lfrag_params.choose(mol3)
    self.assertEqual(Chem.MolToSmiles(lfrag3), "O=C(O)c1ccc2nc(-c3cc(Cl)cc(Cl)c3)oc2c1")

    smi4 = "CC.O=[Pb]=O"

    lfParams = rdMolStandardize.CleanupParameters()
    lfrag_params = rdMolStandardize.LargestFragmentChooser(lfParams)
    mol4 = Chem.MolFromSmiles(smi4)
    lfrag4 = lfrag_params.choose(mol4)
    self.assertEqual(Chem.MolToSmiles(lfrag4), "CC")

    lfParams = rdMolStandardize.CleanupParameters()
    lfParams.largestFragmentChooserCountHeavyAtomsOnly = True
    lfrag_params = rdMolStandardize.LargestFragmentChooser(lfParams)
    mol4 = Chem.MolFromSmiles(smi4)
    lfrag4 = lfrag_params.choose(mol4)
    self.assertEqual(Chem.MolToSmiles(lfrag4), "O=[Pb]=O")

    lfParams = rdMolStandardize.CleanupParameters()
    lfParams.largestFragmentChooserUseAtomCount = False
    lfrag_params = rdMolStandardize.LargestFragmentChooser(lfParams)
    mol4 = Chem.MolFromSmiles(smi4)
    lfrag4 = lfrag_params.choose(mol4)
    self.assertEqual(Chem.MolToSmiles(lfrag4), "O=[Pb]=O")

    lfParams = rdMolStandardize.CleanupParameters()
    lfParams.largestFragmentChooserCountHeavyAtomsOnly = True
    lfParams.preferOrganic = True
    lfrag_params = rdMolStandardize.LargestFragmentChooser(lfParams)
    mol4 = Chem.MolFromSmiles(smi4)
    lfrag4 = lfrag_params.choose(mol4)
    self.assertEqual(Chem.MolToSmiles(lfrag4), "CC")

    lfParams = rdMolStandardize.CleanupParameters()
    lfParams.largestFragmentChooserUseAtomCount = False
    lfParams.preferOrganic = True
    lfrag_params = rdMolStandardize.LargestFragmentChooser(lfParams)
    mol4 = Chem.MolFromSmiles(smi4)
    lfrag4 = lfrag_params.choose(mol4)
    self.assertEqual(Chem.MolToSmiles(lfrag4), "CC")

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

  def test10NormalizeFromData(self):
    data = """//	Name	SMIRKS
Nitro to N+(O-)=O	[N,P,As,Sb;X3:1](=[O,S,Se,Te:2])=[O,S,Se,Te:3]>>[*+1:1]([*-1:2])=[*:3]
Sulfone to S(=O)(=O)	[S+2:1]([O-:2])([O-:3])>>[S+0:1](=[O-0:2])(=[O-0:3])
Pyridine oxide to n+O-	[n:1]=[O:2]>>[n+:1][O-:2]
// Azide to N=N+=N-	[*,H:1][N:2]=[N:3]#[N:4]>>[*,H:1][N:2]=[N+:3]=[N-:4]
"""
    normalizer1 = rdMolStandardize.Normalizer()
    params = rdMolStandardize.CleanupParameters()
    normalizer2 = rdMolStandardize.NormalizerFromData(data, params)

    imol = Chem.MolFromSmiles("O=N(=O)CCN=N#N", sanitize=False)
    mol1 = normalizer1.normalize(imol)
    mol2 = normalizer2.normalize(imol)
    self.assertEqual(Chem.MolToSmiles(imol), "N#N=NCCN(=O)=O")
    self.assertEqual(Chem.MolToSmiles(mol1), "[N-]=[N+]=NCC[N+](=O)[O-]")
    self.assertEqual(Chem.MolToSmiles(mol2), "N#N=NCC[N+](=O)[O-]")

  def test11FragmentParams(self):
    data = """//   Name	SMARTS
fluorine	[F]
chlorine	[Cl]
        """
    fragremover = rdMolStandardize.FragmentRemoverFromData(data)
    mol = Chem.MolFromSmiles("CN(C)C.Cl.Cl.Br")
    nm = fragremover.remove(mol)
    self.assertEqual(Chem.MolToSmiles(nm), "Br.CN(C)C")

  def test12ChargeParams(self):
    params = """// The default list of AcidBasePairs, sorted from strongest to weakest.
// This list is derived from the Food and Drug: Administration Substance
// Registration System Standard Operating Procedure guide.
//
//	Name	Acid	Base
-SO2H	[!O][SD3](=O)[OH]	[!O][SD3](=O)[O-]
-SO3H	[!O]S(=O)(=O)[OH]	[!O]S(=O)(=O)[O-]
"""
    mol = Chem.MolFromSmiles("C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O")
    # instantiate with default acid base pair library
    reionizer = rdMolStandardize.ReionizerFromData(params, [])
    print("done")
    nm = reionizer.reionize(mol)
    self.assertEqual(Chem.MolToSmiles(nm), "O=S([O-])c1ccc(S(=O)(=O)O)cc1")

  def test13Tautomers(self):
    enumerator = rdMolStandardize.TautomerEnumerator()
    m = Chem.MolFromSmiles("C1(=CCCCC1)O")
    ctaut = enumerator.Canonicalize(m)
    self.assertEqual(Chem.MolToSmiles(ctaut), "O=C1CCCCC1")

    params = rdMolStandardize.CleanupParameters()
    enumerator = rdMolStandardize.TautomerEnumerator(params)
    m = Chem.MolFromSmiles("C1(=CCCCC1)O")
    ctaut = enumerator.Canonicalize(m)
    self.assertEqual(Chem.MolToSmiles(ctaut), "O=C1CCCCC1")

    taut_res = enumerator.Enumerate(m)
    self.assertEqual(len(taut_res), 2)
    ctauts = list(sorted(Chem.MolToSmiles(x) for x in taut_res))
    self.assertEqual(ctauts, ['O=C1CCCCC1', 'OC1=CCCCC1'])
    self.assertEqual(list(taut_res.smiles), ['O=C1CCCCC1', 'OC1=CCCCC1'])
    # this tests the non-templated overload
    self.assertEqual(Chem.MolToSmiles(enumerator.PickCanonical(taut_res)), "O=C1CCCCC1")
    # this tests the templated overload
    self.assertEqual(Chem.MolToSmiles(enumerator.PickCanonical(set(taut_res()))), "O=C1CCCCC1")
    with self.assertRaises(TypeError):
      enumerator.PickCanonical(1)
    with self.assertRaises(TypeError):
      enumerator.PickCanonical([0, 1])
    self.assertEqual(
      Chem.MolToSmiles(
        enumerator.PickCanonical(Chem.MolFromSmiles(x) for x in ['O=C1CCCCC1', 'OC1=CCCCC1'])),
      "O=C1CCCCC1")

    def scorefunc1(mol):
      ' stupid tautomer scoring function '
      p = Chem.MolFromSmarts('[OH]')
      return len(mol.GetSubstructMatches(p))

    def scorefunc2(mol):
      ' stupid tautomer scoring function '
      p = Chem.MolFromSmarts('O=C')
      return len(mol.GetSubstructMatches(p))

    m = Chem.MolFromSmiles("C1(=CCCCC1)O")
    ctaut = enumerator.Canonicalize(m, scorefunc1)
    self.assertEqual(Chem.MolToSmiles(ctaut), "OC1=CCCCC1")
    ctaut = enumerator.Canonicalize(m, scorefunc2)
    self.assertEqual(Chem.MolToSmiles(ctaut), "O=C1CCCCC1")
    # make sure lambdas work
    ctaut = enumerator.Canonicalize(m,
                                    lambda x: len(x.GetSubstructMatches(Chem.MolFromSmarts('C=O'))))
    self.assertEqual(Chem.MolToSmiles(ctaut), "O=C1CCCCC1")

    # make sure we behave if we return something bogus from the scoring function
    with self.assertRaises(TypeError):
      ctaut = enumerator.Canonicalize(m, lambda x: 'fail')

    self.assertEqual(enumerator.ScoreTautomer(Chem.MolFromSmiles('N=c1[nH]cccc1')), 99)
    self.assertEqual(enumerator.ScoreTautomer(Chem.MolFromSmiles('Nc1ncccc1')), 100)

    def scorefunc2(mol):
      ' stupid tautomer scoring function '
      p = Chem.MolFromSmarts('O=C')
      return len(mol.GetSubstructMatches(p))

    m = Chem.MolFromSmiles("C1(=CCCCC1)O")
    ctaut = enumerator.Canonicalize(m, scorefunc1)
    self.assertEqual(Chem.MolToSmiles(ctaut), "OC1=CCCCC1")
    ctaut = enumerator.Canonicalize(m, scorefunc2)
    self.assertEqual(Chem.MolToSmiles(ctaut), "O=C1CCCCC1")
    # make sure lambdas work
    ctaut = enumerator.Canonicalize(m,
                                    lambda x: len(x.GetSubstructMatches(Chem.MolFromSmarts('C=O'))))
    self.assertEqual(Chem.MolToSmiles(ctaut), "O=C1CCCCC1")

    # make sure we behave if we return something bogus from the scoring function
    with self.assertRaises(TypeError):
      ctaut = enumerator.Canonicalize(m, lambda x: 'fail')

    self.assertEqual(enumerator.ScoreTautomer(Chem.MolFromSmiles('N=c1[nH]cccc1')), 99)
    self.assertEqual(enumerator.ScoreTautomer(Chem.MolFromSmiles('Nc1ncccc1')), 100)

    res = enumerator.Enumerate(m)
    # this test the specialized overload
    ctaut = enumerator.PickCanonical(res, scorefunc1)
    self.assertEqual(Chem.MolToSmiles(ctaut), "OC1=CCCCC1")
    ctaut = enumerator.PickCanonical(res, scorefunc2)
    self.assertEqual(Chem.MolToSmiles(ctaut), "O=C1CCCCC1")
    # make sure lambdas work
    ctaut = enumerator.PickCanonical(
      res, lambda x: len(x.GetSubstructMatches(Chem.MolFromSmarts('C=O'))))
    self.assertEqual(Chem.MolToSmiles(ctaut), "O=C1CCCCC1")
    # make sure we behave if we return something bogus from the scoring function
    with self.assertRaises(TypeError):
      ctaut = enumerator.PickCanonical(res, lambda x: 'fail')
    # this test the non-specialized overload
    ctaut = enumerator.PickCanonical(set(res()), scorefunc1)
    self.assertEqual(Chem.MolToSmiles(ctaut), "OC1=CCCCC1")
    ctaut = enumerator.PickCanonical(set(res()), scorefunc2)
    self.assertEqual(Chem.MolToSmiles(ctaut), "O=C1CCCCC1")
    # make sure lambdas work
    ctaut = enumerator.PickCanonical(
      set(res()), lambda x: len(x.GetSubstructMatches(Chem.MolFromSmarts('C=O'))))
    self.assertEqual(Chem.MolToSmiles(ctaut), "O=C1CCCCC1")
    # make sure we behave if we return something bogus from the scoring function
    with self.assertRaises(TypeError):
      ctaut = enumerator.PickCanonical(set(res()), lambda x: 'fail')

  def test14TautomerDetails(self):
    enumerator = rdMolStandardize.TautomerEnumerator()
    m = Chem.MolFromSmiles("c1ccccc1CN=c1[nH]cccc1")
    taut_res = enumerator.Enumerate(m)
    self.assertEqual(len(taut_res.tautomers), 2)
    self.assertEqual(taut_res.modifiedAtoms, (7, 9))
    self.assertEqual(len(taut_res.modifiedBonds), 7)
    self.assertEqual(taut_res.modifiedBonds, (7, 8, 9, 10, 11, 12, 14))

    taut_res = enumerator.Enumerate(m)
    self.assertEqual(len(taut_res.tautomers), 2)
    self.assertEqual(taut_res.modifiedAtoms, (7, 9))

    taut_res = enumerator.Enumerate(m)
    self.assertEqual(len(taut_res.tautomers), 2)
    self.assertEqual(len(taut_res.modifiedBonds), 7)
    self.assertEqual(taut_res.modifiedBonds, (7, 8, 9, 10, 11, 12, 14))

  def test15EnumeratorParams(self):
    # Test a structure with hundreds of tautomers.
    smi68 = "[H][C](CO)(NC(=O)C1=C(O)C(O)=CC=C1)C(O)=O"
    m68 = Chem.MolFromSmiles(smi68)

    enumerator = rdMolStandardize.TautomerEnumerator()
    res68 = enumerator.Enumerate(m68)
    self.assertEqual(len(res68), 72)
    self.assertEqual(len(res68.tautomers), len(res68))
    self.assertEqual(res68.status, rdMolStandardize.TautomerEnumeratorStatus.Completed)

    enumerator = rdMolStandardize.GetV1TautomerEnumerator()
    res68 = enumerator.Enumerate(m68)
    self.assertEqual(len(res68), 292)
    self.assertEqual(len(res68.tautomers), len(res68))
    self.assertEqual(res68.status, rdMolStandardize.TautomerEnumeratorStatus.MaxTransformsReached)

    params = rdMolStandardize.CleanupParameters()
    params.maxTautomers = 50
    enumerator = rdMolStandardize.TautomerEnumerator(params)
    res68 = enumerator.Enumerate(m68)
    self.assertEqual(len(res68), 50)
    self.assertEqual(res68.status, rdMolStandardize.TautomerEnumeratorStatus.MaxTautomersReached)

    sAlaSmi = "C[C@H](N)C(=O)O"
    sAla = Chem.MolFromSmiles(sAlaSmi)
    # test remove (S)-Ala stereochemistry
    self.assertEqual(sAla.GetAtomWithIdx(1).GetChiralTag(), Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
    self.assertEqual(sAla.GetAtomWithIdx(1).GetProp("_CIPCode"), "S")
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = True
    enumerator = rdMolStandardize.TautomerEnumerator(params)
    res = enumerator.Enumerate(sAla)
    for taut in res:
      self.assertEqual(taut.GetAtomWithIdx(1).GetChiralTag(), Chem.ChiralType.CHI_UNSPECIFIED)
      self.assertFalse(taut.GetAtomWithIdx(1).HasProp("_CIPCode"))
    for taut in res.tautomers:
      self.assertEqual(taut.GetAtomWithIdx(1).GetChiralTag(), Chem.ChiralType.CHI_UNSPECIFIED)
      self.assertFalse(taut.GetAtomWithIdx(1).HasProp("_CIPCode"))
    for i, taut in enumerate(res):
      self.assertEqual(Chem.MolToSmiles(taut), Chem.MolToSmiles(res.tautomers[i]))
    self.assertEqual(len(res), len(res.smiles))
    self.assertEqual(len(res), len(res.tautomers))
    self.assertEqual(len(res), len(res()))
    self.assertEqual(len(res), len(res.smilesTautomerMap))
    for i, taut in enumerate(res.tautomers):
      self.assertEqual(Chem.MolToSmiles(taut), Chem.MolToSmiles(res[i]))
      self.assertEqual(Chem.MolToSmiles(taut), res.smiles[i])
      self.assertEqual(Chem.MolToSmiles(taut),
                       Chem.MolToSmiles(res.smilesTautomerMap.values()[i].tautomer))
    for i, k in enumerate(res.smilesTautomerMap.keys()):
      self.assertEqual(k, res.smiles[i])
    for i, v in enumerate(res.smilesTautomerMap.values()):
      self.assertEqual(Chem.MolToSmiles(v.tautomer), Chem.MolToSmiles(res[i]))
    for i, (k, v) in enumerate(res.smilesTautomerMap.items()):
      self.assertEqual(k, res.smiles[i])
      self.assertEqual(Chem.MolToSmiles(v.tautomer), Chem.MolToSmiles(res[i]))
    for i, smiles in enumerate(res.smiles):
      self.assertEqual(smiles, Chem.MolToSmiles(res[i]))
      self.assertEqual(smiles, res.smilesTautomerMap.keys()[i])
    self.assertEqual(Chem.MolToSmiles(res.tautomers[-1]), Chem.MolToSmiles(res[-1]))
    self.assertEqual(Chem.MolToSmiles(res[-1]), Chem.MolToSmiles(res[len(res) - 1]))
    self.assertEqual(Chem.MolToSmiles(res.tautomers[-1]),
                     Chem.MolToSmiles(res.tautomers[len(res) - 1]))
    with self.assertRaises(IndexError):
      res[len(res)]
    with self.assertRaises(IndexError):
      res[-len(res) - 1]
    with self.assertRaises(IndexError):
      res.tautomers[len(res)]
    with self.assertRaises(IndexError):
      res.tautomers[-len(res.tautomers) - 1]

    # test retain (S)-Ala stereochemistry
    self.assertEqual(sAla.GetAtomWithIdx(1).GetChiralTag(), Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
    self.assertEqual(sAla.GetAtomWithIdx(1).GetProp("_CIPCode"), "S")
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    enumerator = rdMolStandardize.TautomerEnumerator(params)
    res = enumerator.Enumerate(sAla)
    for taut in res:
      tautAtom = taut.GetAtomWithIdx(1)
      if (tautAtom.GetHybridization() == Chem.HybridizationType.SP3):
        self.assertEqual(tautAtom.GetChiralTag(), Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        self.assertTrue(tautAtom.HasProp("_CIPCode"))
        self.assertEqual(tautAtom.GetProp("_CIPCode"), "S")
      else:
        self.assertFalse(tautAtom.HasProp("_CIPCode"))
        self.assertEqual(tautAtom.GetChiralTag(), Chem.ChiralType.CHI_UNSPECIFIED)

    eEnolSmi = "C/C=C/O"
    eEnol = Chem.MolFromSmiles(eEnolSmi)
    self.assertEqual(eEnol.GetBondWithIdx(1).GetStereo(), Chem.BondStereo.STEREOE)
    # test remove enol E stereochemistry
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveBondStereo = True
    enumerator = rdMolStandardize.TautomerEnumerator(params)
    res = enumerator.Enumerate(eEnol)
    for taut in res.tautomers:
      self.assertEqual(taut.GetBondWithIdx(1).GetStereo(), Chem.BondStereo.STEREONONE)
    # test retain enol E stereochemistry
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveBondStereo = False
    enumerator = rdMolStandardize.TautomerEnumerator(params)
    res = enumerator.Enumerate(eEnol)
    for taut in res.tautomers:
      if (taut.GetBondWithIdx(1).GetBondType() == Chem.BondType.DOUBLE):
        self.assertEqual(taut.GetBondWithIdx(1).GetStereo(), Chem.BondStereo.STEREOE)

    zEnolSmi = "C/C=C\\O"
    zEnol = Chem.MolFromSmiles(zEnolSmi)
    self.assertEqual(zEnol.GetBondWithIdx(1).GetStereo(), Chem.BondStereo.STEREOZ)
    # test remove enol Z stereochemistry
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveBondStereo = True
    enumerator = rdMolStandardize.TautomerEnumerator(params)
    res = enumerator.Enumerate(zEnol)
    for taut in res:
      self.assertEqual(taut.GetBondWithIdx(1).GetStereo(), Chem.BondStereo.STEREONONE)
    # test retain enol Z stereochemistry
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveBondStereo = False
    enumerator = rdMolStandardize.TautomerEnumerator(params)
    res = enumerator.Enumerate(zEnol)
    for taut in res:
      if (taut.GetBondWithIdx(1).GetBondType() == Chem.BondType.DOUBLE):
        self.assertEqual(taut.GetBondWithIdx(1).GetStereo(), Chem.BondStereo.STEREOZ)

    chembl2024142Smi = "[2H]C1=C(C(=C2C(=C1[2H])C(=O)C(=C(C2=O)C([2H])([2H])[2H])C/C=C(\\C)/CC([2H])([2H])/C=C(/CC/C=C(\\C)/CCC=C(C)C)\\C([2H])([2H])[2H])[2H])[2H]"
    chembl2024142 = Chem.MolFromSmiles(chembl2024142Smi)
    params = Chem.RemoveHsParameters()
    params.removeAndTrackIsotopes = True
    chembl2024142 = Chem.RemoveHs(chembl2024142, params)
    self.assertTrue(chembl2024142.GetAtomWithIdx(12).HasProp("_isotopicHs"))
    # test remove isotopic Hs involved in tautomerism
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveIsotopicHs = True
    enumerator = rdMolStandardize.TautomerEnumerator(params)
    res = enumerator.Enumerate(chembl2024142)
    for taut in res:
      self.assertFalse(taut.GetAtomWithIdx(12).HasProp("_isotopicHs"))
    # test retain isotopic Hs involved in tautomerism
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveIsotopicHs = False
    enumerator = rdMolStandardize.TautomerEnumerator(params)
    res = enumerator.Enumerate(chembl2024142)
    for taut in res:
      self.assertTrue(taut.GetAtomWithIdx(12).HasProp("_isotopicHs"))

  def test16EnumeratorCallback(self):

    class MyTautomerEnumeratorCallback(rdMolStandardize.TautomerEnumeratorCallback):

      def __init__(self, parent, timeout_ms):
        super().__init__()
        self._parent = parent
        self._timeout = timedelta(milliseconds=timeout_ms)
        self._start_time = datetime.now()

      def __call__(self, mol, res):
        self._parent.assertTrue(isinstance(mol, Chem.Mol))
        self._parent.assertTrue(isinstance(res, rdMolStandardize.TautomerEnumeratorResult))
        return (datetime.now() - self._start_time < self._timeout)

    class MyBrokenCallback(rdMolStandardize.TautomerEnumeratorCallback):
      pass

    class MyBrokenCallback2(rdMolStandardize.TautomerEnumeratorCallback):
      __call__ = 1

    # Test a structure with hundreds of tautomers.
    smi68 = "[H][C](CO)(NC(=O)C1=C(O)C(O)=CC=C1)C(O)=O"
    m68 = Chem.MolFromSmiles(smi68)

    def createV1Enumerator():
      enumerator = rdMolStandardize.GetV1TautomerEnumerator()
      enumerator.SetMaxTransforms(10000)
      enumerator.SetMaxTautomers(10000)
      return enumerator

    params = rdMolStandardize.CleanupParameters()
    params.maxTransforms = 10000
    params.maxTautomers = 10000

    enumerator = createV1Enumerator()
    enumerator.SetCallback(MyTautomerEnumeratorCallback(self, 50.0))
    res68 = enumerator.Enumerate(m68)
    # either the enumeration was canceled due to timeout
    # or it has completed very quickly
    hasReachedTimeout = (len(res68.tautomers) < 375
                         and res68.status == rdMolStandardize.TautomerEnumeratorStatus.Canceled)
    hasCompleted = (len(res68.tautomers) == 375
                    and res68.status == rdMolStandardize.TautomerEnumeratorStatus.Completed)
    if hasReachedTimeout:
      print("Enumeration was canceled due to timeout (50 ms)", file=sys.stderr)
    if hasCompleted:
      print("Enumeration has completed", file=sys.stderr)
    self.assertTrue(hasReachedTimeout or hasCompleted)
    self.assertTrue(hasReachedTimeout ^ hasCompleted)

    enumerator = createV1Enumerator()
    enumerator.SetCallback(MyTautomerEnumeratorCallback(self, 10000.0))
    res68 = enumerator.Enumerate(m68)
    # either the enumeration completed
    # or it ran very slowly and was canceled due to timeout
    hasReachedTimeout = (len(res68.tautomers) < 375
                         and res68.status == rdMolStandardize.TautomerEnumeratorStatus.Canceled)
    hasCompleted = (len(res68.tautomers) == 375
                    and res68.status == rdMolStandardize.TautomerEnumeratorStatus.Completed)
    if hasReachedTimeout:
      print("Enumeration was canceled due to timeout (10 s)", file=sys.stderr)
    if hasCompleted:
      print("Enumeration has completed", file=sys.stderr)
    self.assertTrue(hasReachedTimeout or hasCompleted)
    self.assertTrue(hasReachedTimeout ^ hasCompleted)

    enumerator = rdMolStandardize.TautomerEnumerator(params)
    with self.assertRaises(AttributeError):
      enumerator.SetCallback(MyBrokenCallback())
    with self.assertRaises(AttributeError):
      enumerator.SetCallback(MyBrokenCallback2())

    # GitHub #4736
    enumerator = createV1Enumerator()
    enumerator.SetCallback(MyTautomerEnumeratorCallback(self, 50.0))
    enumerator_copy = rdMolStandardize.TautomerEnumerator(enumerator)
    res68 = enumerator.Enumerate(m68)
    res68_copy = enumerator_copy.Enumerate(m68)
    self.assertTrue(res68.status == res68_copy.status)

  def test17PickCanonicalCIPChangeOnChiralCenter(self):

    def get_canonical_taut(res):
      best_idx = max([(rdMolStandardize.TautomerEnumerator.ScoreTautomer(t), i)
                      for i, t in enumerate(res.tautomers)])[1]
      return res.tautomers[best_idx]

    smi = "CC\\C=C(/O)[C@@H](C)C(C)=O"
    mol = Chem.MolFromSmiles(smi)
    self.assertIsNotNone(mol)
    self.assertEqual(mol.GetAtomWithIdx(5).GetChiralTag(), Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    self.assertEqual(mol.GetAtomWithIdx(5).GetProp("_CIPCode"), "R")

    # here the chirality disappears as the chiral center is itself involved in tautomerism
    te = rdMolStandardize.TautomerEnumerator()
    can_taut = te.Canonicalize(mol)
    self.assertIsNotNone(can_taut)
    self.assertEqual(can_taut.GetAtomWithIdx(5).GetChiralTag(), Chem.ChiralType.CHI_UNSPECIFIED)
    self.assertFalse(can_taut.GetAtomWithIdx(5).HasProp("_CIPCode"))
    self.assertEqual(Chem.MolToSmiles(can_taut), "CCCC(=O)C(C)C(C)=O")

    # here the chirality stays even if the chiral center is itself involved in tautomerism
    # because of the tautomerRemoveSp3Stereo parameter being set to false
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    te = rdMolStandardize.TautomerEnumerator(params)
    can_taut = te.Canonicalize(mol)
    self.assertIsNotNone(can_taut)
    self.assertEqual(can_taut.GetAtomWithIdx(5).GetChiralTag(), Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    self.assertEqual(can_taut.GetAtomWithIdx(5).GetProp("_CIPCode"), "S")
    self.assertEqual(Chem.MolToSmiles(can_taut), "CCCC(=O)[C@@H](C)C(C)=O")

    # here the chirality disappears as the chiral center is itself involved in tautomerism
    # the reassignStereo setting has no influence
    te = rdMolStandardize.TautomerEnumerator()
    res = te.Enumerate(mol)
    self.assertEqual(res.status, rdMolStandardize.TautomerEnumeratorStatus.Completed)
    self.assertEqual(len(res.tautomers), 8)
    best_taut = get_canonical_taut(res)
    self.assertIsNotNone(best_taut)
    self.assertEqual(best_taut.GetAtomWithIdx(5).GetChiralTag(), Chem.ChiralType.CHI_UNSPECIFIED)
    self.assertFalse(best_taut.GetAtomWithIdx(5).HasProp("_CIPCode"))
    self.assertEqual(Chem.MolToSmiles(best_taut), "CCCC(=O)C(C)C(C)=O")

    # here the chirality disappears as the chiral center is itself involved in tautomerism
    # the reassignStereo setting has no influence
    params = rdMolStandardize.CleanupParameters()
    params.tautomerReassignStereo = False
    te = rdMolStandardize.TautomerEnumerator(params)
    res = te.Enumerate(mol)
    self.assertEqual(res.status, rdMolStandardize.TautomerEnumeratorStatus.Completed)
    self.assertEqual(len(res.tautomers), 8)
    best_taut = get_canonical_taut(res)
    self.assertIsNotNone(best_taut)
    self.assertEqual(best_taut.GetAtomWithIdx(5).GetChiralTag(), Chem.ChiralType.CHI_UNSPECIFIED)
    self.assertFalse(best_taut.GetAtomWithIdx(5).HasProp("_CIPCode"))
    self.assertEqual(Chem.MolToSmiles(best_taut), "CCCC(=O)C(C)C(C)=O")

    # here the chirality stays even if the chiral center is itself involved in tautomerism
    # because of the tautomerRemoveSp3Stereo parameter being set to false
    # as reassignStereo by default is true, the CIP code has  been recomputed
    # and therefore it is now S (correct)
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    te = rdMolStandardize.TautomerEnumerator(params)
    res = te.Enumerate(mol)
    self.assertEqual(res.status, rdMolStandardize.TautomerEnumeratorStatus.Completed)
    self.assertEqual(len(res.tautomers), 8)
    best_taut = get_canonical_taut(res)
    self.assertIsNotNone(best_taut)
    self.assertEqual(best_taut.GetAtomWithIdx(5).GetChiralTag(), Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    self.assertEqual(best_taut.GetAtomWithIdx(5).GetProp("_CIPCode"), "S")
    self.assertEqual(Chem.MolToSmiles(best_taut), "CCCC(=O)[C@@H](C)C(C)=O")

    # here the chirality stays even if the chiral center is itself involved in tautomerism
    # because of the tautomerRemoveSp3Stereo parameter being set to false
    # as reassignStereo is false, the CIP code has not been recomputed
    # and therefore it is still R (incorrect)
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    params.tautomerReassignStereo = False
    te = rdMolStandardize.TautomerEnumerator(params)
    res = te.Enumerate(mol)
    self.assertEqual(res.status, rdMolStandardize.TautomerEnumeratorStatus.Completed)
    self.assertEqual(len(res.tautomers), 8)
    best_taut = get_canonical_taut(res)
    self.assertIsNotNone(best_taut)
    self.assertEqual(best_taut.GetAtomWithIdx(5).GetChiralTag(), Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    self.assertEqual(best_taut.GetAtomWithIdx(5).GetProp("_CIPCode"), "R")
    self.assertEqual(Chem.MolToSmiles(best_taut), "CCCC(=O)[C@@H](C)C(C)=O")

    smi = "CC\\C=C(/O)[C@@](CC)(C)C(C)=O"
    mol = Chem.MolFromSmiles(smi)
    self.assertIsNotNone(mol)
    self.assertEqual(mol.GetAtomWithIdx(5).GetProp("_CIPCode"), "S")
    self.assertEqual(mol.GetAtomWithIdx(5).GetChiralTag(), Chem.ChiralType.CHI_TETRAHEDRAL_CW)

    # here the chirality stays no matter how tautomerRemoveSp3Stereo
    # is set as the chiral center is not involved in tautomerism
    te = rdMolStandardize.TautomerEnumerator()
    can_taut = te.Canonicalize(mol)
    self.assertIsNotNone(can_taut)
    self.assertEqual(can_taut.GetAtomWithIdx(5).GetChiralTag(), Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    self.assertEqual(can_taut.GetAtomWithIdx(5).GetProp("_CIPCode"), "R")
    self.assertEqual(Chem.MolToSmiles(can_taut), "CCCC(=O)[C@](C)(CC)C(C)=O")

    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    te = rdMolStandardize.TautomerEnumerator(params)
    can_taut = te.Canonicalize(mol)
    self.assertIsNotNone(can_taut)
    self.assertEqual(can_taut.GetAtomWithIdx(5).GetChiralTag(), Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    self.assertEqual(can_taut.GetAtomWithIdx(5).GetProp("_CIPCode"), "R")
    self.assertEqual(Chem.MolToSmiles(can_taut), "CCCC(=O)[C@](C)(CC)C(C)=O")

    # as reassignStereo by default is true, the CIP code has been recomputed
    # and therefore it is now R (correct)
    te = rdMolStandardize.TautomerEnumerator()
    res = te.Enumerate(mol)
    self.assertEqual(res.status, rdMolStandardize.TautomerEnumeratorStatus.Completed)
    self.assertEqual(len(res.tautomers), 4)
    best_taut = get_canonical_taut(res)
    self.assertIsNotNone(best_taut)
    self.assertEqual(best_taut.GetAtomWithIdx(5).GetChiralTag(), Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    self.assertEqual(best_taut.GetAtomWithIdx(5).GetProp("_CIPCode"), "R")
    self.assertEqual(Chem.MolToSmiles(best_taut), "CCCC(=O)[C@](C)(CC)C(C)=O")

    # as reassignStereo is false, the CIP code has not been recomputed
    # and therefore it is still S (incorrect)
    params = rdMolStandardize.CleanupParameters()
    params.tautomerReassignStereo = False
    te = rdMolStandardize.TautomerEnumerator(params)
    res = te.Enumerate(mol)
    self.assertEqual(res.status, rdMolStandardize.TautomerEnumeratorStatus.Completed)
    self.assertEqual(len(res.tautomers), 4)
    best_taut = get_canonical_taut(res)
    self.assertIsNotNone(best_taut)
    self.assertEqual(best_taut.GetAtomWithIdx(5).GetChiralTag(), Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    self.assertEqual(best_taut.GetAtomWithIdx(5).GetProp("_CIPCode"), "S")
    self.assertEqual(Chem.MolToSmiles(best_taut), "CCCC(=O)[C@](C)(CC)C(C)=O")

    # as reassignStereo by default is true, the CIP code has  been recomputed
    # and therefore it is now R (correct)
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    te = rdMolStandardize.TautomerEnumerator(params)
    res = te.Enumerate(mol)
    self.assertEqual(res.status, rdMolStandardize.TautomerEnumeratorStatus.Completed)
    self.assertEqual(len(res.tautomers), 4)
    best_taut = get_canonical_taut(res)
    self.assertIsNotNone(best_taut)
    self.assertEqual(best_taut.GetAtomWithIdx(5).GetChiralTag(), Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    self.assertEqual(best_taut.GetAtomWithIdx(5).GetProp("_CIPCode"), "R")
    self.assertEqual(Chem.MolToSmiles(best_taut), "CCCC(=O)[C@](C)(CC)C(C)=O")

    # here the chirality stays even if the tautomerRemoveSp3Stereo parameter
    # is set to false as the chiral center is not involved in tautomerism
    # as reassignStereo is false, the CIP code has not been recomputed
    # and therefore it is still S (incorrect)
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    params.tautomerReassignStereo = False
    te = rdMolStandardize.TautomerEnumerator(params)
    res = te.Enumerate(mol)
    self.assertEqual(res.status, rdMolStandardize.TautomerEnumeratorStatus.Completed)
    self.assertEqual(len(res.tautomers), 4)
    best_taut = get_canonical_taut(res)
    self.assertIsNotNone(best_taut)
    self.assertEqual(best_taut.GetAtomWithIdx(5).GetChiralTag(), Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    self.assertEqual(best_taut.GetAtomWithIdx(5).GetProp("_CIPCode"), "S")
    self.assertEqual(Chem.MolToSmiles(best_taut), "CCCC(=O)[C@](C)(CC)C(C)=O")

  def test18TautomerEnumeratorResultIter(self):
    smi = "Cc1nnc(NC(=O)N2CCN(Cc3ccc(F)cc3)C(=O)C2)s1"
    mol = Chem.MolFromSmiles(smi)
    self.assertIsNotNone(mol)
    te = rdMolStandardize.TautomerEnumerator()
    res = te.Enumerate(mol)
    res_it = iter(res)
    i = 0
    while 1:
      try:
        t = next(res_it)
      except StopIteration:
        break
      self.assertEqual(Chem.MolToSmiles(t), Chem.MolToSmiles(res[i]))
      i += 1
    self.assertEqual(i, len(res))
    res_it = iter(res)
    i = -len(res)
    while 1:
      try:
        t = next(res_it)
      except StopIteration:
        break
      self.assertEqual(Chem.MolToSmiles(t), Chem.MolToSmiles(res[i]))
      i += 1
    self.assertEqual(i, 0)

  def test19NormalizeFromParams(self):
    params = rdMolStandardize.CleanupParameters()
    params.normalizationsFile = "ThisFileDoesNotExist.txt"
    with self.assertRaises(OSError):
      rdMolStandardize.NormalizerFromParams(params)

  def test20NoneHandling(self):
    with self.assertRaises(ValueError):
      rdMolStandardize.ChargeParent(None)
    with self.assertRaises(ValueError):
      rdMolStandardize.Cleanup(None)
    with self.assertRaises(ValueError):
      rdMolStandardize.FragmentParent(None)
    with self.assertRaises(ValueError):
      rdMolStandardize.Normalize(None)
    with self.assertRaises(ValueError):
      rdMolStandardize.Reionize(None)

  def test21UpdateFromJSON(self):
    params = rdMolStandardize.CleanupParameters()
    # note: these actual parameters aren't useful... they are for testing
    rdMolStandardize.UpdateParamsFromJSON(
      params, """{
    "normalizationData":[
      {"name":"silly 1","smarts":"[Cl:1]>>[F:1]"},
      {"name":"silly 2","smarts":"[Br:1]>>[F:1]"}
    ],
    "acidbaseData":[
      {"name":"-CO2H","acid":"C(=O)[OH]","base":"C(=O)[O-]"},
      {"name":"phenol","acid":"c[OH]","base":"c[O-]"}
    ],
    "fragmentData":[
      {"name":"hydrogen", "smarts":"[H]"}, 
      {"name":"fluorine", "smarts":"[F]"}, 
      {"name":"chlorine", "smarts":"[Cl]"}
    ],
    "tautomerTransformData":[
      {"name":"1,3 (thio)keto/enol f","smarts":"[CX4!H0]-[C]=[O,S,Se,Te;X1]","bonds":"","charges":""},
      {"name":"1,3 (thio)keto/enol r","smarts":"[O,S,Se,Te;X2!H0]-[C]=[C]"}
    ]}""")

    m = Chem.MolFromSmiles("CCC=O")
    te = rdMolStandardize.TautomerEnumerator(params)
    tauts = [Chem.MolToSmiles(x) for x in te.Enumerate(m)]
    self.assertEqual(tauts, ["CC=CO", "CCC=O"])
    self.assertEqual(Chem.MolToSmiles(rdMolStandardize.CanonicalTautomer(m, params)), "CCC=O")
    # now with defaults
    te = rdMolStandardize.TautomerEnumerator()
    tauts = [Chem.MolToSmiles(x) for x in te.Enumerate(m)]
    self.assertEqual(tauts, ["CC=CO", "CCC=O"])
    self.assertEqual(Chem.MolToSmiles(rdMolStandardize.CanonicalTautomer(m)), "CCC=O")

    m = Chem.MolFromSmiles('ClCCCBr')
    nm = rdMolStandardize.Normalize(m, params)
    self.assertEqual(Chem.MolToSmiles(nm), "FCCCF")
    # now with defaults
    nm = rdMolStandardize.Normalize(m)
    self.assertEqual(Chem.MolToSmiles(nm), "ClCCCBr")

    m = Chem.MolFromSmiles('c1cc([O-])cc(C(=O)O)c1')
    nm = rdMolStandardize.Reionize(m, params)
    self.assertEqual(Chem.MolToSmiles(nm), "O=C([O-])c1cccc(O)c1")
    # now with defaults
    nm = rdMolStandardize.Reionize(m)
    self.assertEqual(Chem.MolToSmiles(nm), "O=C([O-])c1cccc(O)c1")

    m = Chem.MolFromSmiles('C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O')
    nm = rdMolStandardize.Reionize(m, params)
    self.assertEqual(Chem.MolToSmiles(nm), "O=S([O-])c1ccc(S(=O)(=O)O)cc1")
    # now with defaults
    nm = rdMolStandardize.Reionize(m)
    self.assertEqual(Chem.MolToSmiles(nm), "O=S(O)c1ccc(S(=O)(=O)[O-])cc1")

    m = Chem.MolFromSmiles('[F-].[Cl-].[Br-].CC')
    nm = rdMolStandardize.RemoveFragments(m, params)
    self.assertEqual(Chem.MolToSmiles(nm), "CC.[Br-]")
    # now with defaults
    nm = rdMolStandardize.RemoveFragments(m)
    self.assertEqual(Chem.MolToSmiles(nm), "CC")


if __name__ == "__main__":
  unittest.main()
