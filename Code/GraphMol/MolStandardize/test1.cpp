//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MolStandardize.h"
#include "Metal.h"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ROMol.h>

#include <iostream>

using namespace RDKit;
using namespace std;

void testCleanup() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n test cleanup" << std::endl;
  string smi1, smi2, smi3, smi4;
  MolStandardize::CleanupParameters params;

  // Test covalent metal is disconnected during standardize.
  smi1 = "CCC(=O)O[Na]";
  unique_ptr<RWMol> m1(SmilesToMol(smi1));
  unique_ptr<RWMol> res1(MolStandardize::cleanup(*m1, params));
  TEST_ASSERT(MolToSmiles(*res1) == "CCC(=O)[O-].[Na+]");

  // Test metal ion is untouched during standardize.
  smi2 = "CCC(=O)[O-].[Na+]";
  unique_ptr<RWMol> m2(SmilesToMol(smi2));
  unique_ptr<RWMol> res2(MolStandardize::cleanup(*m2, params));
  TEST_ASSERT(MolToSmiles(*res2) == "CCC(=O)[O-].[Na+]");

  // Test Hg is disconnected from O during standardize.
  smi3 = "CCC(=O)O[Hg]";
  unique_ptr<RWMol> m3(SmilesToMol(smi3));
  unique_ptr<RWMol> res3(MolStandardize::cleanup(*m3, params));
  TEST_ASSERT(MolToSmiles(*res3) == "CCC(=O)[O-].[Hg+]")

  // Test dimethylmercury is not disconnected during standardize.
  smi4 = "C[Hg]C";
  unique_ptr<RWMol> m4(SmilesToMol(smi4));
  unique_ptr<RWMol> res4(MolStandardize::cleanup(*m4, params));
  TEST_ASSERT(MolToSmiles(*res4) == "C[Hg]C")
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testStandardizeSm() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n test standardize smiles"
                       << std::endl;

  // check aromaticity
  std::string smi1 = "C1=CC=CC=C1";
  std::string ss1 = MolStandardize::standardizeSmiles(smi1);
  TEST_ASSERT(ss1 == "c1ccccc1");

  // both rings should be aromatic
  std::string smi2 = "C[N]1C=NC2=C1C(=O)N(C)C(=O)N2C";
  std::string ss2 = MolStandardize::standardizeSmiles(smi2);
  TEST_ASSERT(ss2 == "Cn1c(=O)c2c(ncn2C)n(C)c1=O");

  // both rings should be aromatic
  std::string smi3 = "C=Cc1ccc2c(c1)[nH]c(=O)/c/2=C\\c1ccc[nH]1";
  std::string ss3 = MolStandardize::standardizeSmiles(smi3);
  TEST_ASSERT(ss3 == "C=Cc1ccc2c(c1)NC(=O)/C2=C\\c1ccc[nH]1");

  // check stereochemistry is correctly perceived
  std::string smi4 = "Cl\\C=C/Cl";
  std::string ss4 = MolStandardize::standardizeSmiles(smi4);
  TEST_ASSERT(ss4 == "Cl/C=C\\Cl");

  // Break metal-organic covalent bonds
  std::string smi5 = "[Na]OC(=O)c1ccccc1";
  std::string ss5 = MolStandardize::standardizeSmiles(smi5);
  TEST_ASSERT(ss5 == "O=C([O-])c1ccccc1.[Na+]");

  // SMILES parsing error should stop tests
  //	std::string smi6 = "C1CCC1C(=O)O.Na";
  //	std::string ss6 = MolStandardize::standardizeSmiles(smi6);
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testMetalDisconnector() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n test metal disconnector"
                       << std::endl;

  MolStandardize::MetalDisconnector md;

  // testing overloaded function
  string smi1 = "CCC(=O)O[Na]";
  unique_ptr<ROMol> m1(SmilesToMol(smi1));
  TEST_ASSERT(m1);
  unique_ptr<ROMol> nm(md.disconnect(*m1));
  TEST_ASSERT(MolToSmiles(*nm) == "CCC(=O)[O-].[Na+]");

  string smi2 = "[Na]OC(=O)CCC(=O)O[Na]";
  unique_ptr<RWMol> m2(SmilesToMol(smi2));
  TEST_ASSERT(m2);
  md.disconnect(*m2);
  TEST_ASSERT(MolToSmiles(*m2) == "O=C([O-])CCC(=O)[O-].[Na+].[Na+]");

  string smi3 = "c1ccccc1[Mg]Br";
  unique_ptr<RWMol> m3(SmilesToMol(smi3));
  TEST_ASSERT(m3);
  md.disconnect(*m3);
  TEST_ASSERT(MolToSmiles(*m3) == "Br[Mg]c1ccccc1");

  string smi4 = "Br[Mg]c1ccccc1CCC(=O)O[Na]";
  unique_ptr<RWMol> m4(SmilesToMol(smi4));
  TEST_ASSERT(m4);
  md.disconnect(*m4);
  TEST_ASSERT(MolToSmiles(*m4) == "O=C([O-])CCc1ccccc1[Mg]Br.[Na+]");

  // test input own metal_non, metal_nof
  // missing out Na
  unique_ptr<ROMol> metal_nof(
      SmartsToMol("[Li,K,Rb,Cs,Fr,Be,Mg,Ca,Sr,Ba,Ra,Sc,Ti,V,Cr,Mn,Fe,Co,"
                  "Ni,Cu,Zn,Al,Ga,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Hf,Ta,"
                  "W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi]~[N,O,F]"));
  md.setMetalNof(*metal_nof);
  string smi5 = "CCC(=O)O[Na]";
  unique_ptr<ROMol> m5(SmilesToMol(smi5));
  TEST_ASSERT(m5);
  unique_ptr<ROMol> nm5(md.disconnect(*m5));
  TEST_ASSERT(MolToSmiles(*nm5) == "CCC(=O)O[Na]");  // not disconnected
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testNormalize() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n test normalize"
                       << std::endl;

  string smi1, smi2, smi3, smi4, smi5, smi6, smi7, smi8, smi9, smi10, smi11,
      smi12, smi13, smi14, smi15, smi16, smi17;
  MolStandardize::CleanupParameters params;

  std::string rdbase = getenv("RDBASE");
  std::string transformFile =
      rdbase + "/Data/MolStandardize/normalizations.txt";
  params.normalizations = transformFile;

  // Normalize nitro group.
  smi1 = "C1(=CC=CC=C1)[N+](=O)[O-]";
  unique_ptr<RWMol> m1(SmilesToMol(smi1));
  unique_ptr<RWMol> res1(MolStandardize::cleanup(*m1, params));
  TEST_ASSERT(MolToSmiles(*res1) == "O=[N+]([O-])c1ccccc1");

  // Normalize nitro group.
  smi2 = "O=[N](=O)c1ccccc1";
  unique_ptr<RWMol> m2(SmilesToMol(smi2));
  unique_ptr<RWMol> res2(MolStandardize::cleanup(*m2, params));
  TEST_ASSERT(MolToSmiles(*res2) == "O=[N+]([O-])c1ccccc1");

  // Normalize nitro group.
  smi3 = "[O-][N+](=O)c1ccccc1";
  unique_ptr<RWMol> m3(SmilesToMol(smi3));
  unique_ptr<RWMol> res3(MolStandardize::cleanup(*m3, params));
  TEST_ASSERT(MolToSmiles(*res3) == "O=[N+]([O-])c1ccccc1");

  // Normalize nitro group.
  smi4 = "[N](=O)(=O)O";
  unique_ptr<RWMol> m4(SmilesToMol(smi4));
  unique_ptr<RWMol> res4(MolStandardize::cleanup(*m4, params));
  TEST_ASSERT(MolToSmiles(*res4) == "O=[N+]([O-])O");

  // Normalize nitro group.
  smi5 = "O[N+](=O)[O-]";
  unique_ptr<RWMol> m5(SmilesToMol(smi5));
  unique_ptr<RWMol> res5(MolStandardize::cleanup(*m5, params));
  TEST_ASSERT(MolToSmiles(*res5) == "O=[N+]([O-])O");

  // Normalize pyridine oxide.
  smi6 = "C1=[N](C=CC=C1)=O";
  unique_ptr<RWMol> m6(SmilesToMol(smi6));
  unique_ptr<RWMol> res6(MolStandardize::cleanup(*m6, params));
  TEST_ASSERT(MolToSmiles(*res6) == "[O-][n+]1ccccc1");

  // Normalize pyridine oxide.
  smi7 = "O=n1ccccc1";
  unique_ptr<RWMol> m7(SmilesToMol(smi7));
  unique_ptr<RWMol> res7(MolStandardize::cleanup(*m7, params));
  TEST_ASSERT(MolToSmiles(*res7) == "[O-][n+]1ccccc1");

  // normalize sulfone.
  smi8 = "C[S+2]([O-])([O-])C";
  unique_ptr<RWMol> m8(SmilesToMol(smi8));
  unique_ptr<RWMol> res8(MolStandardize::cleanup(*m8, params));
  TEST_ASSERT(MolToSmiles(*res8) == "CS(C)(=O)=O");

  // normalize sulfone.
  smi9 = "C[S+2]([O-])([O-])O";
  unique_ptr<RWMol> m9(SmilesToMol(smi9));
  unique_ptr<RWMol> res9(MolStandardize::cleanup(*m9, params));
  TEST_ASSERT(MolToSmiles(*res9) == "CS(=O)(=O)O");

  // normalize sulfoxide..
  smi10 = "CS(=O)C";
  unique_ptr<RWMol> m10(SmilesToMol(smi10));
  unique_ptr<RWMol> res10(MolStandardize::cleanup(*m10, params));
  TEST_ASSERT(MolToSmiles(*res10) == "C[S+](C)[O-]");

  // normalize sulfoxide.
  smi11 = "COC1=CC2=C(C=C1)[N]C(=N2)[S](=O)CC3=C(C(=C(C=N3)C)OC)C";
  unique_ptr<RWMol> m11(SmilesToMol(smi11));
  unique_ptr<RWMol> res11(MolStandardize::cleanup(*m11, params));
  TEST_ASSERT(MolToSmiles(*res11) ==
              "COc1ccc2c(c1)N=C([S+]([O-])Cc1ncc(C)c(OC)c1C)[N]2");

  // Normalize sulfoxide.
  smi12 = "COc1ccc2c(c1)nc([nH]2)S(=O)Cc1ncc(c(c1C)OC)C";
  unique_ptr<RWMol> m12(SmilesToMol(smi12));
  unique_ptr<RWMol> res12(MolStandardize::cleanup(*m12, params));
  TEST_ASSERT(MolToSmiles(*res12) ==
              "COc1ccc2[nH]c([S+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1");

  // Normalize azide.
  smi13 = "C1(=CC=C(C=C1)N)N=[N]#N";
  unique_ptr<RWMol> m13(SmilesToMol(smi13));
  unique_ptr<RWMol> res13(MolStandardize::cleanup(*m13, params));
  TEST_ASSERT(MolToSmiles(*res13) == "[N-]=[N+]=Nc1ccc(N)cc1");

  // Normalize diazo.
  smi14 = "[N](#N)=C1C(NC(N=C1)=O)=O";
  unique_ptr<RWMol> m14(SmilesToMol(smi14));
  unique_ptr<RWMol> res14(MolStandardize::cleanup(*m14, params));
  TEST_ASSERT(MolToSmiles(*res14) == "[N-]=[N+]=C1C=NC(=O)NC1=O");

  // Normalize phosphate.
  smi15 = "C1=NC=C([N]1)CO[P+]([O-])([O-])[O-]";
  unique_ptr<RWMol> m15(SmilesToMol(smi15));
  unique_ptr<RWMol> res15(MolStandardize::cleanup(*m15, params));
  TEST_ASSERT(MolToSmiles(*res15) == "O=P([O-])([O-])OCC1=CN=C[N]1");

  // Normalize hydrazine-diazonium.
  smi16 = "CNNC[N+]#N";
  unique_ptr<RWMol> m16(SmilesToMol(smi16));
  unique_ptr<RWMol> res16(MolStandardize::cleanup(*m16, params));
  TEST_ASSERT(MolToSmiles(*res16) == "CN=[NH+]CN=N");

  // Normalize amidinium.
  smi17 = "[C+](C)(N)N";
  unique_ptr<RWMol> m17(SmilesToMol(smi17));
  unique_ptr<RWMol> res17(MolStandardize::cleanup(*m17, params));
  TEST_ASSERT(MolToSmiles(*res17) == "CC(N)=[NH2+]");
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testNormalizeMultiFrags() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n test normalize multiple frags" << std::endl;

  string smi1, smi2, smi3, smi4, smi5, smi6, smi7, smi8, smi9, smi10, smi11,
      smi12, smi13, smi14, smi15, smi16, smi17;
  MolStandardize::CleanupParameters params;

  std::string rdbase = getenv("RDBASE");
  std::string transformFile =
      rdbase + "/Data/MolStandardize/normalizations.txt";
  params.normalizations = transformFile;

  // All fragments should stay if one gets transformed by normalization.
  smi1 = "[Na]OC(=O)c1ccc(C[S+2]([O-])([O-]))cc1";
  unique_ptr<RWMol> m1(SmilesToMol(smi1));
  unique_ptr<RWMol> res1(MolStandardize::cleanup(*m1, params));
  TEST_ASSERT(MolToSmiles(*res1) == "O=C([O-])c1ccc(C[S](=O)=O)cc1.[Na+]");

  // All fragments should stay if one gets transformed by normalization.
  smi2 = "[Na+].[O-]C(=O)c1ccc(C[S+2]([O-])([O-]))cc1";
  unique_ptr<RWMol> m2(SmilesToMol(smi2));
  unique_ptr<RWMol> res2(MolStandardize::cleanup(*m2, params));
  TEST_ASSERT(MolToSmiles(*res2) == "O=C([O-])c1ccc(C[S](=O)=O)cc1.[Na+]");

  // Recombine non-aromatic 1,3-separated charges.
  smi3 = "C[N-]C(C)=[N+](C)C";
  unique_ptr<RWMol> m3(SmilesToMol(smi3));
  unique_ptr<RWMol> res3(MolStandardize::cleanup(*m3, params));
  TEST_ASSERT(MolToSmiles(*res3) == "CN=C(C)N(C)C");

  // Recombine aromatic 1,3-separated charges.
  smi4 = "[n-]1c(=[N+](C)C)cccc1";
  unique_ptr<RWMol> m4(SmilesToMol(smi4));
  unique_ptr<RWMol> res4(MolStandardize::cleanup(*m4, params));
  TEST_ASSERT(MolToSmiles(*res4) == "CN(C)c1ccccn1");

  // Recombine aromatic 1,3-separated charges.
  smi5 = "C[n+]1c([N-](C))cccc1";
  unique_ptr<RWMol> m5(SmilesToMol(smi5));
  unique_ptr<RWMol> res5(MolStandardize::cleanup(*m5, params));
  TEST_ASSERT(MolToSmiles(*res5) == "CN=c1ccccn1C");

  // Recombine aromatic 1,3-separated charges to form pyrimidone.
  smi6 = "[O-]c1[n+](C)cccc1";
  unique_ptr<RWMol> m6(SmilesToMol(smi6));
  unique_ptr<RWMol> res6(MolStandardize::cleanup(*m6, params));
  TEST_ASSERT(MolToSmiles(*res6) == "Cn1ccccc1=O");

  // Recombine aromatic 1,3-separated charges to form pyrimidone.
  smi7 = "COc1cc2ccc3c4c(OC)cc(OC)c(OC)c4c([O-])[n+](C)c3c2cc1OC";
  unique_ptr<RWMol> m7(SmilesToMol(smi7));
  unique_ptr<RWMol> res7(MolStandardize::cleanup(*m7, params));
  TEST_ASSERT(MolToSmiles(*res7) ==
              "COc1cc2ccc3c4c(OC)cc(OC)c(OC)c4c(=O)n(C)c3c2cc1OC");

  // Recombine non-aromatic 1,5-separated charges.
  smi8 = "C[N-]C=CC=[N+](C)C";
  unique_ptr<RWMol> m8(SmilesToMol(smi8));
  unique_ptr<RWMol> res8(MolStandardize::cleanup(*m8, params));
  TEST_ASSERT(MolToSmiles(*res8) == "CN=CC=CN(C)C");

  // Recombine aromatic 1,5-separated charges.
  smi9 = "[n-]1ccc(=[N+](C)C)cc1";
  unique_ptr<RWMol> m9(SmilesToMol(smi9));
  unique_ptr<RWMol> res9(MolStandardize::cleanup(*m9, params));
  TEST_ASSERT(MolToSmiles(*res9) == "CN(C)c1ccncc1");

  // Recombine aromatic 1,5-separated charges.
  smi10 = "C[n+]1ccc([N-]C)cc1";
  unique_ptr<RWMol> m10(SmilesToMol(smi10));
  unique_ptr<RWMol> res10(MolStandardize::cleanup(*m10, params));
  TEST_ASSERT(MolToSmiles(*res10) == "CN=c1ccn(C)cc1");

  // Shift positive charge from nonprotonated to protonated atom.
  smi11 = "CNC=[N+](C)C";
  unique_ptr<RWMol> m11(SmilesToMol(smi11));
  unique_ptr<RWMol> res11(MolStandardize::cleanup(*m11, params));
  TEST_ASSERT(MolToSmiles(*res11) == "C[NH+]=CN(C)C");

  // Shift positive charge from nonprotonated to protonated atom."
  smi12 = "CNC=CC=[N+](C)C";
  unique_ptr<RWMol> m12(SmilesToMol(smi12));
  unique_ptr<RWMol> res12(MolStandardize::cleanup(*m12, params));
  TEST_ASSERT(MolToSmiles(*res12) == "C[NH+]=CC=CN(C)C");

  // Shift positive charge from nonprotonated to protonated atom."
  smi13 = "[nH]1ccc(=[N+](C)C)cc1";
  unique_ptr<RWMol> m13(SmilesToMol(smi13));
  unique_ptr<RWMol> res13(MolStandardize::cleanup(*m13, params));
  TEST_ASSERT(MolToSmiles(*res13) == "CN(C)c1cc[nH+]cc1");

  // Ensure no transforms inadvertently breaks open rings.
  smi14 = "[O-]C1=CC=CC2=CC=CC=[N+]12";
  unique_ptr<RWMol> m14(SmilesToMol(smi14));
  unique_ptr<RWMol> res14(MolStandardize::cleanup(*m14, params));
  TEST_ASSERT(MolToSmiles(*res14) == "O=c1cccc2ccccn12");

  // Shift positive charge from nonprotonated to protonated atom.
  smi15 = "[nH]1c(=[N+](C)C)cccc1";
  unique_ptr<RWMol> m15(SmilesToMol(smi15));
  unique_ptr<RWMol> res15(MolStandardize::cleanup(*m15, params));
  TEST_ASSERT(MolToSmiles(*res15) == "CN(C)c1cccc[nH+]1");
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testCharge() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n test charges" << std::endl;

  std::string smi1, smi2, smi3, smi4;

  // Reionization should not infinitely loop forever on these molecules.
  smi1 = "CCCCCCCCCCCCCCCCCC(=O)CC(=C)C(=O)O[Ti](=O)(OC(C)C)C(C)C";
  std::string ss = MolStandardize::standardizeSmiles(smi1);
  TEST_ASSERT(
      ss ==
      "C=C(CC(=O)[CH-]CCCCCCCCCCCCCCCC)C(=O)[O-].CC(C)[O-].CCC.[O-2].[Ti+5]");

  // Reionization should not infinitely loop forever on these molecules.
  smi2 =
      "OP(=O)(O)[O-].OP(=O)([O-])[O-].[O-]S(=O)(=O)[O-].[Na+].[Na+].[Na+].[Mg+"
      "2].[Cl-].[Cl-].[K+].[K+]";
  std::string ss2 = MolStandardize::standardizeSmiles(smi2);
  TEST_ASSERT(ss2 ==
              "O=P([O-])(O)O.O=P([O-])([O-])O.O=S(=O)([O-])[O-].[Cl-].[Cl-].[K+"
              "].[K+].[Mg+2].[Na+].[Na+].[Na+]");

  // Charge parent!!
  // Test reionizer moves proton to weaker acid.
  smi3 = "[Na].[Na]";
  std::string ss3 = MolStandardize::standardizeSmiles(smi3);
  TEST_ASSERT(ss3 == "[Na+].[Na+]");

  // TODO: Arguably should become selenite ion... O=[Se]([O-])[O-].
  // Need an AcidBasePair?
  smi4 = "[Na].[Na].O[Se](O)=O";
  std::string ss4 = MolStandardize::standardizeSmiles(smi4);
  TEST_ASSERT(ss4 == "O=[Se](O)O.[Na+].[Na+]");
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testEnumerateTautomerSmiles() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n test enumerate tautomer smiles"
      << std::endl;
  MolStandardize::CleanupParameters params;
  std::string smi1 = "c1(ccccc1)/C=C(/O)\\C";
  std::vector<std::string> tsmiles =
      MolStandardize::enumerateTautomerSmiles(smi1, params);
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

int main() {
  testCleanup();
  testStandardizeSm();
  testMetalDisconnector();
  testNormalize();
  testNormalizeMultiFrags();
  testCharge();
  //	testEnumerateTautomerSmiles();
  return 0;
}
