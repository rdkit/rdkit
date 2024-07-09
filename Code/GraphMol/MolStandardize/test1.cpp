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
#include "Normalize.h"
#include "Fragment.h"
#include "Charge.h"
#include "Metal.h"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ROMol.h>

#include <iostream>

using namespace RDKit;

void testCleanup() {
  BOOST_LOG(rdDebugLog) << "-----------------------\n test cleanup"
                        << std::endl;
  MolStandardize::CleanupParameters params;

  // Test covalent metal is disconnected during standardize.
  {
    RWMOL_SPTR m = "CCC(=O)O[Na]"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "CCC(=O)[O-].[Na+]");
  }
  {
    // Github 5997
    auto m = "CC(=O)O[Mg]OC(=O)C"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "CC(=O)[O-].CC(=O)[O-].[Mg+2]");
  }

  // Test metal ion is untouched during standardize.
  {
    RWMOL_SPTR m = "CCC(=O)[O-].[Na+]"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "CCC(=O)[O-].[Na+]");
  }

  // Test Hg is disconnected from O during standardize.
  {
    RWMOL_SPTR m = "CCC(=O)O[Hg]"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "CCC(=O)[O-].[Hg+]")
  }

  // Test dimethylmercury is not disconnected during standardize.
  {
    RWMOL_SPTR m = "C[Hg]C"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "C[Hg]C")
  }
  BOOST_LOG(rdDebugLog) << "Finished" << std::endl;
}

void testStandardizeSm() {
  BOOST_LOG(rdDebugLog) << "-----------------------\n test standardize smiles"
                        << std::endl;

  // check aromaticity
  {
    std::string smi = "C1=CC=CC=C1";
    std::string ss = MolStandardize::standardizeSmiles(smi);
    TEST_ASSERT(ss == "c1ccccc1");
  }

  // both rings should be aromatic
  {
    std::string smi = "C[N]1C=NC2=C1C(=O)N(C)C(=O)N2C";
    std::string ss = MolStandardize::standardizeSmiles(smi);
    TEST_ASSERT(ss == "Cn1c(=O)c2c(ncn2C)n(C)c1=O");
  }

  // both rings should be aromatic
  {
    std::string smi = "C=Cc1ccc2c(c1)[nH]c(=O)/c/2=C\\c1ccc[nH]1";
    std::string ss = MolStandardize::standardizeSmiles(smi);
    TEST_ASSERT(ss == "C=Cc1ccc2c(c1)NC(=O)/C2=C\\c1ccc[nH]1");
  }

  // check stereochemistry is correctly perceived
  {
    std::string smi = "Cl\\C=C/Cl";
    std::string ss = MolStandardize::standardizeSmiles(smi);
    TEST_ASSERT(ss == "Cl/C=C\\Cl");
  }

  // Break metal-organic covalent bonds
  {
    std::string smi = "[Na]OC(=O)c1ccccc1";
    std::string ss = MolStandardize::standardizeSmiles(smi);
    TEST_ASSERT(ss == "O=C([O-])c1ccccc1.[Na+]");
  }

  // SMILES parsing error should stop tests
  //	std::string smi ="C1CCC1C(=O)O.Na";
  //	std::string ss6 = MolStandardize::standardizeSmiles(smi);
  BOOST_LOG(rdDebugLog) << "Finished" << std::endl;
}

void testMetalDisconnector() {
  BOOST_LOG(rdDebugLog) << "-----------------------\n test metal disconnector"
                        << std::endl;

  MolStandardize::MetalDisconnector md;
  unsigned int failedOp;
  {
    RWMOL_SPTR m(SmilesToMol("[O-]C(=O)C.[Mg+2][O-]C(=O)C", 0, false));
    MolOps::sanitizeMol(*m, failedOp, MolOps::SANITIZE_CLEANUP);
    TEST_ASSERT(m);
    md.disconnect(*m);
    TEST_ASSERT(MolToSmiles(*m) == "CC(=O)[O-].CC(=O)[O-].[Mg+2]");
  }

  // testing overloaded function
  {
    ROMOL_SPTR m("CCC(=O)O[Na]"_smiles);
    TEST_ASSERT(m);
    ROMOL_SPTR nm(md.disconnect(*m));
    TEST_ASSERT(MolToSmiles(*nm) == "CCC(=O)[O-].[Na+]");
  }

  {
    RWMOL_SPTR m("[Na]OC(=O)CCC(=O)O[Na]"_smiles);
    TEST_ASSERT(m);
    md.disconnect(*m);
    TEST_ASSERT(MolToSmiles(*m) == "O=C([O-])CCC(=O)[O-].[Na+].[Na+]");
  }

  {
    RWMOL_SPTR m("c1ccccc1[Mg]Br"_smiles);
    TEST_ASSERT(m);
    md.disconnect(*m);
    TEST_ASSERT(MolToSmiles(*m) == "Br[Mg]c1ccccc1");
  }

  {
    RWMOL_SPTR m("C1(CCCCC1)[Zn]Br"_smiles);
    TEST_ASSERT(m);
    md.disconnect(*m);
    TEST_ASSERT(MolToSmiles(*m) == "[Br-].[CH-]1CCCCC1.[Zn+2]");
  }

  {
    RWMOL_SPTR m("Br[Mg]c1ccccc1CCC(=O)O[Na]"_smiles);
    TEST_ASSERT(m);
    md.disconnect(*m);
    TEST_ASSERT(MolToSmiles(*m) == "O=C([O-])CCc1ccccc1[Mg]Br.[Na+]");
  }

  // test input own dp_metal_non, dp_metal_nof
  // missing out Na
  {
    MolStandardize::MetalDisconnector md2;
    ROMOL_SPTR metal_nof(
        SmartsToMol("[Li,K,Rb,Cs,Fr,Be,Mg,Ca,Sr,Ba,Ra,Sc,Ti,V,Cr,Mn,Fe,Co,"
                    "Ni,Cu,Zn,Al,Ga,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Hf,Ta,"
                    "W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi]~[#7,#8,F]"));
    md2.setMetalNof(*metal_nof);
    ROMOL_SPTR m("CCC(=O)O[Na]"_smiles);
    TEST_ASSERT(m);
    ROMOL_SPTR nm(md2.disconnect(*m));
    TEST_ASSERT(MolToSmiles(*nm) == "CCC(=O)O[Na]");  // not disconnected
  }

  // test that metals are not assigned excess positive charge
  {
    RWMOL_SPTR m(SmilesToMol("[Be](F)(F)(F)OC", 0, false));
    MolOps::sanitizeMol(*m, failedOp, MolOps::SANITIZE_CLEANUP);
    TEST_ASSERT(m);
    md.disconnect(*m);
    TEST_ASSERT(MolToSmiles(*m) == "C[O-].[Be+2].[F-].[F-].[F-]");
  }

  // test that badly written complexes with 4ary nitrogen
  // are not assigned a negative charge nor the metal is
  // stolen electrons from
  {
    RWMOL_SPTR m(SmilesToMol("[Ru](SC)(SC)(SC)N(C)(C)C", 0, false));
    MolOps::sanitizeMol(*m, failedOp, MolOps::SANITIZE_CLEANUP);
    TEST_ASSERT(m);
    md.disconnect(*m);
    TEST_ASSERT(MolToSmiles(*m) == "CN(C)C.C[S-].C[S-].C[S-].[Ru+3]");
  }

  // test that badly written salts are not assigned excess formal charges
  {
    RWMOL_SPTR m(SmilesToMol("[Na+][O-]C(=O)C", 0, false));
    MolOps::sanitizeMol(*m, failedOp, MolOps::SANITIZE_CLEANUP);
    TEST_ASSERT(m);
    md.disconnect(*m);
    TEST_ASSERT(MolToSmiles(*m) == "CC(=O)[O-].[Na+]");
  }
  // make sure that -1 as a valence is dealt with (Github 5997)
  {
    RWMOL_SPTR m(SmilesToMol("[O-]C(=O)C.[Mg+2][O-]C(=O)C", 0, false));
    MolOps::sanitizeMol(*m, failedOp, MolOps::SANITIZE_CLEANUP);
    TEST_ASSERT(m);
    md.disconnect(*m);
    TEST_ASSERT(MolToSmiles(*m) == "CC(=O)[O-].CC(=O)[O-].[Mg+2]");
  }

  // test that badly specified dative bonds are perceived as such
  {
    RWMOL_SPTR m(
        SmilesToMol("CC1=CC=CC2=[N]1[Cu+2]3[N](=C2)NC(=[S]3)N(C)C", 0, false));
    MolOps::sanitizeMol(*m, failedOp, MolOps::SANITIZE_CLEANUP);
    TEST_ASSERT(m);
    md.disconnect(*m);
    TEST_ASSERT(MolToSmiles(*m) == "CC1=CC=CC(C=NNC(=S)N(C)C)=N1.[Cu+2]");
  }

  // test that carbonyl complexes are not assigned excess formal charges
  {
    RWMOL_SPTR m(SmilesToMol("[Ni+2]([C]=O)([C]=O)([C]=O)([C]=O)([C]=O)[C]=O",
                             0, false));
    MolOps::sanitizeMol(
        *m, failedOp, MolOps::SANITIZE_CLEANUP | MolOps::SANITIZE_FINDRADICALS);
    TEST_ASSERT(m);
    md.disconnect(*m);
    TEST_ASSERT(MolToSmiles(*m) ==
                "[C]=O.[C]=O.[C]=O.[C]=O.[C]=O.[C]=O.[Ni+2]");
  }

  // test that dative bonds are handled appropriately
  {
    RWMOL_SPTR m("O->[Fe](<-O)(O)O"_smiles);
    TEST_ASSERT(m);
    md.disconnect(*m);
    TEST_ASSERT(MolToSmiles(*m) == "O.O.[Fe+2].[OH-].[OH-]");
  }

  // test that dative bonds are handled appropriately
  {
    RWMOL_SPTR m("[OH-]->[Co+3](<-[OH-])(<-O)<-O"_smiles);
    TEST_ASSERT(m);
    md.disconnect(*m);
    TEST_ASSERT(MolToSmiles(*m) == "O.O.[Co+3].[OH-].[OH-]");
  }

  // test that pre-existing formal charges on metals are honored (github #3625)
  {
    RWMOL_SPTR m("[Pd+]Cl"_smiles);
    TEST_ASSERT(m);
    md.disconnect(*m);
    TEST_ASSERT(MolToSmiles(*m) == "[Cl-].[Pd+2]");
  }

  {
    RWMOL_SPTR m("[Pd+2]<-[Cl-]"_smiles);
    TEST_ASSERT(m);
    md.disconnect(*m);
    TEST_ASSERT(MolToSmiles(*m) == "[Cl-].[Pd+2]");
  }

  {
    RWMOL_SPTR m("[Al](Cl)(Cl)Cl"_smiles);
    TEST_ASSERT(m);
    md.disconnect(*m);
    TEST_ASSERT(MolToSmiles(*m) == "[Al+3].[Cl-].[Cl-].[Cl-]");
  }

  {
    RWMOL_SPTR m("[Al+](<-[Cl-])(Cl)Cl"_smiles);
    TEST_ASSERT(m);
    md.disconnect(*m);
    TEST_ASSERT(MolToSmiles(*m) == "[Al+3].[Cl-].[Cl-].[Cl-]");
  }

  {
    RWMOL_SPTR m("[Al+](Cl)Cl"_smiles);
    TEST_ASSERT(m);
    md.disconnect(*m);
    TEST_ASSERT(MolToSmiles(*m) == "[Al+3].[Cl-].[Cl-]");
  }

  {
    RWMOL_SPTR m("[Al+2]Cl"_smiles);
    TEST_ASSERT(m);
    md.disconnect(*m);
    TEST_ASSERT(MolToSmiles(*m) == "[Al+3].[Cl-]");
  }

  BOOST_LOG(rdDebugLog) << "Finished" << std::endl;
}

void testNormalize() {
  BOOST_LOG(rdDebugLog) << "-----------------------\n test normalize"
                        << std::endl;

  MolStandardize::CleanupParameters params;

  std::string rdbase = getenv("RDBASE");
  std::string transformFile =
      rdbase + "/Code/GraphMol/MolStandardize/test_data/normalizations.txt";
  params.normalizations = transformFile;

  // Normalize nitro group.
  {
    RWMOL_SPTR m = "C1(=CC=CC=C1)[N+](=O)[O-]"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "O=[N+]([O-])c1ccccc1");
  }

  // Normalize nitro group.
  {
    RWMOL_SPTR m = "O=[N](=O)c1ccccc1"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "O=[N+]([O-])c1ccccc1");
  }

  // Normalize nitro group.
  {
    RWMOL_SPTR m = "[O-][N+](=O)c1ccccc1"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "O=[N+]([O-])c1ccccc1");
  }

  // Normalize nitro group.
  {
    RWMOL_SPTR m = "[N](=O)(=O)O"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "O=[N+]([O-])O");
  }

  // Normalize nitro group.
  {
    RWMOL_SPTR m = "O[N+](=O)[O-]"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "O=[N+]([O-])O");
  }

  // Normalize pyridine oxide.
  {
    RWMOL_SPTR m = "C1=[N](C=CC=C1)=O"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "[O-][n+]1ccccc1");
  }

  // Normalize pyridine oxide.
  {
    RWMOL_SPTR m = "O=n1ccccc1"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "[O-][n+]1ccccc1");
  }

  // normalize sulfone.
  {
    RWMOL_SPTR m = "C[S+2]([O-])([O-])C"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "CS(C)(=O)=O");
  }

  // normalize sulfone.
  {
    RWMOL_SPTR m = "C[S+2]([O-])([O-])O"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "CS(=O)(=O)O");
  }

  // normalize sulfoxide..
  {
    RWMOL_SPTR m = "CS(=O)C"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "C[S+](C)[O-]");
  }

  // normalize sulfoxide.
  {
    RWMOL_SPTR m =
        "COC1=CC2=C(C=C1)[N]C(=N2)[S](=O)CC3=C(C(=C(C=N3)C)OC)C"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) ==
                "COc1ccc2c(c1)N=C([S+]([O-])Cc1ncc(C)c(OC)c1C)[N]2");
  }

  // Normalize sulfoxide.
  {
    RWMOL_SPTR m = "COc1ccc2c(c1)nc([nH]2)S(=O)Cc1ncc(c(c1C)OC)C"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) ==
                "COc1ccc2[nH]c([S+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1");
  }

  // Normalize azide.
  {
    RWMOL_SPTR m = "C1(=CC=C(C=C1)N)N=[N]#N"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "[N-]=[N+]=Nc1ccc(N)cc1");
  }

  // Normalize diazo.
  {
    RWMOL_SPTR m = "[N](#N)=C1C(NC(N=C1)=O)=O"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "[N-]=[N+]=C1C=NC(=O)NC1=O");
  }

  // Normalize phosphate.
  {
    RWMOL_SPTR m = "C1=NC=C([N]1)CO[P+]([O-])([O-])[O-]"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "O=P([O-])([O-])OCC1=CN=C[N]1");
  }

  // Normalize hydrazine-diazonium.
  {
    RWMOL_SPTR m = "CNNC[N+]#N"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "CN=[NH+]CN=N");
  }

  // Normalize amidinium.
  {
    RWMOL_SPTR m = "[C+](C)(N)N"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "CC(N)=[NH2+]");
    BOOST_LOG(rdDebugLog) << "Finished" << std::endl;
  }
}

void testNormalizeMultiFrags() {
  BOOST_LOG(rdDebugLog)
      << "-----------------------\n test normalize multiple frags" << std::endl;

  MolStandardize::CleanupParameters params;

  std::string rdbase = getenv("RDBASE");
  std::string transformFile =
      rdbase + "/Code/GraphMol/MolStandardize/test_data/normalizations.txt";
  params.normalizations = transformFile;

  // All fragments should stay if one gets transformed by normalization.
  {
    RWMOL_SPTR m = "[Na]OC(=O)c1ccc(C[S+2]([O-])([O-]))cc1"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "O=C([O-])c1ccc(C[S](=O)=O)cc1.[Na+]");
  }

  // All fragments should stay if one gets transformed by normalization.
  {
    RWMOL_SPTR m = "[Na+].[O-]C(=O)c1ccc(C[S+2]([O-])([O-]))cc1"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "O=C([O-])c1ccc(C[S](=O)=O)cc1.[Na+]");
  }

  // Recombine non-aromatic 1,3-separated charges.
  {
    RWMOL_SPTR m = "C[N-]C(C)=[N+](C)C"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "CN=C(C)N(C)C");
  }

  // Recombine aromatic 1,3-separated charges.
  {
    RWMOL_SPTR m = "[n-]1c(=[N+](C)C)cccc1"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "CN(C)c1ccccn1");
  }

  // Recombine aromatic 1,3-separated charges.
  {
    RWMOL_SPTR m = "C[n+]1c([N-](C))cccc1"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "CN=c1ccccn1C");
  }

  // Recombine aromatic 1,3-separated charges to form pyrimidone.
  {
    RWMOL_SPTR m = "[O-]c1[n+](C)cccc1"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "Cn1ccccc1=O");
  }

  // Recombine aromatic 1,3-separated charges to form pyrimidone.
  {
    RWMOL_SPTR m =
        "COc1cc2ccc3c4c(OC)cc(OC)c(OC)c4c([O-])[n+](C)c3c2cc1OC"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) ==
                "COc1cc2ccc3c4c(OC)cc(OC)c(OC)c4c(=O)n(C)c3c2cc1OC");
  }

  // Recombine non-aromatic 1,5-separated charges.
  {
    RWMOL_SPTR m = "C[N-]C=CC=[N+](C)C"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "CN=CC=CN(C)C");
  }

  // Recombine aromatic 1,5-separated charges.
  {
    RWMOL_SPTR m = "[n-]1ccc(=[N+](C)C)cc1"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "CN(C)c1ccncc1");
  }

  // Recombine aromatic 1,5-separated charges.
  {
    RWMOL_SPTR m = "C[n+]1ccc([N-]C)cc1"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "CN=c1ccn(C)cc1");
  }

  // Shift positive charge from nonprotonated to protonated atom.
  {
    RWMOL_SPTR m = "CNC=[N+](C)C"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "C[NH+]=CN(C)C");
  }

  // Shift positive charge from nonprotonated to protonated atom."
  {
    RWMOL_SPTR m = "CNC=CC=[N+](C)C"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "C[NH+]=CC=CN(C)C");
  }

  // Shift positive charge from nonprotonated to protonated atom."
  {
    RWMOL_SPTR m = "[nH]1ccc(=[N+](C)C)cc1"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "CN(C)c1cc[nH+]cc1");
  }

  // Ensure no transforms inadvertently breaks open rings.
  {
    RWMOL_SPTR m = "[O-]C1=CC=CC2=CC=CC=[N+]12"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "O=c1cccc2ccccn12");
  }

  // Shift positive charge from nonprotonated to protonated atom.
  {
    RWMOL_SPTR m = "[nH]1c(=[N+](C)C)cccc1"_smiles;
    RWMOL_SPTR res(MolStandardize::cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "CN(C)c1cccc[nH+]1");
    BOOST_LOG(rdDebugLog) << "Finished" << std::endl;
  }
}

void testCharge() {
  BOOST_LOG(rdDebugLog) << "-----------------------\n test charges"
                        << std::endl;

  // Reionization should not infinitely loop forever on these molecules.
  {
    std::string smi = "CCCCCCCCCCCCCCCCCC(=O)CC(=C)C(=O)O[V](=O)(OC(C)C)C(C)C";
    std::string ss = MolStandardize::standardizeSmiles(smi);
    TEST_ASSERT(
        ss ==
        "C=C(CC(=O)[CH-]CCCCCCCCCCCCCCCC)C(=O)[O-].CC(C)[O-].CCC.[O-2].[V+5]");
  }

  // Reionization should not infinitely loop forever on these molecules.
  {
    std::string smi =
        "OP(=O)(O)[O-].OP(=O)([O-])[O-].[O-]S(=O)(=O)[O-].[Na+].[Na+].[Na+].["
        "Mg+"
        "2].[Cl-].[Cl-].[K+].[K+]";
    std::string ss = MolStandardize::standardizeSmiles(smi);
    TEST_ASSERT(
        ss ==
        "O=P([O-])(O)O.O=P([O-])([O-])O.O=S(=O)([O-])[O-].[Cl-].[Cl-].[K+"
        "].[K+].[Mg+2].[Na+].[Na+].[Na+]");
  }

  // Charge parent!!
  // Test reionizer moves proton to weaker acid.
  {
    std::string smi = "[Na].[Na]";
    std::string ss = MolStandardize::standardizeSmiles(smi);
    TEST_ASSERT(ss == "[Na+].[Na+]");
  }

  // TODO: Arguably should become selenite ion... O=[Se]([O-])[O-].
  // Need an AcidBasePair?
  {
    std::string smi = "[Na].[Na].O[Se](O)=O";
    std::string ss = MolStandardize::standardizeSmiles(smi);
    TEST_ASSERT(ss == "O=[Se](O)O.[Na+].[Na+]");
    BOOST_LOG(rdDebugLog) << "Finished" << std::endl;
  }

  // Test that tetrazolate salts are not neutralised
  {
    std::string smi = "c1nn[n-]n1.[Na+]";
    std::string ss = MolStandardize::standardizeSmiles(smi);
    TEST_ASSERT(ss == "[Na+].c1nn[n-]n1");
    BOOST_LOG(rdDebugLog) << "Finished" << std::endl;
  }

  // Test that tetrazolate zwitterions are not neutralised
  {
    std::string smi = "CCCc1cc(-c2ccccc2)cc(CCC)[n+]1-c1nn[n-]n1";
    std::string ss = MolStandardize::standardizeSmiles(smi);
    TEST_ASSERT(ss == "CCCc1cc(-c2ccccc2)cc(CCC)[n+]1-c1nn[n-]n1");
    BOOST_LOG(rdDebugLog) << "Finished" << std::endl;
  }
}

void testEnumerateTautomerSmiles() {
  BOOST_LOG(rdDebugLog)
      << "-----------------------\n test enumerate tautomer smiles"
      << std::endl;
  MolStandardize::CleanupParameters params;
  std::string smi = "c1(ccccc1)/C=C(/O)\\C";
  std::vector<std::string> tsmiles =
      MolStandardize::enumerateTautomerSmiles(smi, params);
  BOOST_LOG(rdDebugLog) << "Finished" << std::endl;
}

void testMetalDisconnectorLigandExpo() {
  BOOST_LOG(rdDebugLog) << "-----------------------\n test metal disconnector "
                           "on LigandExpo ligands"
                        << std::endl;

  std::list<std::pair<std::string, std::string>> ligandExpoSmiles{
      {"[Be](O[P@@](=O)(O)O[P@](=O)(O)OCCCNc1ccc(cc1[N+](=O)[O-])[N+](=O)[O-])("
       "F)(F)F",
       "O=[N+]([O-])c1ccc(NCCCO[P@@](=O)(O)OP(=O)(O)O)c([N+](=O)[O-])c1"},
      {"[Be](OP(=O)(O)OP(=O)(O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)n2cnc3c2ncnc3N)O)"
       "O)(F)(F)F",
       "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O"},
      {"[Be-](O[P@@](=O)(O)O[P@](=O)(O)OCCNc1cccc(c1)[N+](=O)[O-])(F)(F)F",
       "O=[N+]([O-])c1cccc(NCCO[P@@](=O)(O)OP(=O)(O)O)c1"},
      {"[Be-](O[P@@](=O)(O)O[P@](=O)(O)OCCNc1ccccc1[N+](=O)[O-])(F)(F)F",
       "O=[N+]([O-])c1ccccc1NCCO[P@@](=O)(O)OP(=O)(O)O"},
      {"[Be](O[P@@](=O)(O)O[P@](=O)(O)OCCNc1ccc(cc1)[N+](=O)[O-])(F)(F)F",
       "O=[N+]([O-])c1ccc(NCCO[P@@](=O)(O)OP(=O)(O)O)cc1"},
      {"[Be-](O[P@@](=O)(O)O[P@](=O)(O)OCCNc1ccc(cc1[N+](=O)[O-])[N+](=O)[O-])("
       "F)(F)F",
       "O=[N+]([O-])c1ccc(NCCO[P@@](=O)(O)OP(=O)(O)O)c([N+](=O)[O-])c1"},
      {"[Be-](O[P@@](=O)(O)O[P@](=O)(O)OCC[N@@](C)c1ccccc1[N+](=O)[O-])(F)(F)F",
       "CN(CCO[P@@](=O)(O)OP(=O)(O)O)c1ccccc1[N+](=O)[O-]"},
      {"c1c2c3c4c(c1NC(=O)CI)C=CC=[N]4[Ru+2]56([N]3=CC=C2)([N]7=C(C=CC=C7)C8=["
       "N]5C=CC=C8)[N]9=C(C=CC=C9)C1=CC=CC=[N]61",
       "O=C(CI)Nc1cc2cccnc2c2ncccc12"},
      {"c1c2c(cc3c1C4=[N]([Ru]356([S]7CC[S]5CC[S]6CC7)C#[O+])C=CC=C4)C(=O)NC2="
       "O",
       "O=C1NC(=O)c2cc(-c3ccccn3)c([Ru])cc21"},
      {"c1c2c(cc3c1C4=[N]([Ru]356([S]7CC[S]5CC[S]6CC7)N=C=S)C=CC=C4)C(=O)NC2=O",
       "O=C1NC(=O)c2cc(-c3ccccn3)c([Ru])cc21"},
      {"C1=C2C(=C(C3=[N]2[Fe]45n6c1c(c(c6C=C7[N]4=C(C=C8N5C(=C3)C(=C8CCC(=O)O)"
       "CC(=O)O)C(=C7CC(=O)O)CCC(=O)O)CC(=O)O)CCC(=O)O)CCC(=O)O)CC(=O)O",
       "O=C(O)CCC1=C(CC(=O)O)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(CC(=O)O)c5CCC(="
       "O)O)C(CCC(=O)O)=C4CC(=O)O)c(CC(=O)O)c3CCC(=O)O"},
      {"C1=CC23=C4([Ru+]2567892%10([CH]3=[CH]5[CH]6=[CH]74)[CH]3=[CH]8C9[CH]2=["
       "CH]%103)C=C1",
       "c1ccc2ccccc2c1"},
      {"c1cc2c3c4c1C=CC=[N]4[Ru]5678(C3=CNC2=O)(C9=C5[C]6C7=C89)C#O",
       "O=c1[nH]ccc2c1ccc1cccnc12"},
      {"c1cc2c3c4c1N=CC=[N]4[Ru+2]56([N]3=CC=N2)([N]7=C8C(=CC=C7)c9c(nc1cc(c("
       "cc1n9)F)F)C1=CC=C[N]5=C18)[N]1=CC=Nc2c1c1c(cc2)N=CC=[N]61",
       "Fc1cc2nc3c4cccnc4c4ncccc4c3nc2cc1F"},
      {"c1cc2c(cc1c3ccc4c(c3)nc5c(n4)C6=CC=C[N]7=C6C8=[N]([Ru]791([N]3=CC="
       "Cc4c3c3c(cc4)C=CC=[N]93)[N]3=CC=Cc4c3c3c(cc4)C=CC=[N]13)C=CC=C58)nc1c("
       "n2)C2=CC=C[N]3=C2C2=[N]([Ru]334([N]5=CC=Cc6c5c5c(cc6)C=CC=[N]35)[N]3="
       "CC=Cc5c3c3c(cc5)C=CC=[N]43)C=CC=C12",
       "c1cnc2c(c1)c1nc3ccc(-c4ccc5nc6c7cccnc7c7ncccc7c6nc5c4)cc3nc1c1cccnc12"},
      {"c1cc2c(cc1Cl)nc3c(n2)C4=C5C6=[N](C=CC=C36)[Ru+2]78([N]5=CC=C4)([N]9=CC="
       "Nc1c9c2c(cc1)N=CC=[N]72)[N]1=CC=Nc2c1c1c(cc2)N=CC=[N]81",
       "Clc1ccc2nc3c4cccnc4c4ncccc4c3nc2c1"},
      {"c1cc2c(cc1C(N)N)NC3=[N]2[Zn+2][N]4=C(C3=O)Nc5c4ccc(c5)C(N)N",
       "NC(N)c1ccc2nc(C(=O)c3nc4ccc(C(N)N)cc4[nH]3)[nH]c2c1"},
      {"c1cc2c(c(c1)[N+](=O)[OH-])nc3c4c5c6c(c3n2)C=CC=[N]6[Ru]78([N]5=CC=C4)(["
       "N]9=C1C(=CC=C9)C=CC2=CC=C[N]7=C21)[N]1=C2C(=CC=C1)C=CC1=CC=C[N]8=C12",
       "O=[N+]([1O-])c1cccc2nc3c4cccnc4c4ncccc4c3nc12"},
      {"c1cc2c(cc1N(=O)=O)nc3c(n2)C4=CC=C[N]5=C4C6=[N]([Ru]578([N]9=CC="
       "Cc1c9c2c(cc1)C=CC=[N]72)[N]1=CC=Cc2c1c1c(cc2)C=CC=[N]81)C=CC=C36",
       "O=[N+]([O-])c1ccc2nc3c4cccnc4c4ncccc4c3nc2c1"},
      {"c1cc2c(cc1O)c3c4c(c5c6c3n2[Ru]7891([N]6=CC=C5)([CH]2=[CH]7C8[CH]9=[CH]"
       "12)C#O)C(=O)NC4=O",
       "O=C1NC(=O)c2c1c1cccnc1c1[nH]c3ccc(O)cc3c21"},
      {"c1cc2c(cc1O)c3c4n2[Ru]56([N]7=C4C(=CC(=C7)F)C8=C3C(=O)NC8=O)([S]9CC[S]"
       "5CC[S]6CC9)[N]CS",
       "O=C1NC(=O)c2c1c1cc(F)cnc1c1[nH]c3ccc(O)cc3c21"},
      {"c1cc2c(cc1OCCN3CCCCC3)O[Cu]45[N](=C2)c6cc(c(cc6[N]4=Cc7ccc(cc7O5)"
       "OCCN8CCCCC8)F)F",
       "Oc1cc(OCCN2CCCCC2)ccc1C=Nc1cc(F)c(F)cc1N=Cc1ccc(OCCN2CCCCC2)cc1O"},
      {"c1cc2c(cc1OCCN3CCCCC3)O[Ni]45[N](=C2)c6cc(c(cc6[N]4=Cc7ccc(cc7O5)"
       "OCCN8CCCCC8)F)F",
       "Oc1cc(OCCN2CCCCC2)ccc1C=Nc1cc(F)c(F)cc1N=Cc1ccc(OCCN2CCCCC2)cc1O"},
      {"c1cc2c(cc(c(c2cc1)N)/N=N/c3ccc(cc3)c4ccc(cc4)/N=N/"
       "c5c(c6c(c(c5)S(=O)(=O)[O-][Na+])cccc6)N)S(=O)(=O)[O-][Na+]",
       "Nc1c(/N=N/c2ccc(-c3ccc(/N=N/"
       "c4cc(S(=O)(=O)O)c5ccccc5c4N)cc3)cc2)cc(S(=O)(=O)O)c2ccccc12"},
      {"c1cc2n3c1C=C4C=CC5=[N]4[Fe]36[N]7=C(C=CC7=C2)C=C8N6C(=C5)C=C8",
       "C1=Cc2cc3ccc(cc4nc(cc5ccc(cc1n2)[nH]5)C=C4)[nH]3"},
      {"c1cc2n3c1C(=C4C=CC5=[N]4[Zn]36N7C(=C5C(F)(F)F)C=CC7=C(C8=[N]6C(=C2C(F)("
       "F)F)C=C8)C(F)(F)F)C(F)(F)F",
       "FC(F)(F)c1c2nc(c(C(F)(F)F)c3ccc([nH]3)c(C(F)(F)F)c3nc(c(C(F)(F)F)"
       "c4ccc1[nH]4)C=C3)C=C2"},
      {"C1=CC2=[N](C=C1)[Pt]34[N]5=C(C=CC=C5)NC6=CC=CC(=[N]63)C7=[N]4C(=CC=C7)"
       "N2",
       "c1ccc(Nc2cccc(-c3cccc(Nc4ccccn4)n3)n2)nc1"},
      {"c1ccc2c(c1)c3c4c(c5c6c3n2[Ru]7([N]6=CC(=C5)N)(NCC8=[N]7C=CC=C8)(C#O)Cl)"
       "C(=O)NC4=O",
       "Nc1cnc2c(c1)c1c(c3c4ccccc4[nH]c23)C(=O)NC1=O"},
      {"c1ccc2c(c1)-c3c4c(c5c6c3[NH]2[Ru]7891([N]6=CC=C5)([CH]2=[CH]7[CH]8=C9=["
       "CH]12)C#O)C(=O)NC4=O",
       "O=C1NC(=O)c2c1c1cccnc1c1[nH]c3ccccc3c21"},
      {"c1ccc2c(c1)c3c4n2[Ru]5678(C9C5C6C7C89)([N]1=C4C(=CC(=C1)F)C1=C3C(=O)"
       "NC1=O)C#[O+]",
       "O=C1NC(=O)c2c1c1cc(F)cnc1c1[nH]c3ccccc3c21"},
      {"c1ccc2c(c1)c3n4c2N=C5c6ccccc6C7=[N]5[Ga]48N9C(=NC1=[N]8C(=N3)c2c1cccc2)"
       "c1ccccc1C9=N7",
       "c1ccc2c(c1)-c1nc-2nc2[nH]c(nc3nc(nc4[nH]c(n1)c1ccccc41)-c1ccccc1-3)"
       "c1ccccc21"},
      {"c1ccc2c(c1)C(=CC3=[N]2[Ru]4567([C]8[C]4[C]5[C]6[C]78)(OC3=O)[O])N9CCN("
       "CC9)C(=O)CCCC[C@H]1[C@@H]2[C@H](CS1)NC(=O)N2",
       "O=C1N[C@@H]2[C@H](CCCCC(=O)N3CCN(c4cc(C(=O)O)nc5ccccc45)CC3)SC[C@@H]"
       "2N1"},
      {"c1ccc2c(c1)C=[N]3c4ccc(cc4O[Cr+]3O2)CCC(=O)O",
       "O=C(O)CCc1ccc(N=Cc2ccccc2O)c(O)c1"},
      {"c1ccc2c(c1)nc3c4c5c6c(c3n2)C=CC=[N]6[Ru]78([N]5=CC=C4)([N]9=CC="
       "Cc1c9c2c(cc1)C=CC=[N]72)[N]1=CC=Cc2c1c1c(cc2)C=CC=[N]81",
       "c1ccc2nc3c4cccnc4c4ncccc4c3nc2c1"},
      {"c1ccc2c(c1)nc3c4c5c6c(c3n2)C=CC=[N]6[Ru]78([N]5=CC=C4)([N]9=CC="
       "Cc1c9c2c(cc1)C=CC=[N]72)[N]1=CC=Cc2c1c1c(cc2)C=CC=[N]81",
       "c1ccc2nc3c4cccnc4c4ncccc4c3nc2c1"},
      {"c1ccc2c(c1)nc3c(n2)C4=CC=C[N]5=C4C6=[N]([Ru+2]578([N]9=C1C(=NC=C9)C="
       "CC2=NC=C[N]7=C21)[N]1=C2C(=NC=C1)C=CC1=NC=C[N]8=C12)C=CC=C36",
       "c1ccc2nc3c4cccnc4c4ncccc4c3nc2c1"},
      {"c1ccc2c(c1)nc3c(n2)C4=CC=C[N]5=C4C6=[N]([Ru]578([N]9=C(C=CC=C9)C1=CC="
       "CC=[N]71)[N]1=CC=CC=C1C1=CC=CC=[N]81)C=CC=C36",
       "c1ccc2nc3c4cccnc4c4ncccc4c3nc2c1"},
      {"c1ccc2c(c1)nc3c(n2)C4=CC=C[N]5=C4C6=[N]([Ru]578([N]9=C(C=CC=C9)C1=[N]"
       "7C=CC=C1)[N]1=CC=CC=C1C1=CC=CC=[N]81)C=CC=C36",
       "c1ccc2nc3c4cccnc4c4ncccc4c3nc2c1"},
      {"C1=CC(C=C1)[Ru@]2(N3C=CC=C4C3=C5[N]2=C6C=CC(=CC6=C5C7=C4C(=O)NC7=O)O)C#"
       "O",
       "O=C1NC(=O)c2c1c1ccc[nH]c1c1nc3ccc(O)cc3c21"},
      {"c1c(c(cc2c1nc3c4c5c6c(c3n2)C=CC=[N]6[Ru]78([N]5=CC=C4)([N]9=C1C(=NC=C9)"
       "C=CC2=NC=C[N]7=C21)[N]1=C2C(=NC=C1)C=CC1=NC=C[N]8=C12)C#N)C#N",
       "N#Cc1cc2nc3c4cccnc4c4ncccc4c3nc2cc1C#N"},
      {"c1ccc(cc1)C2=C3C=CC4=CC5=[N]6C(=C(c7ccc8n7[Zn]6(N43)[N]9=C2C=CC9=C8)"
       "c1ccc(cc1)C(=O)O)C=C5",
       "O=C(O)c1ccc(-c2c3nc(cc4ccc([nH]4)c(-c4ccccc4)c4nc(cc5ccc2[nH]5)C=C4)C="
       "C3)cc1"},
      {"c1ccc(cc1)C2=C3C=CC4=[N]3[Fe]56n7c2ccc7C=C8[N]5=C(C=C8)C(=C9N6C(="
       "C4c1ccccc1)C=C9)c1ccccc1",
       "C1=Cc2nc1cc1ccc([nH]1)c(-c1ccccc1)c1nc(c(-c3ccccc3)c3ccc([nH]3)c2-"
       "c2ccccc2)C=C1"},
      {"c1ccc(cc1)C2=C3C=CC4=[N]3[Fe]56n7c2ccc7C=C8[N]5=C(C=C8)C(=C9N6C(=C4)C="
       "C9)c1ccccc1",
       "C1=Cc2nc1cc1ccc([nH]1)c(-c1ccccc1)c1nc(cc3ccc([nH]3)c2-c2ccccc2)C=C1"},
      {"c1ccc(cc1)C2=C3C=CC4=[N]3[Fe]56n7c2ccc7N=C8[N]5=C(C=C8)C(=C9N6C(=N4)C="
       "C9)c1ccccc1",
       "C1=Cc2nc1nc1ccc([nH]1)c(-c1ccccc1)c1nc(nc3ccc([nH]3)c2-c2ccccc2)C=C1"},
      {"c1ccc(cc1)C2=C(C3=[N]4C(=C5C=CC6=C(C(=C7C=CC8=[N]7[Co]4(N56)n9c8ccc92)"
       "c1ccccc1)c1ccccc1)C=C3)c1ccccc1",
       "C1=Cc2nc1c(-c1ccccc1)c(-c1ccccc1)c1ccc([nH]1)c1nc(c(-c3ccccc3)c(-"
       "c3ccccc3)c3ccc2[nH]3)C=C1"},
      {"c1ccc(cc1)C[C@H]2C(=O)O[Cu]3[N]2=Cc4ccccc4O3",
       "O=C(O)[C@H](Cc1ccccc1)N=Cc1ccccc1O"},
      {"c1ccc(cc1)[C@H]2COC3=[N]2[Rh]4c5c3cccc5C6=[N]4[C@H](CO6)c7ccccc7",
       "[Rh]c1c(C2=N[C@@H](c3ccccc3)CO2)cccc1C1=N[C@@H](c2ccccc2)CO1"},
      {"c1cc(ccc1C[N@]23CCCN4[Ni]25N(CCCN5CC3)CC4)C[N@]67CCC[N@@]8[Ni@@]69[N@]("
       "CCC[N@]9CC7)CC8",
       "c1cc(CN2CCCNCCNCCCNCC2)ccc1CN1CCCNCCNCCCNCC1"},
      {"c1cc(ccc1C[N@]23CCC[N@H]4[Cu@@]25[N@H](CCC[N@H]5CC3)CC4)C[N@@]67CCC[N@"
       "H]8[Cu]69[N@H](CCC[N@H]9CC7)CC8",
       "c1cc(CN2CCCNCCNCCCNCC2)ccc1CN1CCCNCCNCCCNCC1"},
      {"c1cc(ccc1CNC(=O)C23[CH]4=[CH]5[Re]426([CH]5=[CH]63)(C#O)(C#O)C#O)S(=O)("
       "=O)N",
       "NS(=O)(=O)c1ccc(CNC(=O)C2C=CC=C2)cc1"},
      {"c1cc(ccc1NC(=O)CCCC[C@H]2[C@@H]3[C@H](CS2)NC(=O)N3)S(=O)(=O)[N-]4CC["
       "NH2][Ir+3]4",
       "NCCNS(=O)(=O)c1ccc(NC(=O)CCCC[C@@H]2SC[C@@H]3NC(=O)N[C@H]23)cc1"},
      {"c1cc(ccc1NC(=O)CCCC[C@H]2[C@@H]3[C@H](CS2)NC(=O)N3)S(=O)(=O)[N@@]4CCN["
       "Ru]456789([CH]1=[CH]5[CH]6=[CH]7[CH]8=[CH]91)Cl",
       "NCCNS(=O)(=O)c1ccc(NC(=O)CCCC[C@@H]2SC[C@@H]3NC(=O)N[C@H]23)cc1"},
      {"C1CCC(CC1)[NH2][Pt+2][NH3]", "NC1CCCCC1"},
      {"c1ccc(cc1)[P](c2ccccc2)(c3ccccc3)[Pt+2]4(c5ccccc5C6=CC=CC=[N]64)Cl",
       "c1ccc(P(c2ccccc2)c2ccccc2)cc1"},
      {"C1CC[C@@H]2[C@@H](C1)[N]34CC5=[N]([Fe]367([N]2(CC8=CC=CC=[N]68)CC(=O)"
       "O7)OC(=O)C4)C=CC=C5",
       "O=C(O)CN(Cc1ccccn1)[C@@H]1CCCC[C@H]1N(CC(=O)O)Cc1ccccn1"},
      {"C1CC[C@@H]2[C@@H](C1)[NH2][Pt+2][NH2]2", "N[C@@H]1CCCC[C@H]1N"},
      {"C1=CC=[N]2C(=C1)C3=CC=CC4=[N]3[Pt+]2([N]5=C4C=CC=C5)Cl",
       "c1ccc(-c2cccc(-c3ccccn3)n2)nc1"},
      {"C1=CC=[N]2C(=C1)C3=CC=CC=[N]3[Ru+2]245([N]6=C7C(=CC=C6)c8c(nc9c1c2c3c("
       "c9n8)C=CC=[N]3[Ru+2]36([N]2=CC=C1)([N]1=C(C=CC=C1)C1=[N]3C=CC=C1)[N]1="
       "C(C=CC=C1)C1=[N]6C=CC=C1)C1=CC=C[N]4=C17)[N]1=CC=CC=C1C1=CC=CC=[N]51",
       "c1cnc2c(c1)c1nc3c4cccnc4c4ncccc4c3nc1c1cccnc12"},
      {"C1=CC=[N]2C(=C1)C3=CC=CC=[N]3[Ru+2]245([N]6=C(C=CC=C6)C7=[N]4C=CC=C7)["
       "N]8=C9C(=CC=C8)c1c(nc2c3c4c6c(c2n1)C=CC=[N]6[Ru+2]12([N]4=CC=C3)([N]3="
       "C(C=CC=C3)C3=[N]1C=CC=C3)[N]1=C(C=CC=C1)C1=[N]2C=CC=C1)C1=C9[N]5=CC=C1",
       "c1cnc2c(c1)c1nc3c4cccnc4c4ncccc4c3nc1c1cccnc12"},
      {"C1=CC=[N]2C(=C1)C3=[N]4[Ru+2]25([N]6=C(C4=CC=C3)C=CC=C6)[N]7=C8C(=CC="
       "C7)C=CC9=CC=C[N]5=C98",
       "c1ccc(-c2cccc(-c3ccccn3)n2)nc1"},
      {"C1=CC=[N]2C(=C1)N(C3=[N]([Ru]2Cl)C=CC=C3)CCCCCN",
       "NCCCCCN(c1ccccn1)c1ccccn1"},
      {"c1ccnc(c1)C2=[N](C=CC=C2)[Ru]([N]3=CNC=C3)[N]4=C(C=CC=C4)c5ccccn5",
       "c1ccc(-c2ccccn2)nc1"},
      {"C1=CC=[N](C=C1)[Pt+2]([NH3])([NH3])Cl", "c1ccncc1"},
      {"C1CC(=O)N(C1=O)CCC23=[CH]4[Fe]2567891([CH]3=[CH]5[C@H]64)[CH]2=[CH]7C8["
       "CH]9=[CH]12",
       "O=C1CCC(=O)N1CCC1=CCC=C1"},
      {"C1CC(=O)N(C1=O)CCCCCN2C3=[N](C=CC=C3)[Ru]([N]4=CC=CC=C42)Cl",
       "O=C1CCC(=O)N1CCCCCN(c1ccccn1)c1ccccn1"},
      {"C1[C@H]2[C@@H]([C@@H](S1)CCCCC(=O)NCCC3O[Co]45(O6[Co]7(O4[Co]8(O3)(O7["
       "Co]6(O58)([N]9=CC=CC=C9)[O])([N]1=CC=CC=C1)[O])([N]1=CC=CC=C1)[O])[O])"
       "NC(=O)N2",
       "O=C(CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12)NCCC(O)O"},
      {"C1[C@H]2[C@@H]([C@@H](S1)CCCCC(=O)NCCC3O[Co]45(O6[Co]7(O4[Co]8(O3)(O7["
       "Co]6(O58)([N]9=CC=CC=C9)([O])[O])([N]1=CC=CC=C1)[O])[N]1=CC=CC=C1)([N]"
       "1=CC=CC=C1)[O])NC(=O)N2",
       "O=C(CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12)NCCC(O)O"},
      {"C1[C@H]2[C@@H]([C@@H](S1)CCCCC(=O)NCCCO[C@]34C5=[N]6C(=CC=C5)C7=[N](["
       "Co]68([N]9=CC=CC=C9C1=CC=CC3=[N]81)([N]1=CC=CC=C41)[O])C=CC=C7)NC(=O)"
       "N2",
       "O=C(CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12)NCCCOC(c1ccccn1)(c1cccc(-"
       "c2ccccn2)n1)c1cccc(-c2ccccn2)n1"},
      {"C1[C@@H]([C@H](O[C@H]1N2C=[N](C3=C2N=C(NC3=O)N)[Pt](N)(N)[N]4=C5C=CC="
       "CC5=C6C=CC=CC6=C4)COP(=O)(O)O)O",
       "Nc1nc2c(ncn2[C@H]2C[C@H](O)[C@@H](COP(=O)(O)O)O2)c(=O)[nH]1"},
      {"C1=C[N]2=C3C(=C1)C=CC4=CC=C[N](=C43)[Re]2(C#O)(C#O)C#O",
       "c1cnc2c(c1)ccc1cccnc12"},
      {"C1C[N]2(CCC3=CC=CC=[N]3[Cu]2[N]4=C1C=CC=C4)CCCNC(=O)CCCC[C@H]5[C@@H]6["
       "C@H](CS5)NC(=O)N6",
       "O=C(CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12)NCCCN(CCc1ccccn1)CCc1ccccn1"},
      {"C1C[N]2(CCC3=CC=CC=[N]3[Cu]2([N]4=C1C=CC=C4)N=[N+]=[NH-])CCNC(=O)CCCC["
       "C@H]5[C@@H]6[C@H](CS5)NC(=O)N6",
       "O=C(CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12)NCCN(CCc1ccccn1)CCc1ccccn1"},
      {"C1C[N]2(CCC3=CC=CC=[N]3[Cu]2([N]4=C1C=CC=C4)(N=[N+]=[NH-])[O])CCCNC(=O)"
       "CCCC[C@H]5[C@@H]6[C@H](CS5)NC(=O)N6",
       "O=C(CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12)NCCCN(CCc1ccccn1)CCc1ccccn1"},
      {"C1C[N]2(CCC3=CC=CC=[N]3[Cu]2([N]4=C1C=CC=C4)[O])CCCNC(=O)CCCC[C@H]5[C@@"
       "H]6[C@H](CS5)NC(=O)N6",
       "O=C(CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12)NCCCN(CCc1ccccn1)CCc1ccccn1"},
      {"C1C[N]2(CCC3=CC=CC=[N]3[Cu]2([N]4=C1C=CC=C4)[O])CCNC(=O)CCCC[C@H]5[C@@"
       "H]6[C@H](CS5)NC(=O)N6",
       "O=C(CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12)NCCN(CCc1ccccn1)CCc1ccccn1"},
      {"C1C[N]2(CCC3=CC=CC=[N]3[Cu]2([N]4=C1C=CC=C4)([O])[O])CCNC(=O)CCCC[C@H]"
       "5[C@@H]6[C@H](CS5)NC(=O)N6",
       "O=C(CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12)NCCN(CCc1ccccn1)CCc1ccccn1"},
      {"C1C[N]2(CCC3=CC=CC=[N]3[Cu]2([N]4=C1C=CC=C4)OO)CCNC(=O)CCCC[C@H]5[C@@H]"
       "6[C@H](CS5)NC(=O)N6",
       "O=C(CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12)NCCN(CCc1ccccn1)CCc1ccccn1"},
      {"C1C[N]2(CCC3=CC=CC=[N]3[Cu]2([N]4=C1C=CC=C4)([O])OO)CCNC(=O)CCCC[C@H]5["
       "C@@H]6[C@H](CS5)NC(=O)N6",
       "O=C(CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12)NCCN(CCc1ccccn1)CCc1ccccn1"},
      {"C1C[N]2(CCC3=CC=CC=[N]3[Cu]2[N]4=CC=CC=C41)CCCCNC(=O)CCCC[C@H]5[C@@H]6["
       "C@H](CS5)NC(=O)N6",
       "O=C(CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12)NCCCCN(CCc1ccccn1)CCc1ccccn1"},
      {"C1C[N]2(CCC3=[N]([Cu]2[N]4=C1C=CC=C4)C=CC=C3)CCNC(=O)CCCC[C@H]5[C@@H]6["
       "C@H](CS5)NC(=O)N6",
       "O=C(CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12)NCCN(CCc1ccccn1)CCc1ccccn1"},
      {"C1=C[N](=CN1)[Cu+]([N]2=CNC=C2)([N]3=CNC=C3)([N]4=CNC=C4)([OH2])[OH2]",
       "c1c[nH]cn1"},
      {"c1cn(cn1)[Os+2]23([N]4=C(C=CC=C4)C5=CC=CC=[N]52)[N]6=C(C=CC=C6)C7=CC="
       "CC=[N]37",
       "c1ccc(-c2ccccn2)nc1"},
      {"c1cn(cn1)[Os+2]23([N]4=C(C=CC=C4)C5=CC=CC=[N]52)[N]6=CC=CC=C6C7=CC=CC=["
       "N]37",
       "c1ccc(-c2ccccn2)nc1"},
      {"C1=C[N](=CN1)[Re+](C#O)(C#O)(C#O)[OH2]", "c1c[nH]cn1"},
      {"c1cn(cn1)[Ru+2]23([N]4=C(C=CC=C4)C5=CC=CC=[N]52)[N]6=CC=CC=C6C7=CC=CC=["
       "N]37",
       "c1ccc(-c2ccccn2)nc1"},
      {"c1cn(cn1)[Ru+2]23([N]4=C(C=CC=C4)C5=[N]2C=CC=C5)[N]6=CC=CC=C6C7=CC=CC=["
       "N]37",
       "c1ccc(-c2ccccn2)nc1"},
      {"C1C[N@@H]2CC[N@@H]3[Cu]24[N@@H](C1)CC[N@@H]4CCC3", "C1CNCCNCCCNCCNC1"},
      {"C1C[NH2][Pt+]2([NH2]1)I[Pt+]3(I2)[NH2]CC[NH2]3", "NCCN"},
      {"C1C[N][Ru]23456([N]1)(C7=C2C38=C4(C5=C67)CC9=C(C8)CC=CC9)Cl",
       "C1=CCC2=C(C1)Cc1ccccc1C2"},
      {"C1C([O@@H][Fe]234[O@H]C(C[C@@]56[OH]2[Fe]7([OH]3[C@]1(CC([O@@H]7)[O])C("
       "[O@@H]4)[O])([O@H]C(C5)[O])[O@H]C6[O])[O])[O]",
       "[O]C(O)CC(O)(CC([O])O)C([O])O"},
      {"C1C[O@@H][Y]234[N]1(CC[O@@H]2)C(C[OH]3)(C[O@@H]4)CO",
       "OCCN(CCO)C(CO)(CO)CO"},
      {"C[As-](C)(O)OC1=[O+][Na]2345OC(=O)C[O+]2c6c7cc(cc6Cc8cc(cc(c8[O+]3CC(="
       "O)O)Cc9cc(cc(c9[O+]4C1)Cc1cc(cc(c1[O+]5CC(=O)O)C7)S(=O)(=O)O)S(=O)(=O)"
       "O)S(=O)(=O)O)S(=O)(=O)O",
       "C[AsH](C)(O)OC(=O)COc1c2cc(S(=O)(=O)O)cc1Cc1cc(S(=O)(=O)O)cc(c1OCC(=O)"
       "O)Cc1cc(S(=O)(=O)O)cc(c1OCC(=O)O)Cc1cc(S(=O)(=O)O)cc(c1OCC(=O)O)C2"},
      {"C[C-]12C3(=C4([Ir+3]1356(C2(=C54C)C)([N-](CC7=[N]6C=CC(=C7)c8ccc(cc8)S("
       "=O)(=O)N)S(=O)(=O)c9ccccc9)[Cl-])C)C",
       "NS(=O)(=O)c1ccc(-c2ccnc(CNS(=O)(=O)c3ccccc3)c2)cc1"},
      {"CC12C3([Ir]1456(C2(C4(C53C)C)C)([N]7=C(C=CC=C7)C(=O)N6CCc8ccc(cc8)S(=O)"
       "(=O)N)Cl)C",
       "NS(=O)(=O)c1ccc(CCNC(=O)c2ccccn2)cc1"},
      {"CC12C3([Ir]1456(C2(C4(C53C)C)C)[N]7=CC=C(C=C7C(=O)[N]6(CCCc8ccc(cc8)S(="
       "O)(=O)N)Cl)O)C",
       "NS(=O)(=O)c1ccc(CCCN(Cl)C(=O)c2cc(O)ccn2)cc1"},
      {"CC12C3([Ir]1456(C2(C4(C53C)C)C)([N]7=CC=C(C=C7C(=O)N6CC(=O)Nc8ccc(cc8)"
       "S(=O)(=O)N)O)Cl)C",
       "NS(=O)(=O)c1ccc(NC(=O)CNC(=O)c2cc(O)ccn2)cc1"},
      {"CC12C3([Ir]1456(C2(C4(C53C)C)C)([N]7=C(C=C(C=C7)O)C(=O)N6CCc8ccc(cc8)S("
       "=O)(=O)N)Cl)C",
       "NS(=O)(=O)c1ccc(CCNC(=O)c2cc(O)ccn2)cc1"},
      {"CC12=C3([Ir]145C2(=C4(C53C)C)C)C", "CC1=C(C)C(C)C(C)=C1C"},
      {"CC12=C3([Rh+]145(C2(=C4(C53C)CCNC(=O)CCCCC6C7C(CS6)NC(=O)N7)C)(Cl)(Cl)"
       "Cl)C",
       "CC1=C(C)C(C)C(CCNC(=O)CCCCC2SCC3NC(=O)NC32)=C1C"},
      {"C[C@@]12CC(=O)N[C@@]13C[C@H]4[C@H]([C@](C5=[N]4[Ni+]67[N]3=C([C@H]2CCC("
       "=O)O)C=C8N6C(=C9C1=[N]7[C@H](C5)[C@H]([C@@H]1C[C@@H](C9=O)SC)CC(=O)O)["
       "C@H]([C@@H]8CC(=O)O)CCC(=O)O)(C)CC(=O)N)CCC(=O)O",
       "CS[C@H]1C[C@@H]2C3=N[C@H](CC4=N[C@@H](C[C@@]56N=C(C=C7NC(=C3C1=O)[C@@H]"
       "(CCC(=O)O)[C@@H]7CC(=O)O)[C@@H](CCC(=O)O)[C@]5(C)CC(=O)N6)[C@@H](CCC(="
       "O)O)[C@]4(C)CC(N)=O)[C@H]2CC(=O)O"},
      {"C[C@@]12CC(=O)N[C@@]13C[C@H]4[C@H]([C@](C5=[N]4[Ni+]67[N]3=C([C@H]2CCC("
       "=O)O)C=C8N6C(=C9C(=O)CC[C@@H]1C9=[N]7[C@H](C5)[C@H]1CC(=O)O)[C@H]([C@@"
       "H]8CC(=O)O)CCC(=O)O)(C)CC(=O)N)CCC(=O)O",
       "C[C@@]1(CC(N)=O)C2=N[C@@H](C[C@@]34N=C(C=C5NC(=C6C(=O)CC[C@@H]7C6=N[C@"
       "H](C2)[C@H]7CC(=O)O)[C@@H](CCC(=O)O)[C@@H]5CC(=O)O)[C@@H](CCC(=O)O)[C@]"
       "3(C)CC(=O)N4)[C@H]1CCC(=O)O"},
      {"Cc1c2cc3[n+]4c(cc5c(c(c6n5[Ga+3]47n2c(c1CCC(=O)O)cc8[n+]7c(c6)C(=C8CCC("
       "=O)O)C)C=C)C)C(=C3C)C=C",
       "C=CC1=C(C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5C=C)C(C)=C4CCC(=O)O)c("
       "CCC(=O)O)c3C"},
      {"Cc1c2cc3[n+]4c(cc5c(c(c6n5[Mg@]47n2c(c1CCC(=O)O)cc8[n+]7c(c6)C(=C8CCC(="
       "O)O)C)C)C=C)C(=C3C=C)C",
       "C=CC1=C(C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5CCC(=O)O)C(CCC(=O)O)="
       "C4C)c(C)c3C=C"},
      {"CC1=C2C=CC3=C4C2=[N](C=C1)[Re]([N]4=CC=C3C)(C#O)(C#O)C#O",
       "Cc1ccnc2c1ccc1c(C)ccnc12"},
      {"CC1=C2[C@]([C@H]([C@H]3[N-]2[Co+2]45[N]6=C1[C@H](C(C6=CC7=[N]4C(=C(C8=["
       "N]5[C@@]3([C@@]([C@@H]8CCC(=O)N)(C)CC(=O)N)C)C)[C@@]([C@@H]7CCC(=O)N)("
       "C)CC(=O)N)(C)C)CCC(=O)N)CC(=O)N)(C)CCC(=O)NC[C@@H](C)O",
       "CC1=C2N=C(C=C3N=C(C(C)=C4N[C@H]([C@H](CC(N)=O)[C@@]4(C)CCC(=O)NC[C@@H]("
       "C)O)[C@]4(C)N=C1[C@@H](CCC(N)=O)[C@]4(C)CC(N)=O)[C@@H](CCC(N)=O)C3(C)C)"
       "[C@@H](CCC(N)=O)[C@]2(C)CC(N)=O"},
      {"Cc1c2n3c(c1C=C)C=C4C(=C(C5=CC6=[N]7[Fe@]3(N45)[N]8=C(C=C7C(=C6CCC(=O)O)"
       "C=O)C(=C(C8=C2)[C@H](CC\\C=C(/C)\\CC\\C=C(/C)\\CC\\C=C(/"
       "C)\\CCC=C(C)C)O)C)CCC(=O)O)C",
       "C=Cc1c(C)c2cc3nc(cc4nc(cc5[nH]c(cc1[nH]2)c(C)c5CCC(=O)O)C(CCC(=O)O)="
       "C4C=O)C(C)=C3[C@@H](O)CC/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CCC=C(C)C"},
      {"Cc1c2n3c(c1C=C)C=C4C(=C(C5=[N]4[Cr]36[N]7=C(C=C8N6C(=C5)C(=C8CCC(=O)O)"
       "C)C(=C(C7=C2)C)CCC(=O)O)C=C)C",
       "C=CC1=C(C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5CCC(=O)O)C(CCC(=O)O)="
       "C4C)c(C)c3C=C"},
      {"Cc1c2n3c(c1C=C)C=C4C(=C(C5=[N]4[Fe@@]36[N]7=C(C=C8N6C(=C5)C(=C8C)CCC(="
       "O)O)C(=C(C7C2)C=C)C)CCC(=O)O)C",
       "C=CC1=C(C)C2=NC1Cc1[nH]c(c(C=C)c1C)C=C1N=C(C=c3[nH]c(c(C)c3CCC(=O)O)="
       "C2)C(CCC(=O)O)=C1C"},
      {"Cc1c2n3c(c1C=C)C=C4C(=C(C5=[N]4[Fe]36[N]7=C(C=C8N6C(=C5)C(=C8C)CCC(=O)"
       "O)C(=C(C7=C2)[C@H](C)OO)C)CCC(=O)O)C",
       "C=Cc1c(C)c2cc3nc(cc4[nH]c(cc5nc(cc1[nH]2)C(C)=C5CCC(=O)O)c(CCC(=O)O)"
       "c4C)C(C)=C3[C@H](C)OO"},
      {"Cc1c2n3c(c1C=C)C=C4C(=C(C5=[N]4[Fe]36[N]7=C(C=C8N6C(=C5)C(=C8C)CCC(=O)"
       "O)C(=C(C7=C2)C=O)C)CCC(=O)O)C",
       "C=Cc1c(C)c2cc3nc(cc4[nH]c(cc5nc(cc1[nH]2)C(C)=C5CCC(=O)O)c(CCC(=O)O)"
       "c4C)C(C)=C3C=O"},
      {"Cc1c2n3c(c1C=C)C=C4C(=C(C5=[N]4[Fe]36([N]7=C(C=C8N6C(=C5)C(=C8CCC(=O)O)"
       "C)C(=C(C7=C2)C)CCC(=O)O)c9ccccc9)C=C)C",
       "C=CC1=C(C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5CCC(=O)O)C(CCC(=O)O)="
       "C4C)c(C)c3C=C"},
      {"Cc1c2n3c(c1C=C)C=C4C(=C(C5=[N]4[Fe]36([N]7=C(C=C8N6C(=C5)C(=C8CCC(=O)O)"
       "C)C(=C(C7=C2)C)CCC(=O)O)c9ccc(cc9)Cl)C=C)C",
       "C=CC1=C(C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5CCC(=O)O)C(CCC(=O)O)="
       "C4C)c(C)c3C=C"},
      {"Cc1c2n3c(c1C=C)C=C4C(=C(C5=[N]4[Fe]36[N]7=C(C=C8N6C(=C5)C(=C8C=O)CCC(="
       "O)O)C(=C(C7=C2)[C@H](CC/C=C(\\C)/CC/C=C(\\C)/CCC=C(C)C)O)C)CCC(=O)O)C",
       "C=Cc1c(C)c2cc3nc(cc4[nH]c(cc5nc(cc1[nH]2)C(C)=C5CCC(=O)O)c(CCC(=O)O)"
       "c4C=O)C(C)=C3[C@@H](O)CC/C=C(\\C)CC/C=C(\\C)CCC=C(C)C"},
      {"Cc1c2n3c(c1C=C)C=C4C(=C(C5=[N]4[Zn]36[N]7=C(C=C8N6C(=C5)C(=C8C)CCC(=O)"
       "O)C(=C(C7=C2)C=C)C)CCC(=O)O)C",
       "C=CC1=C(C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5C=C)C(C)=C4CCC(=O)O)c("
       "CCC(=O)O)c3C"},
      {"Cc1c2n3c(c1CCC(=O)O)C=C4C(=C(C5=CC6=[N]7[Fe]3(N45)[N]8=C(C=C7[C@@](C6="
       "O)(C)CC(=O)O)C(=O)[C@](C8=C2)(C)CC(=O)O)C)CCC(=O)O",
       "Cc1c(CCC(=O)O)c2cc3[nH]c(cc4nc(cc5nc(cc1[nH]2)C(=O)[C@]5(C)CC(=O)O)C(="
       "O)[C@]4(C)CC(=O)O)c(C)c3CCC(=O)O"},
      {"Cc1c2n3c(c1CCC(=O)O)C=C4C(=C(C5=[N]4[Co]36[N]7=C(C=C8N6C(=C5)C(=C8C=C)"
       "C)C(=C(C7=C2)C=C)C)C)CCC(=O)O",
       "C=CC1=C(C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5CCC(=O)O)C(CCC(=O)O)="
       "C4C)c(C)c3C=C"},
      {"Cc1c2n3c(c1CCC(=O)O)C=C4C(=C(C5=[N]4[Fe+2]36[N]7=C(C(C(C7=C2)C=C)C)C(="
       "C8N6C(=C5)C(=C8C=C)C)c9ccccc9)C)CCC(=O)O",
       "C=Cc1c(C)c2cc3nc(cc4[nH]c(cc5nc(c(-c6ccccc6)c1[nH]2)C(C)C5C=C)c(C)"
       "c4CCC(=O)O)C(CCC(=O)O)=C3C"},
      {"Cc1c2n3c(c1CCC(=O)O)C=C4C(=C(C5=[N]4[Fe]36[N]7=C(C=C8N6C(=C5)C(=C8C=C)"
       "C)[C@]9(C(C7=C2)C=CS9)C)C)CCC(=O)O",
       "C=Cc1c(C)c2cc3nc(cc4[nH]c(cc5nc(cc1[nH]2)[C@]1(C)SC=CC51)c(C)c4CCC(=O)"
       "O)C(CCC(=O)O)=C3C"},
      {"Cc1c2n3c(c1CCC(=O)O)C=C4C(=C(C5=[N]4[Fe]36[N]7=C(C=C8N6C(=C5)C(=C8C)C="
       "C)C(=C(C7=C2)C)C=C)C)CCC(=O)O",
       "C=CC1=C(C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5C=C)C(C)=C4CCC(=O)O)c("
       "CCC(=O)O)c3C"},
      {"Cc1c2n3c(c1CCC(=O)O)C=C4C(=C(C5=[N]4[Fe]36[N]7=C(C=C8N6C(=C5)C(=C8C)C="
       "C)C(=C(C7=C2)C)/C=C/[N+](=O)[O-])C)CCC(=O)O",
       "C=Cc1c(C)c2cc3nc(cc4[nH]c(cc5nc(cc1[nH]2)C(C)=C5CCC(=O)O)c(CCC(=O)O)"
       "c4C)C(C)=C3/C=C/[N+](=O)[O-]"},
      {"Cc1c2n3c(c1CCC(=O)O)C=C4C(=C(C5=[N]4[Fe]36[N]7=C(C=C8N6C(=C5)C(=C8C)"
       "CCC(=O)O)C(=C(C7=C2)C)C=C)C)CCC(=O)O",
       "C=CC1=C(C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5CCC(=O)O)C(C)=C4CCC(=O)"
       "O)c(CCC(=O)O)c3C"},
      {"Cc1c2n3c(c1CCC(=O)O)C=C4C(=C(C5=[N]4[In]36[N]7=C(C=C8N6C(=C5)C(=C8C)C="
       "C)C(=C(C7=C2)C)C=C)C)CCC(=O)O",
       "C=CC1=C(C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5C=C)C(C)=C4CCC(=O)O)c("
       "CCC(=O)O)c3C"},
      {"Cc1c2n3c(c1CCC(=O)O)C=C4C(=C(C5=[N]4[Mg]36[N]7=C(C=C8N6C(=C5)C(=C8)C)C("
       "=CC7=C2)C)C)CCC(=O)O",
       "CC1=Cc2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)cc5C)C(C)=C4CCC(=O)O)c(CCC(=O)O)"
       "c3C"},
      {"Cc1c2n3c(c1CCC(=O)O)C=C4C(=C(C5=[N]4[Mn]36[N]7=C(C=C8N6C(=C5)C(=C8C)"
       "CCC(=O)O)C(=C(C7=C2)CCC(=O)O)C)CCC(=O)O)C",
       "CC1=C(CCC(=O)O)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5CCC(=O)O)C(CCC(=O)"
       "O)=C4C)c(CCC(=O)O)c3C"},
      {"Cc1c2n3c(c1CCC(=O)O)C=C4C(=C(C5=[N]4[Ni]36[N]7=C(C=C8N6C(=C5)C(=C8C=C)"
       "C)C(=C(C7=C2)C=C)C)C)CCC(=O)O",
       "C=CC1=C(C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5CCC(=O)O)C(CCC(=O)O)="
       "C4C)c(C)c3C=C"},
      {"Cc1c2n3c(c1CCC(=O)O)C=C4C(=C(C5=[N]4[Zn]36[N]7=C(C=C(C7=C2)C)C=C8N6C(="
       "C5)C=C8C)C)CCC(=O)O",
       "CC1=Cc2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5CCC(=O)O)C(CCC(=O)O)=C4C)"
       "cc3C"},
      {"Cc1c2n3c(c1C=C)[O+]=C4C(=C(C5=[N]4[Fe]36[N]7=C(C=C(C7=C2)C)C=C8N6C(=C5)"
       "C(=C8CCC(=O)O)C)C=C)C",
       "C=CC1=C(C)c2nc1cc1[nH]c(cc3nc(cc4[n-]c([o+]2)c(C=C)c4C)C(C)=C3)c(CCC(="
       "O)O)c1C"},
      {"Cc1c2n3c(c1C(=C)OO)C=C4C(=C(C5=[N]4[Fe]36[N]7=C(C=C8N6C(=C5)C(=C8CCC(="
       "O)O)C)C(=C(C7=C2)C)CCC(=O)O)C=C)C",
       "C=CC1=C(C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5CCC(=O)O)C(CCC(=O)O)="
       "C4C)c(C)c3C(=C)OO"},
      {"Cc1c2n3c(c1C(=O)C)C=C4C(=C(C5=[N]4[Fe]36[N]7=C(C=C8N6C(=C5)C(=C8C)CCC(="
       "O)O)C(=C(C7=C2)C(=C)O)C)CCC(=O)O)C",
       "C=C(O)C1=C(C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5C(C)=O)C(C)=C4CCC(="
       "O)O)c(CCC(=O)O)c3C"},
      {"CC1=C([C@]2([C@@]3(C(=C(C4=[N]3[Co+]56[N]2=C1C=C7N5C(=CC8=[N]6C(=C4)C(="
       "C8CCC(=O)O)C)C(=C7C)CCC(=O)O)C)C)C)C)C",
       "CC1=C(CCC(=O)O)C2=NC1=CC1=N[C@@](C)(C(C)=C1C)[C@@]1(C)N=C(C=c3[nH]c(c("
       "CCC(=O)O)c3C)=C2)C(C)=C1C"},
      {"CC1=C([C@@]2([C@]3(C(=C(C4=[N]3[Co+]56[N]2=C1C=C7N5C(=CC8=[N]6C(=C4)C(="
       "C8CCC(=O)O)C)C(=C7C)CCC(=O)O)C)C)C)C)C",
       "CC1=C(CCC(=O)O)C2=NC1=CC1=N[C@](C)(C(C)=C1C)[C@]1(C)N=C(C=c3[nH]c(c("
       "CCC(=O)O)c3C)=C2)C(C)=C1C"},
      {"Cc1cc2c(cc1C)N3C=[N]2[Co]456([N]7=C8C(=C9[N]4=C(C=C1[N]5=C(C(=C2N6C("
       "C7C(C8(CCC(=O)NCC(OP(=O)(OC4C(OC3C4O)CO)O)C)C)CC(=O)N)(C(C2CCC(=O)N)(C)"
       "CC(=O)N)C)C)C(C1CCC(=O)N)(C)CC(=O)N)C(C9CCC(=O)N)(C)C)C)CC1C(C(C(O1)"
       "n1cnc2c1ncnc2N)O)O",
       "CC1=C2N=C(C=C3N=C(C(C)=C4NC(C)(C5N=C1C(C)(CCC(=O)NCC(C)OP(=O)(O)OC1C("
       "CO)OC(n6cnc7cc(C)c(C)cc76)C1O)C5CC(N)=O)C(C)(CC(N)=O)C4CCC(N)=O)C(C)("
       "CC(N)=O)C3CCC(N)=O)C(C)(C)C2CCC(N)=O"},
      {"Cc1cc2c(cc1C)N3C=[N]2[Co]456(N7[C@@H]8[C@@H]([C@](C7=C(C9=[N]4C(=CC1=["
       "N]5C(=C(C2=[N]6[C@@]8([C@@]([C@@H]2CCC(=O)N)(C)CC(=O)N)C)C)[C@@]([C@@H]"
       "1/C=C/"
       "C(=O)N)(C)CC(=O)N)C([C@@H]9CCC(=O)N)(C)C)C)(CCC(=O)NC[C@H](OP(=O)(O[C@@"
       "H]1[C@H](O[C@H]3[C@@H]1O)CO)O)C)C)CC(=O)N)O",
       "CC1=C2N=C(C=C3N=C(C(C)=C4N[C@H]([C@H](CC(N)=O)[C@@]4(C)CCC(=O)NC[C@@H]("
       "C)OP(=O)(O)O[C@@H]4[C@@H](CO)O[C@H](n5cnc6cc(C)c(C)cc65)[C@@H]4O)[C@]4("
       "C)N=C1[C@@H](CCC(N)=O)[C@]4(C)CC(N)=O)[C@@H](CCC(N)=O)C3(C)C)[C@@H](/"
       "C=C/C(N)=O)[C@]2(C)CC(N)=O"},
      {"Cc1cc2c(cc1C)n(cn2)C3C(C(C(O3)CO)OP(=O)(O)OC(C)CNC(=O)CCC4(C(C5C6(C(C("
       "C7=[N]6[Co+2]89[N]5=C4C(=C1[NH]8C(=CC2=[N]9C(=C7C)C(C2CCC(=O)N)(C)CC(="
       "O)N)C(C1CCC(=O)N)(C)C)C)CCC(=O)N)(C)CC(=O)N)C)CC(=O)N)C)O",
       "CC1=C2NC(=CC3=NC(=C(C)C4=NC(C)(C5N=C1C(C)(CCC(=O)NCC(C)OP(=O)(O)OC1C("
       "CO)OC(n6cnc7cc(C)c(C)cc76)C1O)C5CC(N)=O)C(C)(CC(N)=O)C4CCC(N)=O)C(C)("
       "CC(N)=O)C3CCC(N)=O)C(C)(C)C2CCC(N)=O"},
      {"Cc1cc2c(cc1C)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)CO)O[P@@](=O)(O)O[C@H]"
       "(C)CNC(=O)CC[C@@]4([C@H](C5=[N-]6C4=C(C7=[N]8[Co+2]69([N]1=C(C=C8C([C@@"
       "H]7CCC(=O)N)(C)C)[C@H]([C@](C1=C(C1=[N]9[C@@]5([C@@]([C@@H]1CCC(=O)N)("
       "C)CC(=O)N)C)C)(C)CC(=O)N)CCC(=O)N)C)C)CC(=O)N)C)O",
       "CC1=C2N=C(C=C3N=C(C(C)=C4N=C([C@H](CC(N)=O)[C@@]4(C)CCC(=O)NC[C@@H](C)"
       "O[P@](=O)(O)O[C@@H]4[C@@H](CO)O[C@H](n5cnc6cc(C)c(C)cc65)[C@@H]4O)[C@]"
       "4(C)N=C1[C@@H](CCC(N)=O)[C@]4(C)CC(N)=O)[C@@H](CCC(N)=O)C3(C)C)[C@@H]("
       "CCC(N)=O)[C@]2(C)CC(N)=O"},
      {"Cc1cc2c(cc1C)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)CO)O[P@@](=O)(O)O[C@H]"
       "(C)CNC(=O)CC[C@@]4([C@H](C5=[N-]6C4=C(C7=[N]8[Co+2]69([N]1=C(C=C8C([C@@"
       "H]7CCC(=O)N)(C)C)[C@H]([C@](C1=C(C1=[N]9[C@@]5([C@@]([C@@H]1CCC(=O)N)("
       "C)CC(=O)N)C)C)(C)CC(=O)N)CCC(=O)N)C#N)C)CC(=O)N)C)O",
       "CC1=C2N=C(C=C3N=C(C(C)=C4N=C([C@H](CC(N)=O)[C@@]4(C)CCC(=O)NC[C@@H](C)"
       "O[P@](=O)(O)O[C@@H]4[C@@H](CO)O[C@H](n5cnc6cc(C)c(C)cc65)[C@@H]4O)[C@]"
       "4(C)N=C1[C@@H](CCC(N)=O)[C@]4(C)CC(N)=O)[C@@H](CCC(N)=O)C3(C)C)[C@@H]("
       "CCC(N)=O)[C@]2(C)CC(N)=O"},
      {"Cc1cc2c(cc1C)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)CO)O[P@@](=O)(O)O[C@H]"
       "(C)CNC(=O)CC[C@@]4([C@H]([C@@H]5[C@]6([C@@]([C@@H](C7=C([C@H]8[C@@]([C@"
       "@H](C9=CC1C([C@@H](C2=C(C4[N-]5[Co+2](N12)(N89)(N76)"
       "CCCCCn1cnc2c1ncnc2N)C)CCC(=O)N)(C)C)CCC(=O)N)(C)CC(=O)N)C)CCC(=O)N)(C)"
       "CC(=O)N)C)CC(=O)N)C)O",
       "CC1=C2NC(C=C3N[C@@H](C(C)=C4N[C@@](C)([C@@H]5NC1[C@](C)(CCC(=O)NC[C@@H]"
       "(C)O[P@](=O)(O)O[C@@H]1[C@@H](CO)O[C@H](n6cnc7cc(C)c(C)cc76)[C@@H]1O)["
       "C@H]5CC(N)=O)[C@@](C)(CC(N)=O)[C@@H]4CCC(N)=O)[C@@](C)(CC(N)=O)[C@@H]"
       "3CCC(N)=O)C(C)(C)[C@@H]2CCC(N)=O"},
      {"CC1=C(C2=[N]3C1=Cc4c(c(c5n4[Fe]36[N]7=C(C=C8N6C(=C2)C(=C8C)CCC(=O)O)C(="
       "C(C7=C5)C=C)C=C)C=C)C=C)CCC(=O)O",
       "C=CC1=C(C=C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5CCC(=O)O)C(CCC(=O)O)="
       "C4C)c(C=C)c3C=C"},
      {"Cc1cc2n3c1C=C4C=C(C5=[N]4[Fe]36[N]7=C(C=C8N6C(=C5)C(=C8CCC(=O)O)C)C(=C("
       "C7=C2)C)CCC(=O)O)C",
       "CC1=Cc2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5CCC(=O)O)C(CCC(=O)O)=C4C)"
       "cc3C"},
      {"Cc1c(c2n3c1C=C4C(=C(C5=[N]4[Fe]36[N]7=C(C=C8N6C(=C5)C(=C8C(=O)O)C)C(=C("
       "C7=C2)C)C(=O)O)C)C)C",
       "CC1=C(C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5C)C(C)=C4C(=O)O)c(C(=O)O)"
       "c3C"},
      {"Cc1c(c2n3c1C(=C4C(=C(C5=[N]4[Zn]36[N]7=C(C=C8N6C(=C5)C(=C8C=C)C)C(=C("
       "C7=C2)CCC(=O)O)C)C=C)C)O)CCC(=O)O",
       "C=CC1=C(C)c2nc1cc1[nH]c(cc3nc(cc4[nH]c(c(C)c4CCC(=O)O)c2O)C(CCC(=O)O)="
       "C3C)c(C=C)c1C"},
      {"Cc1cc(c2c(c1)nc3c4c5c6c(c3n2)C=CC=[N]6[Ru]78([N]5=CC=C4)([N]9=CC="
       "Nc1c9c2c(cc1)N=CC=[N]72)[N]1=CC=Nc2c1c1c(cc2)N=CC=[N]81)C",
       "Cc1cc(C)c2nc3c4cccnc4c4ncccc4c3nc2c1"},
      {"Cc1ccc2c(c1)nc3c4c5c6c(c3n2)C=CC=[N]6[Ru]78([N]5=CC=C4)([N]9=CC="
       "Nc1c9c2c(cc1)N=CC=[N]72)[N]1=CC=Nc2c1c1c(cc2)N=CC=[N]81",
       "Cc1ccc2nc3c4cccnc4c4ncccc4c3nc2c1"},
      {"CC1=CC=C2C=CC3=CC=C([N]4=C3C2=[N]1[Pt]45([C][C]5)C)C",
       "Cc1ccc2ccc3ccc(C)nc3c2n1"},
      {"CC1=CC=C2C=CC3=CC=C([N]4=C3C2=[N]1[Pt]45(CC5)[N]6=CN(C(=C6)c7cn(nn7)[C@"
       "@H]8[C@H]([C@H]([C@H]([C@H](O8)O)O)O)O)C)C",
       "Cn1cncc1-c1cn([C@H]2O[C@H](O)[C@H](O)[C@H](O)[C@@H]2O)nn1"},
      {"Cc1cccc2c1nc3c(n2)C4=CC=C[N]5=C4C6=C3C=CC=[N]6[Ru+2]578([N]9=CC="
       "Nc1c9c2c(cc1)N=CC=[N]72)[N]1=CC=Nc2c1c1c(cc2)N=CC=[N]81",
       "Cc1cccc2nc3c4cccnc4c4ncccc4c3nc12"},
      {"CC1=CC=CC2=[N]1[Cu]3[N](=C2)NC(=[S]3)N(C)C",
       "Cc1cccc(C=NNC(=S)N(C)C)n1"},
      {"Cc1ccc(cc1)C(C)C.c1ccc(cc1)NC2=[S][Os]([N]3=C2C=CC=C3)Cl",
       "S=C(Nc1ccccc1)c1ccccn1"},
      {"Cc1cccc(c1)[Fe]234n5c6c(c(c5C=C7[N]2=C(C=C8N3C(=CC9=[N]4C(=C6)C(=C9C)C="
       "C)C(=C8CCC(=O)O)C)C(=C7C)CCC(=O)O)C=C)C",
       "C=CC1=C(C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5C=C)C(C)=C4CCC(=O)O)c("
       "CCC(=O)O)c3C"},
      {"CC1=CC=[N]2C(=C1)C3=CC(=CC=[N]3[Fe]2(C#N)(C#N)(C#N)C#N)CNC(=O)CCCC[C@H]"
       "4[C@@H]5[C@H](CS4)NC(=O)N5",
       "Cc1ccnc(-c2cc(CNC(=O)CCCC[C@@H]3SC[C@@H]4NC(=O)N[C@H]34)ccn2)c1"},
      {"CC1=[CH]2C=C3([Ru+]24([CH]1=C3)([NH2][C@@H]5CN(C[C@H]5[NH2]4)C(=O)CCCC["
       "C@H]6[C@@H]7[C@H](CS6)NC(=O)N7)Cl)C(C)C",
       "N[C@@H]1CN(C(=O)CCCC[C@@H]2SC[C@@H]3NC(=O)N[C@H]23)C[C@H]1N"},
      {"Cc1nc2cc3c(cc2c(n1)N)C4=CC=CC=[N]4[Ru]3([C-]#[O+])C5(C=CC=C5)C(=O)Oc6c("
       "cccc6OC)OC",
       "COc1cccc(OC)c1OC(=O)C1C=CC=C1"},
      {"CC1O[Co]23(O4[Co]5(O1)(O2[Co]67(O3[Co]4(O56)(OC(O7)CCNC(=O)CCCC[C@H]8["
       "C@@H]9[C@H](CS8)NC(=O)N9)([N]1=CC=CC=C1)[O])([N]1=CC=CC=C1)[O])([N]1="
       "CC=CC=C1)[O])[N]1=CC=CC=C1",
       "O=C(CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12)NCCC(O)O"},
      {"CC1O[Co]23(O4[Co]5(O1)(O2[Co]67(O3[Co]4(O56)(OC(O7)CCNC(=O)CCCC[C@H]8["
       "C@@H]9[C@H](CS8)NC(=O)N9)([N]1=CC=CC=C1)[O])([N]1=CC=CC=C1)[O])[O])[N]"
       "1=CC=CC=C1",
       "O=C(CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12)NCCC(O)O"},
      {"CC1=O[Ga]2345O=C(N(O2)CCC[C@H]6C(=O)NCC(=O)N[C@H](C(=O)NCC(=O)N[C@@H]("
       "CCCN1O3)C(=O)N[C@@H](CCCN(O4)C(=O5)C)C(=O)N6)CO)C",
       "CC(=O)N(O)CCC[C@@H]1NC(=O)CNC(=O)[C@H](CO)NC(=O)CNC(=O)[C@H](CCCN(O)C("
       "C)=O)NC(=O)[C@H](CCCN(O)C(C)=O)NC1=O"},
      {"C\\C=C/"
       "1\\c2cc3c(c4c5n3[Mg]67[n+]2c(cc8n6c(cc9[n+]7c(c5[C@H](C4=O)C(=O)OC)[C@"
       "H]([C@@H]9C)CCC(=O)OC\\C=C(/"
       "C)\\CCC[C@H](C)CCC[C@H](C)CCCC(C)C)c(c8C(=O)C)C)[C@@H]1C)C",
       "C/C=C1/"
       "c2cc3[nH]c4c(c3C)C(=O)[C@H](C(=O)OC)c4c3nc(cc4[nH]c(cc(n2)[C@@H]1C)c(C("
       "C)=O)c4C)[C@@H](C)[C@@H]3CCC(=O)OC/"
       "C=C(\\C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C"},
      {"C#CC1=C2C=CC3=C(C4=[N]5C(=Cc6ccc7n6[Fe]5(N23)[N]8=C1C=CC8=C7c9ccccc9)C="
       "C4)c1ccccc1",
       "C#Cc1c2nc(c(-c3ccccc3)c3ccc(cc4nc(c(-c5ccccc5)c5ccc1[nH]5)C=C4)[nH]3)C="
       "C2"},
      {"C#CC1=C2C=CC3=[N]2[Fe]45n6c1ccc6C(=C7[N]4=C(C=C7)C(=C8N5C(=C3c9ccccc9)"
       "C=C8)C#C)c1ccccc1",
       "C#Cc1c2nc(c(-c3ccccc3)c3ccc([nH]3)c(C#C)c3nc(c(-c4ccccc4)c4ccc1[nH]4)C="
       "C3)C=C2"},
      {"CCC1=C(C2=CC3=C(C(=C4[N@]3([Cu@@]56[N]2=C1C=C7N5C(=CC8=[N]6C(=C4)C(="
       "C8CCC(=O)O)C)C(=C7C)CCC(=O)O)C)C)CC)C",
       "CCC1=C(C)c2cc3c(CC)c(C)c(cc4nc(cc5[nH]c(cc1n2)c(C)c5CCC(=O)O)C(CCC(=O)"
       "O)=C4C)n3C"},
      {"CCC1=C(C2=Cc3c(c(c4n3[Fe]56[N]2=C1C=C7N5C8=C(CC(=O)C8=C7C)C9=[N]6C(=C4)"
       "[C@H]([C@@H]9CCC(=O)OC)C)C)C=C)C",
       "C=Cc1c(C)c2cc3nc(c4c5[nH]c(cc6nc(cc1[nH]2)C(C)=C6CC)c(C)c5C(=O)C4)[C@@"
       "H](CCC(=O)OC)[C@@H]3C"},
      {"CCC1=C(C2=Cc3c(c(c4n3[Mg]56N2C1=Cc7n5c8c(c7C)C(=O)[C@@H](C8=C9N6C(=C4)["
       "C@H]([C@@H]9CCC(=O)OC/C=C(\\C)/"
       "CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C)C(=O)OC)C)C=O)C",
       "CCc1c(C)c2[nH]c1=Cc1[nH]c3c(c1C)C(=O)[C@H](C(=O)OC)C3=C1NC(=Cc3[nH]c(c("
       "C=O)c3C)C=2)[C@@H](C)[C@@H]1CCC(=O)OC/"
       "C=C(\\C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C"},
      {"CCC1=C(c2cc3c(c(c4n3[Mg]56[n+]2c1cc7n5c8c(c9[n+]6c(c4)C(C9CCC(=O)OC/"
       "C=C(\\C)/"
       "CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C)[C@@H](C(=O)c8c7C)C(=O)OC)C)C=C)C",
       "C=Cc1c(C)c2cc3nc(c4c5[nH]c(cc6nc(cc1[nH]2)C(C)=C6CC)c(C)c5C(=O)[C@H]4C("
       "=O)OC)C(CCC(=O)OC/C=C(\\C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C3C"},
      {"CCC1=C(c2cc3c(c(c4n3[Mg]56[n+]2c1cc7n5c8c(c9[n+]6c(c4)C(C9CCC(=O)OC/"
       "C=C(\\C)/"
       "CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C)[C@H](C(=O)c8c7C)C(=O)OC)C)C=C)C=O",
       "C=Cc1c(C)c2cc3nc(c4c5[nH]c(cc6nc(cc1[nH]2)C(C=O)=C6CC)c(C)c5C(=O)[C@@H]"
       "4C(=O)OC)C(CCC(=O)OC/C=C(\\C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C3C"},
      {"CCC1=C(c2cc3c(c(c4n3[Mg]56[n+]2c1cc7n5c8c(c9[n+]6c(c4)C(=C9CCC(=O)O)C)["
       "C@H](C(=O)c8c7C)C(=O)OC)C)C=C)C",
       "C=Cc1c(C)c2cc3nc(c4c5[nH]c(cc6nc(cc1[nH]2)C(C)=C6CC)c(C)c5C(=O)[C@@H]"
       "4C(=O)OC)C(CCC(=O)O)=C3C"},
      {"CCC1=C(c2cc3c(c(c4n3[Mg]56[n+]2c1cc7n5c8c(c9[n+]6c(c4)[C@H]([C@@H]9CCC("
       "=O)OC/C=C(\\C)/"
       "CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C)[C@H](C(=O)c8c7C)C(=O)OC)C)C=C)C",
       "C=Cc1c(C)c2cc3nc(c4c5[nH]c(cc6nc(cc1[nH]2)C(C)=C6CC)c(C)c5C(=O)[C@@H]"
       "4C(=O)OC)[C@@H](CCC(=O)OC/"
       "C=C(\\C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)[C@@H]3C"},
      {"CCC1=C(c2cc3c(c(c4n3[Mg]56[n+]2c1cc7n5c8c(c9[n+]6c(c4)[C@H]([C@@H]9CCC("
       "=O)OC/C=C(\\C)/"
       "CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C)[C@H](C(=O)c8c7C)C(=O)OC)C)C=C)C",
       "C=Cc1c(C)c2cc3nc(c4c5[nH]c(cc6nc(cc1[nH]2)C(C)=C6CC)c(C)c5C(=O)[C@@H]"
       "4C(=O)OC)[C@@H](CCC(=O)OC/"
       "C=C(\\C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)[C@@H]3C"},
      {"CCC1=C(c2cc3c(c(c4n3[Mg]56[n+]2c1cc7n5c8c(c9[n+]6c(c4)[C@H]([C@@H]9CCC("
       "=O)OC/C=C(\\C)/"
       "CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C)[C@H](C(=O)c8c7C)C(=O)OC)C)C=C)C",
       "C=Cc1c(C)c2cc3nc(c4c5[nH]c(cc6nc(cc1[nH]2)C(C)=C6CC)c(C)c5C(=O)[C@@H]"
       "4C(=O)OC)[C@@H](CCC(=O)OC/"
       "C=C(\\C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)[C@@H]3C"},
      {"CCC1=C(C2=Cc3c(c(c4n3[Zn]56[N]2=C1C=C7N5C8=C(C(=C(C8=C7C)O)C(=O)OC)C9=["
       "N]6C(=C4)C(=C9CCC(=O)OC/C=C(\\C)/"
       "CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C)C)C(=O)C)C",
       "CCC1=C(C)c2cc3[nH]c(cc4nc(c5c6[nH]c(cc1n2)c(C)c6C(O)=C5C(=O)OC)C(CCC(="
       "O)OC/C=C(\\C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)=C4C)c(C)c3C(C)=O"},
      {"CCC1=C(C2=CC3=[N]4C(=CC5=[N]6[Fe]47[N]2=C1C=C8[N]7=C(C=C6C(=C5CCC(=O)O)"
       "C)C(=C8C)CC)C(=C3C)CCC(=O)O)C",
       "CCC1=C(C)C2=CC3=NC(=CC4=NC(=CC5=NC(=CC1=N2)C(C)=C5CC)C(C)=C4CCC(=O)O)C("
       "CCC(=O)O)=C3C"},
      {"CCc1c(c2n3c1C=C4C(=C5C6=[N]4[Mg]37[N]8=C(C=C9N7C(=C6[C@](C5=O)(C(=O)OC)"
       "O)C(=C9C)CCC(=O)OCC=C(C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C(=C(C8=C2)CC)C)"
       "C)C",
       "CCC1=C(C)c2cc3[nH]c(c(CCC(=O)OCC=C(C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)"
       "c3C)c3c4nc(cc5[nH]c(cc1n2)c(C)c5CC)C(C)=C4C(=O)[C@]3(O)C(=O)OC"},
      {"CCc1c(c2n3c1C=C4C(=C5C(=O)[C@@H](C6=C7C(=C(C8=CC9=[N]([Mg]3(N87)[N]4="
       "C65)C(=C2)C(=C9C)CC)C)CCC(=O)OC/C=C(/C)\\CC/C=C(\\C)/"
       "CCC=C(C)C)C(=O)OC)C)C",
       "CCC1=C(C)c2cc3[nH]c(c(CCC(=O)OC/C=C(/C)CC/"
       "C=C(\\C)CCC=C(C)C)c3C)c3c4nc(cc5[nH]c(cc1n2)c(C)c5CC)C(C)=C4C(=O)[C@@H]"
       "3C(=O)OC"},
      {"CCc1c(c2n3c1C=C4C(=C5C(=O)[C@H](C6=C7C(=C(C8=CC9=[N]([Mg]3(N87)[N]4="
       "C65)C(=C2)C(=C9C)CC)C)CCC(=O)OC/C=C(\\C)/CC/C=C(\\C)/"
       "CCC=C(C)C)C(=O)OC)C)C",
       "CCC1=C(C)c2cc3[nH]c(c(CCC(=O)OC/C=C(\\C)CC/"
       "C=C(\\C)CCC=C(C)C)c3C)c3c4nc(cc5[nH]c(cc1n2)c(C)c5CC)C(C)=C4C(=O)[C@H]"
       "3C(=O)OC"},
      {"CCc1c(c2n3c1C=C4C(=C5C(=O)[C@@H](C6=C7C(=C(C8=CC9=[N]([Mg]3(N87)[N]4="
       "C65)C(=C2)C(=C9C=O)C=C)C)CCC(=O)OC/C=C(\\C)/"
       "CCC[C@H](C)CCCC(C)CCCC(C)C)C(=O)OC)C)C",
       "C=CC1=C(C=O)c2cc3[nH]c(c(CCC(=O)OC/"
       "C=C(\\C)CCC[C@H](C)CCCC(C)CCCC(C)C)c3C)c3c4nc(cc5[nH]c(cc1n2)c(C)c5CC)"
       "C(C)=C4C(=O)[C@@H]3C(=O)OC"},
      {"CCC1=C(C2=[N]3C1=Cc4c(c(c5n4[Fe]36[N]7=C(C=C8N6C(=C2)C(=C8CCC(=O)O)C)C("
       "=C(C7=C5)C)CCC(=O)O)C=C)C)C",
       "C=Cc1c(C)c2cc3nc(cc4[nH]c(cc5nc(cc1[nH]2)C(C)=C5CCC(=O)O)c(CCC(=O)O)"
       "c4C)C(C)=C3CC"},
      {"CCc1c(c2n3c1C=C4C(=C(C5=[N]4[Fe]36[N]7=C(C=C8N6C(=C5)C(=C8CCC(=O)O)C)C("
       "=C(C7=C2)C)CCC(=O)O)CC)C(F)(F)F)C",
       "CCC1=C(C(F)(F)F)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5CCC(=O)O)C(CCC(="
       "O)O)=C4C)c(C)c3CC"},
      {"CCC1=C(C2=[N-]3C1=Cc4c(c(c5n4[Mg]36[N-]7=C8C(=C9N6C(=C2)[C@H]([C@@H]"
       "9CCC(=O)OC/C=C(\\C)/CC/C=C(\\C)/"
       "CCC=C(C)C)C)[C@H](C(=O)C8=C(C7=C5)C)C(=O)OC)[C@H](C)O)C)C",
       "CCC1=C(C)c2cc3[nH]c(c4c5nc(cc6[nH]c(cc1n2)c(C)c6[C@H](C)O)C(C)=C5C(=O)["
       "C@@H]4C(=O)OC)[C@@H](CCC(=O)OC/C=C(\\C)CC/C=C(\\C)CCC=C(C)C)[C@@H]3C"},
      {"CCc1c(c2n3c1C=C4C(=C(C5=[N]4[Mo]36(=O)[N]7=C(C=C8N6C(=C5)C(=C8C)CCC(=O)"
       "O)C(=C(C7=C2)CC)C)CCC(=O)O)C)C",
       "CCC1=C(C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5CC)C(C)=C4CCC(=O)O)c("
       "CCC(=O)O)c3C"},
      {"CCc1c(c2n3c1C=C4C(=C(C5=[N]4[Rh]36[N]7=C(C=C8N6C(=C5)C(=C8CCC(=O)O)C)C("
       "=C(C7=C2)C)CCC(=O)O)CC)C)C",
       "CCC1=C(C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5CCC(=O)O)C(CCC(=O)O)="
       "C4C)c(C)c3CC"},
      {"CCc1c(c2n3c1C=C4C(=C(C5=[N]4[Ru@@]36[N]7=C(C=C8N6C(=C5)C(=C8C)CCC(=O)O)"
       "C(=C(C7=C2)CC)C)CCC(=O)O)C)C",
       "CCC1=C(C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5CC)C(C)=C4CCC(=O)O)c("
       "CCC(=O)O)c3C"},
      {"CCc1c(c2n3c1C=C4[C@](C(=O)C5=[N]4[Fe]36[N]7=C(C=C8N6C(=C5)C(=C8CCC(=O)"
       "O)C)C(=C(C7=C2)C)CCC(=O)O)(C)CC)C",
       "CCc1c(C)c2cc3nc(cc4[nH]c(cc5nc(cc1[nH]2)[C@](C)(CC)C5=O)c(C)c4CCC(=O)O)"
       "C(CCC(=O)O)=C3C"},
      {"CCCC12C3([Ir]1456(C3(C4(C52C)C)C)(N(CC7=[N]6C=CC(=C7)c8ccc(cc8)S(=O)(="
       "O)N)S(=O)(=O)c9ccccc9)Cl)C",
       "NS(=O)(=O)c1ccc(-c2ccnc(CNS(=O)(=O)c3ccccc3)c2)cc1"},
      {"CC(C)C12=C3[Os]1456(C2=C4C5(=C63)C)(Cl)Cl", "Cc1ccc(C(C)C)cc1"},
      {"CC(C)C12=C3[Ru]1456(C2=C4C5(=C63)C)([P]78CN9CN(C7)CN(C9)C8)(Cl)Cl",
       "Cc1ccc(C(C)C)cc1"},
      {"CC(C)C12=CC=C([Os]1([N]3=C(C=CC=C3)C(=S)Nc4ccc(cc4)F)Cl)(C=C2)C",
       "Fc1ccc(NC(=S)c2ccccn2)cc1"},
      {"CC(C)C12=[CH]3[Ru]1456([CH]2=[CH]4C5(=[CH]63)C)(Cl)Cl",
       "Cc1ccc(C(C)C)cc1"},
      {"CCCCCCCCCCCC(=O)N1C2=[N](C=CC=C2)[Rh][N]3=CC=CC=C31",
       "CCCCCCCCCCCC(=O)N(c1ccccn1)c1ccccn1"},
      {"CCCCCC[NH2][Pt]([NH3])([NH3])[NH2]CCCCCC", "CCCCCCN"},
      {"CC(C)C[C@H]1C(=O)O[Cu]2[N]1=Cc3ccccc3O2",
       "CC(C)C[C@H](N=Cc1ccccc1O)C(=O)O"},
      {"C(CCC[NH2][Pt]([NH3])([NH3])[NH3])CCN", "NCCCCCCN"},
      {"CC[C@@H]1[C@H](C2=CC3=C(C(=C4[N-]3[Mg@+2]56[N]2=C1C=C7[N-]5C8=C([C@H]("
       "C(=O)C8=C7C)C(=O)OC)C9=[N]6C(=C4)[C@H]([C@@H]9CCC(=O)OC\\C=C(/"
       "C)\\CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C)C)C(=O)C)C",
       "CC[C@H]1c2cc3[nH]c4c(c3C)C(=O)[C@H](C(=O)OC)c4c3nc(cc4[nH]c(cc(n2)[C@@"
       "H]1C)c(C(C)=O)c4C)[C@@H](C)[C@@H]3CCC(=O)OC/"
       "C=C(\\C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C"},
      {"C(CC[NH2][Pt]([NH3])([NH3])Cl)C[NH2][Pt]([NH3])([NH3])Cl", "NCCCCN"},
      {"CC(CO)(CO)NC(=O)C12[CH]3=[CH]4[Ru]3156([CH]4=[CH]52)(n7c8ccc("
       "cc8c9c7c1c(c2c9C(=O)N(C2=O)C)C=C(C=[N]61)F)O)C#O",
       "CN1C(=O)c2c(c3c4cc(O)ccc4[nH]c3c3ncc(F)cc23)C1=O"},
      {"C[CH]12[CH]3([Ir+2]1456([CH]2(C4(=C53C)C)C)([NH2]CC[NH]6S(=O)(=O)c7ccc("
       "cc7)NC(=O)CCCC[C@H]8[C@@H]9[C@H](CS8)NC(=O)N9)Cl)C",
       "NCCNS(=O)(=O)c1ccc(NC(=O)CCCC[C@@H]2SC[C@@H]3NC(=O)N[C@H]23)cc1"},
      {"C([C@@H]1C(=C([C@H]2[C@@H](O1)NC3=C(N2)C(=O)NC(=N3)N)S)S)O[P@](=O)(O)O["
       "Mg]([OH2])([OH2])O[P@@](=O)(O)OC[C@@H]4C(=C([C@H]5[C@@H](O4)NC6=C(N5)C("
       "=O)NC(=N6)N)S[W])S",
       "Nc1nc2c(c(=O)[nH]1)N[C@H]1C(S)=C(S)[C@@H](COP(=O)(O)O)O[C@H]1N2"},
      {"C[C@@H]1C[N]2=CC3=[N]([Ni]24[N]1=CC5=CC=CC=[N]45)C=CC=C3",
       "C[C@H](CN=Cc1ccccn1)N=Cc1ccccn1"},
      {"C([C@H]1C(=O)O[Re+]([NH2]1)(C#O)(C#O)C#O)O", "N[C@@H](CO)C(=O)O"},
      {"CCOC(=O)C12[C@@H]3[Ru]1456([C@@H]3[C@H]4[C@H]52)(n7c8ccc(cc8c9c7c1c("
       "c2c9C(=O)NC2=O)C=CC=[N]61)O)[CH]#O",
       "O=C1NC(=O)c2c1c1cccnc1c1[nH]c3ccc(O)cc3c21"},
      {"CC[O@@H][Ca+2]([OH2])([OH2])([OH2])([OH2])([OH2])[OH2]", "CCO"},
      {"CCO[P@@](=O)(CCCc1cc2c3c(c1)C[N]([Pt]3([N](C2)(C)C)Cl)(C)C)Oc4ccc(cc4)["
       "N+](=O)[O-]",
       "CCO[P@@](=O)(CCCc1cc(CN(C)C)c([Pt+])c(CN(C)C)c1)Oc1ccc([N+](=O)[O-])"
       "cc1"},
      {"[CH]12=[CH]3[Ru]1456[CH]2=[CH]4[CH]5=[CH]63", "c1ccccc1"},
      {"[CH]12=[CH]3[Ru]1456([CH]2=[CH]4[CH]5=[CH]63)(Cl)Cl", "c1ccccc1"},
      {"[CH2]1=[CH2][Pt+2]1([Cl-])([Cl-])[Cl-]", "C=C"},
      {"C[N@]12CC3=CC=CC=[N]3[Fe+2]145[N]6=CC=CC=C6C[N@@]4(CC(=O)O5)[C@@H]7[C@@"
       "H]2CCCC7",
       "CN(Cc1ccccn1)[C@H]1CCCC[C@@H]1N(CC(=O)O)Cc1ccccn1"},
      {"C[N]12CC[N]34[Fe+2]1([N]5=C(C2)C=CC=C5)([N]6=CC=CC=C6C3)OC(=O)C4",
       "CN(CCN(CC(=O)O)Cc1ccccn1)Cc1ccccn1"},
      {"C[N@@]12[C@H]3CS[C@@H]1[C@H]4CSC5=[N]4[Fe]2(Oc6c5cccc6)OC3=O",
       "CN1[C@H](C(=O)O)CS[C@@H]1[C@H]1CSC(c2ccccc2O)=N1"},
      {"CNC1C[S]2CC[S]3[Ru]24(n5c6ccc(cc6c7c5c8c(c9c7C(=O)NC9=O)C=C(C=[N]48)F)"
       "OC)([S](C1)CC3)[N]CS",
       "COc1ccc2[nH]c3c4ncc(F)cc4c4c(c3c2c1)C(=O)NC4=O"},
      {"CNC(=[S][Pt+2]1[NH2]CC[NH2]1)N(C)CCNc2c3ccccc3[nH+]c4c2cccc4",
       "CNC(=S)N(C)CCNc1c2ccccc2nc2ccccc12"},
      {"COc1cc(cc(c1O)OC)[C@@H]2c3cc4c(cc3[C@@H]([C@@H]5[C@@H]2C(=O)OC5)NC(=O)"
       "CC[C@@H]6C[NH2][Pt]([NH2]6)(Cl)Cl)OCO4",
       "COc1cc([C@@H]2c3cc4c(cc3[C@H](NC(=O)CC[C@@H](N)CN)[C@H]3COC(=O)[C@@H]"
       "32)OCO4)cc(OC)c1O"},
      {"COc1cc(cc(c1O)OC)[C@@H]2c3cc4c(cc3[C@@H]([C@@H]5[C@@H]2C(=O)OC5)NC(=O)"
       "CC[C@H]6C[NH2][Pt]([NH2]6)(Cl)Cl)OCO4",
       "COc1cc([C@@H]2c3cc4c(cc3[C@H](NC(=O)CC[C@H](N)CN)[C@H]3COC(=O)[C@@H]32)"
       "OCO4)cc(OC)c1O"},
      {"COc1cc(cc(c1O)OC)[C@@H]2c3cc4c(cc3[C@H]([C@@H]5[C@H]2C(=O)OC5)NC(=O)CC["
       "C@@H]6C[NH2][Pt]([NH2]6)Cl)OCO4",
       "COc1cc([C@@H]2c3cc4c(cc3[C@@H](NC(=O)CC[C@@H](N)CN)[C@H]3COC(=O)[C@H]"
       "32)OCO4)cc(OC)c1O"},
      {"COc1cc(cc(c1O)OC)[C@@H]2c3cc4c(cc3[C@H]([C@@H]5[C@H]2C(=O)OC5)NC(=O)CC["
       "C@H]6C[NH2][Pt]([NH2]6)Cl)OCO4",
       "COc1cc([C@@H]2c3cc4c(cc3[C@@H](NC(=O)CC[C@H](N)CN)[C@H]3COC(=O)[C@H]32)"
       "OCO4)cc(OC)c1O"},
      {"C(=O)[Fe](C=O)(C#N)O(=O)[Ni]", "C=O"},
      {"[C](=O)[Re]12(O3[Re]4(O1[Re]5(O2[Re]3(O45)([C]=O)([C]=O)[C]=O)([C]=O)(["
       "C]=O)[C]=O)([C]=O)([C]=O)[C]=O)([C]=O)[C]=O",
       "[C]=O"},
      {"[H]/N=C(\\c1ccc2c(c1)NC3=[N]2[Zn+2][N]4=C(C3)Nc5c4ccc(c5)C(=N)N)/N",
       "[H]/N=C(/N)c1ccc2nc(Cc3nc4ccc(C(=N)N)cc4[nH]3)[nH]c2c1"},
      {"[NH3][Co+3]([NH3])([NH3])[NH3]", "N"},
      {"[NH3][Co+3]([NH3])([NH3])([NH3])([NH3])[NH3]", "N"},
      {"[NH3][Ir+3]([NH3])([NH3])([NH3])([NH3])[NH3]", "N"},
      {"[NH3][Pt+2]([NH3])[NH3]", "N"},
      {"[NH3][Pt]([NH3])(Cl)Cl", "N"},
      {"[NH3][Rh+3]([NH3])([NH3])([NH3])([NH3])[NH3]", "N"},
      {"[NH3][Ru+3]([NH3])([NH3])([NH3])([NH3])[NH3]", "N"},
      {"[OH2][Ca+2]", "O"},
      {"[OH2][Ca+2][OH2]", "O"},
      {"[OH2][Ca+2]([OH2])[OH2]", "O"},
      {"[OH2][Ca+2]([OH2])([OH2])[OH2]", "O"},
      {"[OH2][Ca+2]([OH2])([OH2])([OH2])[OH2]", "O"},
      {"[OH2][Ca+2]([OH2])([OH2])([OH2])([OH2])[OH2]", "O"},
      {"[OH2][Ca+2]([OH2])([OH2])([OH2])([OH2])([OH2])[OH2]", "O"},
      {"[OH2][Ca+2]([OH2])([OH2])([OH2])([OH2])([OH2])([OH2])[OH2]", "O"},
      {"[OH2][K]([OH2])([OH2])[OH2]", "O"},
      {"[OH2][Mg+2]", "O"},
      {"[OH2][Mg+2][OH2]", "O"},
      {"[OH2][Mg+2]([OH2])[OH2]", "O"},
      {"[OH2][Mg+2]([OH2])([OH2])[OH2]", "O"},
      {"[OH2][Mg+2]([OH2])([OH2])([OH2])[OH2]", "O"},
      {"[OH2][Mg+2]([OH2])([OH2])([OH2])([OH2])[OH2]", "O"},
      {"[OH2-][Na+]", "O"},
      {"[OH2-][Na+][OH2-]", "O"},
      {"[OH2-][Na+]([OH2-])[OH2-]", "O"},
      {"[OH2-][Na+]([OH2-])([OH2-])([OH2-])[OH2-]", "O"},
      {"[OH2-][Na+]([OH2-])([OH2-])([OH2-])([OH2-])[OH2-]", "O"},
  };

  std::string normalizationSmarts = R"DATA(//	Name	SMIRKS
// Opposite of #2.1 in InChI technical manual? Covered by RDKit Sanitization.
Nitro to N+(O-)=O	[N,P,As,Sb;X3:1](=[O,S,Se,Te:2])=[O,S,Se,Te:3]>>[*+1:1]([*-1:2])=[*:3]
Sulfone to S(=O)(=O)	[S+2:1]([O-:2])([O-:3])>>[S+0:1](=[O-0:2])(=[O-0:3])
Pyridine oxide to n+O-	[n:1]=[O:2]>>[n+:1][O-:2]
Azide to N=N+=N-	[N:2]=[N:3]#[N:4]>>[N:2]=[N+:3]=[N-:4]
Broken azide to N=N+=N-	[N:2]=[N:3]=[N:4]>>[NH0:2]=[NH0+:3]=[NH0-:4]
Diazo/azo to =N+=N-	[*:1]=[N:2]#[N:3]>>[*:1]=[N+:2]=[N-:3]
Sulfoxide to -S+(O-)-	[!O:1][S+0;X3:2](=[O:3])[!O:4]>>[*:1][S+1:2]([O-:3])[*:4]
// Equivalent to #1.5 in InChI technical manual
Phosphate to P(O-)=O	[O,S,Se,Te;-1:1][P+;D4:2][O,S,Se,Te;-1:3]>>[*+0:1]=[P+0;D5:2][*-1:3]
// Equivalent to #1.8 in InChI technical manual
C/S+N to C/S=N+	[C,S;X3+1:1]([NX3:2])[NX3!H0:3]>>[*+0:1]([N:2])=[N+:3]
// Equivalent to #1.8 in InChI technical manual
P+N to P=N+	[P;X4+1:1]([NX3:2])[NX3!H0:3]>>[*+0:1]([N:2])=[N+:3]
Normalize hydrazine-diazonium	[CX4:1][NX3H:2]-[NX3H:3][CX4:4][NX2+:5]#[NX1:6]>>[CX4:1][NH0:2]=[NH+:3][C:4][N+0:5]=[NH:6]
// Equivalent to #1.3 in InChI technical manual
Recombine 1,3-separated charges	[N,P,As,Sb,O,S,Se,Te;-1:1]-[A+0:2]=[N,P,As,Sb,O,S,Se,Te;+1:3]>>[*-0:1]=[*:2]-[*+0:3]
Recombine 1,3-separated charges	[n,o,p,s;-1:1]:[a:2]=[N,O,P,S;+1:3]>>[*-0:1]:[*:2]-[*+0:3]
Recombine 1,3-separated charges	[N,O,P,S;-1:1]-[a:2]:[n,o,p,s;+1:3]>>[*-0:1]=[*:2]:[*+0:3]
Recombine 1,5-separated charges	[N,P,As,Sb,O,S,Se,Te;-1:1]-[A+0:2]=[A:3]-[A:4]=[N,P,As,Sb,O,S,Se,Te;+1:5]>>[*-0:1]=[*:2]-[*:3]=[*:4]-[*+0:5]
Recombine 1,5-separated charges	[n,o,p,s;-1:1]:[a:2]:[a:3]:[c:4]=[N,O,P,S;+1:5]>>[*-0:1]:[*:2]:[*:3]:[c:4]-[*+0:5]
Recombine 1,5-separated charges	[N,O,P,S;-1:1]-[c:2]:[a:3]:[a:4]:[n,o,p,s;+1:5]>>[*-0:1]=[c:2]:[*:3]:[*:4]:[*+0:5]
// Conjugated cation rules taken from Francis Atkinson's standardiser. Those that can reduce aromaticity aren't included
Normalize 1,3 conjugated cation	[N,O;+0!H0:1]-[A:2]=[N!$(*[O-]),O;+1H0:3]>>[*+1:1]=[*:2]-[*+0:3]
Normalize 1,3 conjugated cation	[n;+0!H0:1]:[c:2]=[N!$(*[O-]),O;+1H0:3]>>[*+1:1]:[*:2]-[*+0:3]
//Normalization('Normalize 1,3 conjugated cation', '[N,O;+0!H0:1]-[c:2]:[n!$(*[O-]),o;+1H0:3]>>[*+1:1]=[*:2]:[*+0:3]'),
Normalize 1,5 conjugated cation	[N,O;+0!H0:1]-[A:2]=[A:3]-[A:4]=[N!$(*[O-]),O;+1H0:5]>>[*+1:1]=[*:2]-[*:3]=[*:4]-[*+0:5]
Normalize 1,5 conjugated cation	[n;+0!H0:1]:[a:2]:[a:3]:[c:4]=[N!$(*[O-]),O;+1H0:5]>>[n+1:1]:[*:2]:[*:3]:[*:4]-[*+0:5]
// Normalization('Normalize 1,5 conjugated cation', '[N,O;+0!H0:1]-[c:2]:[a:3]:[a:4]:[n!$(*[O-]),o;+1H0:5]>>[*+1:1]=[c:2]:[*:3]:[*:4]:[*+0:5]'),
// Normalization('Normalize 1,5 conjugated cation', '[n;+0!H0:1]1:[a:2]:[a:3]:[a:4]:[n!$(*[O-]);+1H0:5]1>>[n+1:1]1:[*:2]:[*:3]:[*:4]:[n+0:5]1'),
// Normalization('Normalize 1,5 conjugated cation', '[n;+0!H0:1]:[a:2]:[a:3]:[a:4]:[n!$(*[O-]);+1H0:5]>>[n+1:1]:[*:2]:[*:3]:[*:4]:[n+0:5]'),
// Equivalent to #1.6 in InChI technical manual. RDKit Sanitization handles this for perchlorate.
Charge normalization	[F,Cl,Br,I,At;-1:1]=[O:2]>>[*-0:1][O-:2]
Charge recombination	[N,P,As,Sb;-1:1]=[C+;v3:2]>>[*+0:1]#[C+0:2]
Tetravalent [N] to [N+]	[N+0;v4:1]>>[N+1:1]
Tetravalent [N-] to [N]	[N-1;v4:1]>>[N+0:1]
[OH-] to [O-]	[N,P,As,Sb,O,S,Se,Te:1][OH-:2]>>[N,P,As,Sb,O,S,Se,Te:1][OH0;1-:2]
Duplicate C=O	[C:1](=[O:2])=O>>[C:1]=[O:2]
CH-C=[O,S] to C=C-[O,S] in aromatic ring	[#6;H1:1][#6:2]=[#8,#16;+0:3][*:4]>>[CH0:1]=[C:2][*:3][*:4]
CH2-C=[O,S] to CH=C-[O,S] in aromatic ring	[#6;H2:1][#6:2]=[#8,#16;+0:3][*:4]>>[CH1:1]=[C:2][*:3][*:4]
Negatively charged divalent O,S	[O,S;v2;-1:1]>>[*;H0;-0:1]
Neutral trivalent O,S	[O,S;v3;+0:1]>>[*;+1:1]
Broken 1,3-dioxol	[#8,#16;v3,v4:1]=[#6:2][#8,#16:3]>>[*;+0;H0:1][C:2][*;+0;H0:3]
Positively charged tetravalent B	[B;v4;+1:1]>>[*;-1:1])DATA";

  std::stringstream normalizationSmartsStream(normalizationSmarts);
  MolStandardize::CleanupParameters params;
  MolStandardize::Normalizer normalizer(normalizationSmartsStream,
                                        params.maxRestarts);
  MolStandardize::MetalDisconnector metalDisconnector;
  MolStandardize::Reionizer reionizer;
  MolStandardize::LargestFragmentChooser largestFragmentChooser;
  MolStandardize::Uncharger uncharger;

  unsigned int failedOp;
  for (const auto &pair : ligandExpoSmiles) {
    RWMOL_SPTR rwmol(SmilesToMol(pair.first, 0, false));
    MolOps::sanitizeMol(
        *rwmol, failedOp,
        MolOps::SANITIZE_CLEANUP | MolOps::SANITIZE_FINDRADICALS);
    metalDisconnector.disconnect(*rwmol);
    ROMOL_SPTR mol(normalizer.normalize(static_cast<const ROMol &>(*rwmol)));
    rwmol.reset(new RWMol(*mol));
    MolOps::sanitizeMol(*rwmol);
    mol.reset(reionizer.reionize(static_cast<const ROMol &>(*rwmol)));
    MolOps::assignStereochemistry(*mol);
    mol.reset(largestFragmentChooser.choose(*mol));
    mol.reset(uncharger.uncharge(*mol));
    rwmol.reset(new RWMol(*mol));
    MolOps::sanitizeMol(*rwmol);
    std::unique_ptr<ROMol> refmol(SmilesToMol(pair.second));
    auto refsmi = MolToSmiles(*refmol);
    auto prodsmi = MolToSmiles(static_cast<const ROMol &>(*rwmol));

    TEST_ASSERT(prodsmi == refsmi);
  }
  BOOST_LOG(rdDebugLog) << "Finished" << std::endl;
}

void testOrganometallics() {
  BOOST_LOG(rdDebugLog)
      << "-----------------------\n test organometallic disconnections"
      << std::endl;
  std::string rdbase = getenv("RDBASE");
  std::string test_dir = rdbase + "/Code/GraphMol/MolStandardize/test_data";
  std::vector<std::pair<std::string, std::string>> test_files = {
      {"CPLX_0001.mol",
       R"(O=C([O-])CN(CCN(CC(=O)[O-])CC(=O)[O-])CC(=O)[O-].[Mg+2].[Na+].[Na+])"},
      {"CPLX_0002.mol",
       R"(O=C([O-])CN(CCN(CC(=O)[O-])CC(=O)[O-])CC(=O)[O-].[Fe+2].[Na+].[Na+])"},
      {"MOL_00001.mol",
       R"([Pd].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1)"},
      {"MOL_00002.mol",
       R"(ClCCl.[Cl-].[Cl-].[Fe+2].[Pd+2].c1ccc(P(c2ccccc2)[c-]2cccc2)cc1.c1ccc(P(c2ccccc2)[c-]2cccc2)cc1)"},
      {"MOL_00003.mol",
       R"([Cl-].[Cl-].[Fe+2].[Pd+2].c1ccc(P(c2ccccc2)[c-]2cccc2)cc1.c1ccc(P(c2ccccc2)[c-]2cccc2)cc1)"},
      {"MOL_00004.mol",
       R"(CC(C)c1cc(C(C)C)c(-c2ccccc2P(C2CCCCC2)C2CCCCC2)c(C(C)C)c1.Nc1ccccc1-c1[c-]cccc1.[Cl-].[Pd+2])"},
      {"MOL_00005.mol",
       R"(CCCCP(C12CC3CC(CC(C3)C1)C2)C12CC3CC(CC(C3)C1)C2.CS(=O)(=O)[O-].Nc1ccccc1-c1[c-]cccc1.[Pd+2])"},
      {"MOL_00006.mol",
       R"(CCCCP(C12CC3CC(CC(C3)C1)C2)C12CC3CC(CC(C3)C1)C2.CS(=O)(=O)[O-].Nc1ccccc1-c1[c-]cccc1.[Pd+2])"},
      {"MOL_00007.mol",
       R"(CC(C)c1cc(C(C)C)c(-c2ccccc2P(C2CCCCC2)C2CCCCC2)c(C(C)C)c1.CS(=O)(=O)[O-].Nc1ccccc1-c1[c-]cccc1.[Pd+2])"},
      {"MOL_00008.mol",
       R"(CC(C)(C)P([c-]1cccc1)C(C)(C)C.CC(C)(C)P([c-]1cccc1)C(C)(C)C.[Cl-].[Cl-].[Fe+2].[Pd+2])"},
      {"MOL_00009.mol",
       R"(CN(C)c1ccc(P(C(C)(C)C)C(C)(C)C)cc1.CN(C)c1ccc(P(C(C)(C)C)C(C)(C)C)cc1.[Cl-].[Cl-].[Pd+2])"},
      {"MOL_00010.mol",
       R"(Cc1cc(C(c2ccccc2)c2ccccc2)c(-n2[c-][n+](-c3c(C(c4ccccc4)c4ccccc4)cc(C)cc3C(c3ccccc3)c3ccccc3)cc2)c(C(c2ccccc2)c2ccccc2)c1.[CH2-]/C=C/c1ccccc1.[Cl-].[Pd+2])"},
      {"MOL_00011.mol",
       R"([Cl-].[Cl-].[Ni+2].c1ccc(P(CCCP(c2ccccc2)c2ccccc2)c2ccccc2)cc1)"},
      {"MOL_00012.mol", R"(C1=C\CC/C=C\CC/1.CC1=C(C)C(=O)C(C)=C(C)C1=O.[Ni])"},
      {"MOL_00013.mol", R"(CN(C)CCN(C)C.Cc1[c-]cccc1.[Cl-].[Ni+2])"},
      {"MOL_00014.mol",
       R"(Cc1cc(C)cc(/C=C/c2cc(C)cc(C)c2)c1.Cc1cc(C)cc(/C=C/c2cc(C)cc(C)c2)c1.Cc1cc(C)cc(/C=C/c2cc(C)cc(C)c2)c1.[Ni])"},
      {"MOL_00015.mol",
       R"(CC(C)(C)c1ccc(/C=C/c2ccc(C(C)(C)C)cc2)cc1.CC(C)(C)c1ccc(/C=C/c2ccc(C(C)(C)C)cc2)cc1.CC(C)(C)c1ccc(/C=C/c2ccc(C(C)(C)C)cc2)cc1.[Ni])"},
      {"diethylzinc.mol", R"([CH2-]CC.[CH2-]CC.[Zn+2])"},
      {"edta_case.mol",
       R"(O=C([O-])CN(CCN(CC(=O)[O-])CC(=O)[O-])CC(=O)[O-].[Mg+2].[Na+].[Na+])"},
      {"ferrocene.mol", R"([Fe+2].c1cc[cH-]c1.c1cc[cH-]c1)"},
      {"ferrocene_1.mol", R"(C[c-]1cccc1.C[c-]1cccc1.[Fe+2])"},
      {"grignard_1.mol", R"([CH3-].[Cl-].[Mg+2])"},
      {"grignard_2.mol", R"([Cl-].[Mg+2].[c-]1ccccc1)"},
      {"lithium_1.mol", R"([CH3-].[Li+])"},
      {"lithium_2.mol", R"([Li+].[c-]1ccccc1)"},
      {"phenylzincbromide.mol", R"([Br-].[Zn+2].[c-]1ccccc1)"},
      {"ruthenium.mol",
       R"([Cl-].[Cl-].[Cl-].[Cl-].[Ru+2].[Ru+2].c1ccccc1.c1ccccc1)"},
      {"sodium_1.mol", R"([Na+].c1cc[cH-]c1)"},
      {"sodium_2.mol",
       R"(CN(C)CCN(C)CCN(C)C.CN(C)CCN(C)CCN(C)C.[Na+].[Na+].[c-]1ccccc1.[c-]1ccccc1)"},
      {"weirdzinc_1.mol",
       R"([Zn+2].[Zn+2].[c-]1ccccc1.[c-]1ccccc1.[c-]1ccccc1.[c-]1ccccc1)"},
      {"weirdzinc_2.mol",
       R"([Cu+].[Cu+].[Zn+2].[Zn+2].[c-]1ccccc1.[c-]1ccccc1.[c-]1ccccc1.[c-]1ccccc1.[c-]1ccccc1.[c-]1ccccc1)"}};
  for (auto &test_file : test_files) {
    std::string full_file = test_dir + "/" + test_file.first;
    bool takeOwnership = true;
    SDMolSupplier mol_supplier(full_file, takeOwnership);
    std::unique_ptr<ROMol> m(mol_supplier.next());
    TEST_ASSERT(m);
    std::unique_ptr<ROMol> dm(MolStandardize::disconnectOrganometallics(*m));
    //    std::cout << test_file.first << " got : " << MolToSmiles(*dm)
    //              << " expected : " << test_file.second << std::endl;
    TEST_ASSERT(MolToSmiles(*dm) == test_file.second);
    std::unique_ptr<RWMol> em(new RWMol(*m));
    MolStandardize::disconnectOrganometallics(*em);
    TEST_ASSERT(MolToSmiles(*em) == test_file.second);
  }
  {
    auto m("[CH2-](->[K+])c1ccccc1"_smiles);
    TEST_ASSERT(m);
    MolStandardize::disconnectOrganometallics(*m);
    TEST_ASSERT(MolToSmiles(*m) == "[CH2-]c1ccccc1.[K+]");
  }
}

int main() {
  RDLog::InitLogs();
  boost::logging::disable_logs("rdApp.info");
  testCleanup();
  testStandardizeSm();
  testMetalDisconnector();
  testNormalize();
  testNormalizeMultiFrags();
  testCharge();
  testMetalDisconnectorLigandExpo();
  testEnumerateTautomerSmiles();
  testOrganometallics();
  return 0;
}
