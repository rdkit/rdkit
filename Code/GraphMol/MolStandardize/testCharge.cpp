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
#include <GraphMol/MolStandardize/AcidBaseCatalog/AcidBaseCatalogParams.h>
#include <GraphMol/MolStandardize/AcidBaseCatalog/AcidBaseCatalogUtils.h>
#include "Charge.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;
using namespace MolStandardize;

void testReionizer() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n test reionizer"
                       << std::endl;

  std::string smi1, smi2, smi3, smi4, smi5, smi6, smi7;

  Reionizer reionizer;

  // Test table salt.
  smi1 = "[Na].[Cl]";
  std::shared_ptr<ROMol> m1(SmilesToMol(smi1));
  ROMOL_SPTR reionized(reionizer.reionize(*m1));
  TEST_ASSERT(MolToSmiles(*reionized) == "[Cl-].[Na+]");

  // Test forced charge correction maintaining overall neutral charge.
  smi2 = "[Na].O=C(O)c1ccccc1";
  std::shared_ptr<ROMol> m2(SmilesToMol(smi2));
  ROMOL_SPTR reionized2(reionizer.reionize(*m2));
  TEST_ASSERT(MolToSmiles(*reionized2) == "O=C([O-])c1ccccc1.[Na+]");

  // Test reionizer moves proton to weaker acid.
  smi3 = "C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O";
  std::shared_ptr<ROMol> m3(SmilesToMol(smi3));
  ROMOL_SPTR reionized3(reionizer.reionize(*m3));
  TEST_ASSERT(MolToSmiles(*reionized3) == "O=S(O)c1ccc(S(=O)(=O)[O-])cc1");

  // Test reionizer moves proton to weaker acid.
  smi5 = "C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O";
  std::shared_ptr<ROMol> m5(SmilesToMol(smi5));
  ROMOL_SPTR reionized5(reionizer.reionize(*m5));
  TEST_ASSERT(MolToSmiles(*reionized3) == "O=S(O)c1ccc(S(=O)(=O)[O-])cc1");

  // Test charged carbon doesn't get recognised as alpha-carbon-hydrogen-keto.
  smi6 = "CCOC(=O)C(=O)[CH-]C#N";
  std::shared_ptr<ROMol> m6(SmilesToMol(smi6));
  ROMOL_SPTR reionized6(reionizer.reionize(*m6));
  TEST_ASSERT(MolToSmiles(*reionized6) == "CCOC(=O)C(=O)[CH-]C#N");

  // TODO... can't make this work. Python SanitizeMol looks to correct...
  // what is different with MolOps::sanitizeMol?
  smi7 = "C[N+]1=C[CH-]N(C(=N)N)/C1=C/[N+](=O)[O-]";
  std::shared_ptr<ROMol> m7(SmilesToMol(smi7));
  ROMOL_SPTR reionized7(reionizer.reionize(*m7));
  TEST_ASSERT(MolToSmiles(*reionized7) ==
              "C[N+]1=CCN(C(=N)N)/C1=[C-]/[N+](=O)[O-]");
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testChargeParent() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n test charge parent"
                       << std::endl;
  std::string smi1, smi2, smi3, smi4, smi5, smi6, smi7, smi8, smi9, smi10,
      smi11, smi12;
  MolStandardize::CleanupParameters params;
  // initialize CleanupParameters with preferOrganic=true
  MolStandardize::CleanupParameters params_preferorg;
  params_preferorg.preferOrganic = true;

  // Test neutralization of ionized acids and bases.
  smi1 = "C(C(=O)[O-])(Cc1n[n-]nn1)(C[NH3+])(C[N+](=O)[O-])";
  std::unique_ptr<RWMol> m1(SmilesToMol(smi1));
  std::unique_ptr<RWMol> res1(MolStandardize::chargeParent(*m1, params));
  TEST_ASSERT(MolToSmiles(*res1) == "NCC(Cc1nn[nH]n1)(C[N+](=O)[O-])C(=O)O");

  // Test preservation of zwitterion.
  smi2 = "n(C)1cc[n+]2cccc([O-])c12";
  std::unique_ptr<RWMol> m2(SmilesToMol(smi2));
  std::unique_ptr<RWMol> res2(MolStandardize::chargeParent(*m2, params));
  TEST_ASSERT(MolToSmiles(*res2) == "Cn1cc[n+]2cccc([O-])c12");

  // Choline should be left with a positive charge.
  smi3 = "C[N+](C)(C)CCO";
  std::unique_ptr<RWMol> m3(SmilesToMol(smi3));
  std::unique_ptr<RWMol> res3(MolStandardize::chargeParent(*m3, params));
  TEST_ASSERT(MolToSmiles(*res3) == "C[N+](C)(C)CCO");

  // Hydrogen should be removed to give deanol as a charge parent.
  smi4 = "C[NH+](C)CCO";
  std::unique_ptr<RWMol> m4(SmilesToMol(smi4));
  std::unique_ptr<RWMol> res4(MolStandardize::chargeParent(*m4, params));
  TEST_ASSERT(MolToSmiles(*res4) == "CN(C)CCO");

  // Sodium benzoate to benzoic acid.
  smi5 = "[Na+].O=C([O-])c1ccccc1";
  std::unique_ptr<RWMol> m5(SmilesToMol(smi5));
  std::unique_ptr<RWMol> res5(MolStandardize::chargeParent(*m5, params));
  TEST_ASSERT(MolToSmiles(*res5) == "O=C(O)c1ccccc1");

  // Benzoate ion to benzoic acid.
  smi6 = "O=C([O-])c1ccccc1";
  std::unique_ptr<RWMol> m6(SmilesToMol(smi6));
  std::unique_ptr<RWMol> res6(MolStandardize::chargeParent(*m6, params));
  TEST_ASSERT(MolToSmiles(*res6) == "O=C(O)c1ccccc1");

  // Charges in histidine should be neutralized.
  smi7 = "[NH3+]C(Cc1cnc[nH]1)C(=O)[O-]";
  std::unique_ptr<RWMol> m7(SmilesToMol(smi7));
  std::unique_ptr<RWMol> res7(MolStandardize::chargeParent(*m7, params));
  TEST_ASSERT(MolToSmiles(*res7) == "NC(Cc1cnc[nH]1)C(=O)O");

  //
  smi8 = "C[NH+](C)(C).[Cl-]";
  std::unique_ptr<RWMol> m8(SmilesToMol(smi8));
  std::unique_ptr<RWMol> res8(MolStandardize::chargeParent(*m8, params));
  TEST_ASSERT(MolToSmiles(*res8) == "CN(C)C");

  // No organic fragments.
  smi9 = "[N+](=O)([O-])[O-]";
  std::unique_ptr<RWMol> m9(SmilesToMol(smi9));
  std::unique_ptr<RWMol> res9(MolStandardize::chargeParent(*m9, params));
  TEST_ASSERT(MolToSmiles(*res9) == "O=[N+]([O-])[O-]");

  // TODO switch prefer_organic=true
  // No organic fragments.
  smi10 = "[N+](=O)([O-])[O-]";
  std::unique_ptr<RWMol> m10(SmilesToMol(smi10));
  std::unique_ptr<RWMol> res10(
      MolStandardize::chargeParent(*m10, params_preferorg));
  TEST_ASSERT(MolToSmiles(*res10) == "O=[N+]([O-])[O-]");

  // Larger inorganic fragment should be chosen.
  smi11 = "[N+](=O)([O-])[O-].[CH2]";
  std::unique_ptr<RWMol> m11(SmilesToMol(smi11));
  std::unique_ptr<RWMol> res11(MolStandardize::chargeParent(*m11, params));
  TEST_ASSERT(MolToSmiles(*res11) == "O=[N+]([O-])[O-]");

  // TODO prefer_organic=true
  // Smaller organic fragment should be chosen over larger inorganic fragment.
  smi12 = "[N+](=O)([O-])[O-].[CH2]";
  std::unique_ptr<RWMol> m12(SmilesToMol(smi12));
  std::unique_ptr<RWMol> res12(
      MolStandardize::chargeParent(*m12, params_preferorg));
  TEST_ASSERT(MolToSmiles(*res12) == "[CH2]");
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGithub2144() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing github #2144: "
                          "Error when calling ChargeParent twice"
                       << std::endl;

  {
    // Test neutralization of ionized acids and bases.
    auto m1 = "c1ccccn1"_smiles;
    TEST_ASSERT(m1);
    std::unique_ptr<RWMol> res1(MolStandardize::chargeParent(*m1));
    TEST_ASSERT(res1);
    TEST_ASSERT(MolToSmiles(*res1) == MolToSmiles(*m1));

    std::unique_ptr<RWMol> res2(MolStandardize::chargeParent(*res1));
    TEST_ASSERT(res2);
    TEST_ASSERT(MolToSmiles(*res2) == MolToSmiles(*m1));
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

int main() {
  RDLog::InitLogs();
  testReionizer();
  testChargeParent();
  testGithub2144();
  return 0;
}
