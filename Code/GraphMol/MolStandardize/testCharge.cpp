//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include "MolStandardize.h"
#include <GraphMol/MolStandardize/AcidBaseCatalog/AcidBaseCatalogParams.h>
#include <GraphMol/MolStandardize/AcidBaseCatalog/AcidBaseCatalogUtils.h>
#include "Charge.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>

using namespace RDKit;
using namespace MolStandardize;

void testReionizer() {
  BOOST_LOG(rdDebugLog) << "-----------------------\n test reionizer"
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
  BOOST_LOG(rdDebugLog) << "Finished" << std::endl;
}

void testChargeParent() {
  BOOST_LOG(rdDebugLog) << "-----------------------\n test charge parent"
                        << std::endl;
  MolStandardize::CleanupParameters params;
  // initialize CleanupParameters with preferOrganic=true
  MolStandardize::CleanupParameters params_preferorg;
  params_preferorg.preferOrganic = true;

  // Test neutralization of ionized acids and bases.
  auto m1 = "C(C(=O)[O-])(Cc1n[n-]nn1)(C[NH3+])(C[N+](=O)[O-])"_smiles;
  std::unique_ptr<RWMol> res1(MolStandardize::chargeParent(*m1, params));
  TEST_ASSERT(MolToSmiles(*res1) == "NCC(Cc1nn[nH]n1)(C[N+](=O)[O-])C(=O)O");

  // Test preservation of zwitterion.
  auto m2 = "n(C)1cc[n+]2cccc([O-])c12"_smiles;
  std::unique_ptr<RWMol> res2(MolStandardize::chargeParent(*m2, params));
  TEST_ASSERT(MolToSmiles(*res2) == "Cn1cc[n+]2cccc([O-])c12");

  // Choline should be left with a positive charge.
  auto m3 = "C[N+](C)(C)CCO"_smiles;
  std::unique_ptr<RWMol> res3(MolStandardize::chargeParent(*m3, params));
  TEST_ASSERT(MolToSmiles(*res3) == "C[N+](C)(C)CCO");

  // Hydrogen should be removed to give deanol as a charge parent.
  auto m4 = "C[NH+](C)CCO"_smiles;
  std::unique_ptr<RWMol> res4(MolStandardize::chargeParent(*m4, params));
  TEST_ASSERT(MolToSmiles(*res4) == "CN(C)CCO");

  // Sodium benzoate to benzoic acid.
  auto m5 = "[Na+].O=C([O-])c1ccccc1"_smiles;
  std::unique_ptr<RWMol> res5(MolStandardize::chargeParent(*m5, params));
  TEST_ASSERT(MolToSmiles(*res5) == "O=C(O)c1ccccc1");

  // Benzoate ion to benzoic acid.
  auto m6 = "O=C([O-])c1ccccc1"_smiles;
  std::unique_ptr<RWMol> res6(MolStandardize::chargeParent(*m6, params));
  TEST_ASSERT(MolToSmiles(*res6) == "O=C(O)c1ccccc1");

  // Charges in histidine should be neutralized.
  auto m7 = "[NH3+]C(Cc1cnc[nH]1)C(=O)[O-]"_smiles;
  std::unique_ptr<RWMol> res7(MolStandardize::chargeParent(*m7, params));
  TEST_ASSERT(MolToSmiles(*res7) == "NC(Cc1cnc[nH]1)C(=O)O");

  //
  auto m8 = "C[NH+](C)(C).[Cl-]"_smiles;
  std::unique_ptr<RWMol> res8(MolStandardize::chargeParent(*m8, params));
  TEST_ASSERT(MolToSmiles(*res8) == "CN(C)C");

  // No organic fragments.
  auto m9 = "[N+](=O)([O-])[O-]"_smiles;
  std::unique_ptr<RWMol> res9(MolStandardize::chargeParent(*m9, params));
  TEST_ASSERT(MolToSmiles(*res9) == "O=[N+]([O-])O");

  // TODO switch prefer_organic=true
  // No organic fragments.
  auto m10 = "[N+](=O)([O-])[O-]"_smiles;
  std::unique_ptr<RWMol> res10(
      MolStandardize::chargeParent(*m10, params_preferorg));
  TEST_ASSERT(MolToSmiles(*res10) == "O=[N+]([O-])O");

  // Larger inorganic fragment should be chosen.
  auto m11 = "[N+](=O)([O-])[O-].[CH2]"_smiles;
  std::unique_ptr<RWMol> res11(MolStandardize::chargeParent(*m11, params));
  TEST_ASSERT(MolToSmiles(*res11) == "O=[N+]([O-])O");

  // TODO prefer_organic=true
  // Smaller organic fragment should be chosen over larger inorganic fragment.
  auto m12 = "[N+](=O)([O-])[O-].[CH2]"_smiles;
  std::unique_ptr<RWMol> res12(
      MolStandardize::chargeParent(*m12, params_preferorg));
  TEST_ASSERT(MolToSmiles(*res12) == "[CH2]");

  // do not completely neutralize zwitterions
  auto m13 = "C[S+](=O)([O-])NC"_smiles;
  std::unique_ptr<RWMol> res13(MolStandardize::chargeParent(*m13, params));
  TEST_ASSERT(MolToSmiles(*res13) == "CN[S+](C)(=O)[O-]");

  // standalone metal ion
  auto m14 = "[Cu+2]"_smiles;
  std::unique_ptr<RWMol> res14(MolStandardize::chargeParent(*m14));
  TEST_ASSERT(MolToSmiles(*res14) == "[Cu+2]");

  BOOST_LOG(rdDebugLog) << "Finished" << std::endl;
}

void testGithub2144() {
  BOOST_LOG(rdDebugLog) << "-----------------------\n Testing github #2144: "
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
  BOOST_LOG(rdDebugLog) << "Finished" << std::endl;
}

void testGithub2346() {
  BOOST_LOG(rdDebugLog) << "-----------------------\n Testing github #2346: "
                           "uncharger behaves differently on molecules "
                           "constructed from mol blocks and SMILES"
                        << std::endl;

  {
    auto m1 = "[NH3+]CC[O-]"_smiles;
    TEST_ASSERT(m1);
    MolStandardize::Uncharger uncharger;

    std::unique_ptr<ROMol> res1(uncharger.uncharge(*m1));
    TEST_ASSERT(res1);
    TEST_ASSERT(res1->getAtomWithIdx(0)->getFormalCharge() == 0);
    TEST_ASSERT(res1->getAtomWithIdx(1)->getFormalCharge() == 0);

    std::unique_ptr<ROMol> m2(MolBlockToMol(MolToMolBlock(*m1)));
    TEST_ASSERT(m2);
    std::unique_ptr<ROMol> res2(uncharger.uncharge(*m2));
    TEST_ASSERT(res2);
    TEST_ASSERT(res2->getAtomWithIdx(0)->getFormalCharge() == 0);
    TEST_ASSERT(res2->getAtomWithIdx(1)->getFormalCharge() == 0);
  }
  {
    auto m1 = "[O-]C(=O)C([O-])C(=O)[O-]"_smiles;
    TEST_ASSERT(m1);
    MolStandardize::Uncharger uncharger;

    std::unique_ptr<ROMol> res1(uncharger.uncharge(*m1));
    TEST_ASSERT(res1);
    for (auto &atom : res1->atoms()) {
      TEST_ASSERT(atom->getFormalCharge() == 0);
    }

    std::unique_ptr<ROMol> m2(MolBlockToMol(MolToMolBlock(*m1)));
    TEST_ASSERT(m2);
    std::unique_ptr<ROMol> res2(uncharger.uncharge(*m2));
    TEST_ASSERT(res2);
    for (auto &atom : res2->atoms()) {
      TEST_ASSERT(atom->getFormalCharge() == 0);
    }
  }

  BOOST_LOG(rdDebugLog) << "Finished" << std::endl;
}

void testChargedAromatics() {
  BOOST_LOG(rdDebugLog)
      << "-----------------------\n Testing charged aromatics: "
         "need to sanitize after using uncharger"
      << std::endl;

  {
    auto cyclopentadienyl = "[cH-]1cccc1"_smiles;
    TEST_ASSERT(cyclopentadienyl);
    MolStandardize::Uncharger uncharger;

    std::unique_ptr<ROMol> res(uncharger.uncharge(*cyclopentadienyl));
    TEST_ASSERT(res.get());
    TEST_ASSERT(MolToSmiles(*res) == "c1cccc1");
    MolOps::sanitizeMol(*static_cast<RWMol *>(res.get()));
    TEST_ASSERT(MolToSmiles(*res) == "C1=CCC=C1");
  }
  {
    auto tropylium = "[cH+]1cccccc1"_smiles;
    TEST_ASSERT(tropylium);
    MolStandardize::Uncharger uncharger;

    std::unique_ptr<ROMol> res(uncharger.uncharge(*tropylium));
    TEST_ASSERT(res.get());
    TEST_ASSERT(MolToSmiles(*res) == "c1cccccc1");
    MolOps::sanitizeMol(*static_cast<RWMol *>(res.get()));
    TEST_ASSERT(MolToSmiles(*res) == "C1=CC=CCC=C1");
  }
  {
    auto azolium = "[NH2+]1C=CC=C1"_smiles;
    TEST_ASSERT(azolium);
    MolStandardize::Uncharger uncharger;

    std::unique_ptr<ROMol> res(uncharger.uncharge(*azolium));
    TEST_ASSERT(res.get());
    TEST_ASSERT(MolToSmiles(*res) == "C1=CNC=C1");
    MolOps::sanitizeMol(*static_cast<RWMol *>(res.get()));
    TEST_ASSERT(MolToSmiles(*res) == "c1cc[nH]c1");
  }

  BOOST_LOG(rdDebugLog) << "Finished" << std::endl;
}

void testUnchargerProtonationOnly() {
  BOOST_LOG(rdDebugLog)
      << "-----------------------\n Testing removal of formal charges limited to "
         "changing the protonation state (disable the addition/removal of H-)"
      << std::endl;
  {
    // Uncharger options
    bool canonicalOrdering {false};
    bool force;
    bool protonationOnly {true};
    // simple test verifying that for protic compounds the behavior
    // doesn't change if the protonationOnly option is set
    auto m1 = "C[N+](C)(C)CC[O-]"_smiles;
    TEST_ASSERT(m1);
    // with force=false the zwitterion should stay unmodified
    force = false;
    MolStandardize::Uncharger uncharger1(canonicalOrdering, force, protonationOnly);
    std::unique_ptr<ROMol> res1(uncharger1.uncharge(*m1));
    TEST_ASSERT(res1.get());
    TEST_ASSERT(MolToSmiles(*res1) == "C[N+](C)(C)CC[O-]");
    // with force=true the oxygen should be neutralized
    force = true;
    MolStandardize::Uncharger uncharger2(canonicalOrdering, force, protonationOnly);
    std::unique_ptr<ROMol> res2(uncharger2.uncharge(*m1));
    TEST_ASSERT(res2.get());
    TEST_ASSERT(MolToSmiles(*res2) == "C[N+](C)(C)CCO");
  }
  {
    // Uncharger options
    bool canonicalOrdering {false};
    bool force {true};
    bool protonationOnly {true};

    auto tropylium = "[cH+]1cccccc1"_smiles;
    TEST_ASSERT(tropylium);
    // try uncharging as much as possible, but only allow
    // protonating/deprotonating.
    MolStandardize::Uncharger uncharger(canonicalOrdering, force, protonationOnly);
    // tropylium should stay unmodified
    std::unique_ptr<ROMol> res(uncharger.uncharge(*tropylium));
    TEST_ASSERT(res.get());
    TEST_ASSERT(MolToSmiles(*res) == "c1ccc[cH+]cc1");
  }
  {
    // Uncharger options
    bool canonicalOrdering {false};
    bool force {true};
    bool protonationOnly {true};

    auto boronhydride = "[BH4-]"_smiles;
    TEST_ASSERT(boronhydride);
    // try uncharging as much as possible, but only allow
    // protonating/deprotonating.
    MolStandardize::Uncharger uncharger(canonicalOrdering, force, protonationOnly);
    // boronhydride should stay unmodified
    std::unique_ptr<ROMol> res(uncharger.uncharge(*boronhydride));
    TEST_ASSERT(res.get());
    TEST_ASSERT(MolToSmiles(*res) == "[BH4-]");
  }
  {
    // Uncharger options
    bool canonicalOrdering {false};
    bool force;
    bool protonationOnly {true};

    // Test the neutralization of a zwitterion, where the positive charge is a
    // carbocation and not possible to remove when protonationOnly is enabled
    auto m = "[CH2+]c1ccc([O-])cc1"_smiles;
    TEST_ASSERT(m);

    // with force=false the zwitterion should stay unmodified
    force = false;
    MolStandardize::Uncharger uncharger1(canonicalOrdering, force, protonationOnly);
    std::unique_ptr<ROMol> res1(uncharger1.uncharge(*m));
    TEST_ASSERT(res1.get());
    TEST_ASSERT(MolToSmiles(*res1) == "[CH2+]c1ccc([O-])cc1");
    // with force=true the oxygen should be neutralized
    force = true;
    MolStandardize::Uncharger uncharger2(canonicalOrdering, force, protonationOnly);
    std::unique_ptr<ROMol> res2(uncharger2.uncharge(*m));
    TEST_ASSERT(res2.get());
    TEST_ASSERT(MolToSmiles(*res2) == "[CH2+]c1ccc(O)cc1");
  }
  BOOST_LOG(rdDebugLog) << "Finished" << std::endl;
}

void testInorganicAcids() {
  BOOST_LOG(rdDebugLog) << "-----------------------\n Testing inorganic acids"
                        << std::endl;
  MolStandardize::Uncharger uncharger;
  std::vector<std::string> halogens{"Cl", "Br", "I"};
  std::unique_ptr<ROMol> res;
  for (const auto &halogen : halogens) {
    std::unique_ptr<ROMol> hypohalite(SmilesToMol("[" + halogen + "][O-]"));
    TEST_ASSERT(hypohalite);
    res.reset(uncharger.uncharge(*hypohalite));
    TEST_ASSERT(MolToSmiles(*res) == "O" + halogen);
    std::unique_ptr<ROMol> halite(SmilesToMol("[" + halogen + "](=O)[O-]"));
    TEST_ASSERT(halite);
    res.reset(uncharger.uncharge(*halite));
    TEST_ASSERT(MolToSmiles(*res) == "[O-][" + halogen + "+]O");
    std::unique_ptr<ROMol> halate(SmilesToMol("[" + halogen + "](=O)(=O)[O-]"));
    TEST_ASSERT(halate);
    res.reset(uncharger.uncharge(*halate));
    TEST_ASSERT(MolToSmiles(*res) == "[O-][" + halogen + "+2]([O-])O");
    std::unique_ptr<ROMol> perhalate(
        SmilesToMol("[" + halogen + "](=O)(=O)(=O)[O-]"));
    TEST_ASSERT(perhalate);
    res.reset(uncharger.uncharge(*perhalate));
    TEST_ASSERT(MolToSmiles(*res) == "[O-][" + halogen + "+3]([O-])([O-])O");
    // also test uncharging the already neutralized acid
    std::unique_ptr<ROMol> perhalic_acid(
        SmilesToMol("[" + halogen + "](=O)(=O)(=O)O"));
    TEST_ASSERT(perhalic_acid);
    res.reset(uncharger.uncharge(*perhalic_acid));
    TEST_ASSERT(MolToSmiles(*res) == "[O-][" + halogen + "+3]([O-])([O-])O");
  }
  {
    auto hyponitrite = "[O-]N=N[O-]"_smiles;
    TEST_ASSERT(hyponitrite);
    res.reset(uncharger.uncharge(*hyponitrite));
    TEST_ASSERT(MolToSmiles(*res) == "ON=NO");
  }
  {
    auto nitrite = "N(=O)[O-]"_smiles;
    TEST_ASSERT(nitrite);
    res.reset(uncharger.uncharge(*nitrite));
    TEST_ASSERT(MolToSmiles(*res) == "O=NO");
  }
  {
    auto nitrate = "N(=O)(=O)[O-]"_smiles;
    TEST_ASSERT(nitrate);
    res.reset(uncharger.uncharge(*nitrate));
    TEST_ASSERT(MolToSmiles(*res) == "O=[N+]([O-])O");
  }
  {
    auto hyposulfite = "S([O-])[O-]"_smiles;
    TEST_ASSERT(hyposulfite);
    res.reset(uncharger.uncharge(*hyposulfite));
    TEST_ASSERT(MolToSmiles(*res) == "OSO");
  }
  {
    auto sulfite = "S(=O)([O-])[O-]"_smiles;
    TEST_ASSERT(sulfite);
    res.reset(uncharger.uncharge(*sulfite));
    TEST_ASSERT(MolToSmiles(*res) == "O=S(O)O");
  }
  {
    auto sulfate = "S(=O)(=O)([O-])[O-]"_smiles;
    TEST_ASSERT(sulfate);
    res.reset(uncharger.uncharge(*sulfate));
    TEST_ASSERT(MolToSmiles(*res) == "O=S(=O)(O)O");
  }
  {
    auto persulfate = "S(=O)(=O)([O-])OOS(=O)(=O)[O-]"_smiles;
    TEST_ASSERT(persulfate);
    res.reset(uncharger.uncharge(*persulfate));
    TEST_ASSERT(MolToSmiles(*res) == "O=S(=O)(O)OOS(=O)(=O)O");
  }
  {
    auto hypophosphite = "P(=O)[O-]"_smiles;
    TEST_ASSERT(hypophosphite);
    res.reset(uncharger.uncharge(*hypophosphite));
    TEST_ASSERT(MolToSmiles(*res) == "O=PO");
  }
  {
    auto phosphite = "P(=O)([O-])[O-]"_smiles;
    TEST_ASSERT(phosphite);
    res.reset(uncharger.uncharge(*phosphite));
    TEST_ASSERT(MolToSmiles(*res) == "O=[PH](O)O");
  }
  {
    auto phosphate = "P(=O)([O-])([O-])[O-]"_smiles;
    TEST_ASSERT(phosphate);
    res.reset(uncharger.uncharge(*phosphate));
    TEST_ASSERT(MolToSmiles(*res) == "O=P(O)(O)O");
  }
  BOOST_LOG(rdDebugLog) << "Finished" << std::endl;
}

void testReionizerParams() {
  BOOST_LOG(rdDebugLog)
      << "-----------------------\n Testing reionizer parameters" << std::endl;
  {  // defaults
    Reionizer reionizer;
    auto m1 = "c1cc([O-])cc(C(=O)O)c1"_smiles;
    std::unique_ptr<ROMol> reionized1{reionizer.reionize(*m1)};
    TEST_ASSERT(MolToSmiles(*reionized1) == "O=C([O-])c1cccc(O)c1");

    auto m2 = "C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O"_smiles;
    std::unique_ptr<ROMol> reionized2{reionizer.reionize(*m2)};
    TEST_ASSERT(MolToSmiles(*reionized2) == "O=S(O)c1ccc(S(=O)(=O)[O-])cc1");
  }
  {  // parameters via tuple
    std::vector<std::tuple<std::string, std::string, std::string>> params{
        {"-CO2H", "C(=O)[OH]", "C(=O)[O-]"}, {"phenol", "c[OH]", "c[O-]"}};
    Reionizer reionizer(params);
    auto m1 = "c1cc([O-])cc(C(=O)O)c1"_smiles;
    std::unique_ptr<ROMol> reionized1{reionizer.reionize(*m1)};
    TEST_ASSERT(MolToSmiles(*reionized1) == "O=C([O-])c1cccc(O)c1");

    auto m2 = "C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O"_smiles;
    std::unique_ptr<ROMol> reionized2{reionizer.reionize(*m2)};
    TEST_ASSERT(MolToSmiles(*reionized2) == "O=S([O-])c1ccc(S(=O)(=O)O)cc1");
  }
  BOOST_LOG(rdDebugLog) << "Finished" << std::endl;
}
int main() {
  RDLog::InitLogs();
  boost::logging::disable_logs("rdApp.info");
  testReionizer();
  testChargeParent();
  testGithub2144();
  testGithub2346();
  testChargedAromatics();
  testUnchargerProtonationOnly();
  testInorganicAcids();
  testReionizerParams();
  return 0;
}
