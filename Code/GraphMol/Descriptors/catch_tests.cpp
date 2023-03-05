//
//  Copyright (C) 2021 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Descriptors/ConnectivityDescriptors.h>
#include <GraphMol/Descriptors/PMI.h>

using namespace RDKit;

TEST_CASE("Kier kappa2", "[2D]") {
  SECTION("values from https://doi.org/10.1002/qsar.19860050103") {
    std::vector<std::pair<std::string, double>> data = {
      // Table 5 from the paper
      {"c1ccccc15.Cl5", 1.987},
      {"c1ccccc15.F5", 1.735},
      {"c1ccccc15.[N+]5(=O)O", 2.259},
      {"c1ccccc15.C5(=O)C", 2.444},
#if 0
      // expected values from paper (differences are due to hybridization mismatches)
        {"c1ccccc15.N5(C)C", 2.646},
        {"c1ccccc15.C5(=O)N", 2.416},
        {"c1ccccc15.C5(=O)O", 2.416},
        {"c1ccccc15.S5(=O)(=O)C", 2.617},
        {"c1ccccc15.O5", 1.756},
#else
      {"c1ccccc15.N5(C)C", 2.53},
      {"c1ccccc15.C5(=O)N", 2.31},
      {"c1ccccc15.C5(=O)O", 2.31},
      {"c1ccccc15.S5(=O)(=O)C", 2.42},
      {"c1ccccc15.O5", 1.65},
#endif
    };
    for (const auto &pr : data) {
      std::unique_ptr<ROMol> m(SmilesToMol(pr.first));
      REQUIRE(m);
      auto k2 = Descriptors::calcKappa2(*m);
      CHECK(k2 == Approx(pr.second).epsilon(0.01));
    }
  }
}

TEST_CASE("Kier Phi", "[2D]") {
  SECTION("regression-test values from the paper") {
    std::vector<std::pair<std::string, double>> data = {
      // Table 1 from the paper
      {"CCCCCC", 5.00},
      {"CCCCCCC", 6.00},
      {"CCCCCCCC", 7.00},
      {"CCC(C)CC", 3.20},
      {"CC(C)C(C)C", 2.22},
      {"CC(C)(C)CC", 1.63},
      {"C1CCCC1", 0.92},
      {"C1CCCCC1", 1.54},
      {"C1CCCCCC1", 2.25},
      {"CCCCC=C", 4.53},
      {"C=CCCC=C", 4.07},
      {"C#CCCCC", 4.21},
      {"c1ccccc1", 0.91},
      // additional from Table 2
      {"C=CCC=CC", 4.09},
      {"CC=CC=CC", 4.09},
      {"C1=CCCCC1", 1.31},
      {"C1=CC=CCC1", 1.1},
      {"C1=CCC=CC1", 1.1},
      // additional from Table 3
      {"CCCCCCCCCC", 9.00},
      {"CC(C)CCCCC", 5.14},
      {"CCC(C)CCCC", 5.14},
      {"CC(C)CCCC", 4.17},
      {"CCC(C)CCC", 4.17},
      {"CCC(CC)CCC", 5.14},
      {"CCC(CC)CC", 4.17},
      {"CC(C)(C)CCC", 2.34},
      {"CC(C)C(C)CC", 3.06},
      {"CCC(C)(C)CC", 2.34},
      {"CC(C)(C)C(C)C", 1.85},
      // additional from table 4
      {"CCOCC", 3.93},
      {"CCC(=O)CC", 2.73},
      {"CCc1ccc(CC)cc1", 2.49},
      {"CCC(O)CC", 3.14},
      {"CCCC(Cl)(Cl)CCC", 4.69},
      {"CCC(F)C(F)CC", 3.75},
#if 0
      // expected values from paper (differences are due to hybridization mismatches)
        {"CCOC(=O)CC", 3.61},
        {"CCC(=O)Nc1ccc(CC)cc1", 3.65},
#else
      {"CCOC(=O)CC", 3.38},
      {"CCC(=O)Nc1ccc(CC)cc1", 3.50},
#endif

    };
    for (const auto &pr : data) {
      std::unique_ptr<ROMol> m(SmilesToMol(pr.first));
      REQUIRE(m);
      auto val = Descriptors::calcPhi(*m);
      CHECK(val == Approx(pr.second).epsilon(0.01));
    }
  }
}

TEST_CASE("atom counts", "[2D]") {
  {
    SmilesParserParams ps;
    ps.sanitize = true;
    ps.removeHs = false;
    std::unique_ptr<ROMol> m(SmilesToMol("CC[H]", ps));
    REQUIRE(m);
    CHECK(Descriptors::calcNumHeavyAtoms(*m) == 2);
    CHECK(Descriptors::calcNumAtoms(*m) == 8);
  }
}
#ifdef RDK_BUILD_DESCRIPTORS3D
TEST_CASE(
    "Github #4167: SpherocityIndex() not being recalculated for different "
    "conformers",
    "[3D]") {
  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/github4167.sdf";
  bool sanitize = true;
  bool removeHs = false;
  RDKit::SDMolSupplier reader(sdfName, sanitize, removeHs);
  std::unique_ptr<ROMol> m1(reader.next());
  std::unique_ptr<ROMol> m2(reader.next());
  REQUIRE(m1);
  REQUIRE(m2);
  m1->addConformer(new RDKit::Conformer(m2->getConformer()), true);
  REQUIRE(m1->getNumConformers() == 2);

  {
    int confId = 0;
    auto v1_0 = RDKit::Descriptors::spherocityIndex(*m1, confId, true);
    auto v2 = RDKit::Descriptors::spherocityIndex(*m2, confId, true);
    confId = 1;
    auto v1_1 = RDKit::Descriptors::spherocityIndex(*m1, confId, true);
    CHECK(v1_0 != v1_1);
    CHECK(v1_1 == v2);
  }
  {
    std::vector<double (*)(const ROMol &, int, bool, bool)> funcs{
        RDKit::Descriptors::NPR1,
        RDKit::Descriptors::NPR2,
        RDKit::Descriptors::PMI1,
        RDKit::Descriptors::PMI2,
        RDKit::Descriptors::PMI3,
        RDKit::Descriptors::radiusOfGyration,
        RDKit::Descriptors::inertialShapeFactor,
        RDKit::Descriptors::eccentricity,
        RDKit::Descriptors::asphericity};
    for (const auto func : funcs) {
      bool useAtomMasses = true;
      bool force = true;
      int confId = 0;
      auto v1_0 = func(*m1, confId, useAtomMasses, force);
      auto v2 = func(*m2, confId, useAtomMasses, force);
      confId = 1;
      auto v1_1 = func(*m1, confId, useAtomMasses, force);
      CHECK(v1_0 != v1_1);
      CHECK(v1_1 == v2);
    }
  }
  // { // surrogate for NPR1, NPR2, PMI1, PM2,
  //   int confId = 0;
  //   auto v1_0 = RDKit::Descriptors::spherocityIndex(*m1, confId, true);
  //   auto v2 = RDKit::Descriptors::spherocityIndex(*m2, confId, true);
  //   confId = 1;
  //   auto v1_1 = RDKit::Descriptors::spherocityIndex(*m1, confId, true);
  //   CHECK(v1_0 != v1_1);
  //   CHECK(v1_1 == v2);
  // }
}
#endif

TEST_CASE(
    "Github #5104: NumRotatableBonds() incorrect for partially sanitized "
    "molecule") {
  SECTION("basics") {
    auto m1 = "c1ccccc1c1ccc(CCC)cc1"_smiles;
    CHECK(Descriptors::calcNumRotatableBonds(
              *m1, Descriptors::NumRotatableBondsOptions::NonStrict) == 3);
    CHECK(Descriptors::calcNumRotatableBonds(
              *m1, Descriptors::NumRotatableBondsOptions::Strict) == 3);
    CHECK(Descriptors::calcNumRotatableBonds(
              *m1, Descriptors::NumRotatableBondsOptions::StrictLinkages) == 2);
    m1->getBondBetweenAtoms(5, 6)->setBondType(Bond::BondType::AROMATIC);
    CHECK(Descriptors::calcNumRotatableBonds(
              *m1, Descriptors::NumRotatableBondsOptions::NonStrict) == 3);
    CHECK(Descriptors::calcNumRotatableBonds(
              *m1, Descriptors::NumRotatableBondsOptions::Strict) == 3);
    CHECK(Descriptors::calcNumRotatableBonds(
              *m1, Descriptors::NumRotatableBondsOptions::StrictLinkages) == 2);
  }
  SECTION("as reported") {
    SmilesParserParams ps;
    ps.sanitize = false;
    std::unique_ptr<RWMol> m1(SmilesToMol("c1ccccc1c1ccc(CCC)cc1", ps));
    REQUIRE(m1);
    unsigned int whatFailed;
    MolOps::sanitizeMol(*m1, whatFailed,
                        MolOps::SanitizeFlags::SANITIZE_ALL ^
                            MolOps::SanitizeFlags::SANITIZE_KEKULIZE);
    CHECK(Descriptors::calcNumRotatableBonds(
              *m1, Descriptors::NumRotatableBondsOptions::NonStrict) == 3);
    CHECK(Descriptors::calcNumRotatableBonds(
              *m1, Descriptors::NumRotatableBondsOptions::Strict) == 3);
    CHECK(Descriptors::calcNumRotatableBonds(
              *m1, Descriptors::NumRotatableBondsOptions::StrictLinkages) == 2);
  }
  SECTION("linkages") {
    auto m1 = "c1cc[nH]c1c1[nH]c(CCC)cc1"_smiles;
    CHECK(Descriptors::calcNumRotatableBonds(
              *m1, Descriptors::NumRotatableBondsOptions::StrictLinkages) == 3);
    m1->getBondBetweenAtoms(4, 5)->setBondType(Bond::BondType::AROMATIC);
    CHECK(Descriptors::calcNumRotatableBonds(
              *m1, Descriptors::NumRotatableBondsOptions::StrictLinkages) == 3);
  }
}

TEST_CASE("TPSA caching ignores options") {
  SECTION("with force") {
    auto m = "OCCS"_smiles;
    REQUIRE(m);
    bool force = true;
    double tv1 = Descriptors::calcTPSA(*m, force, false);
    double tv2 = Descriptors::calcTPSA(*m, force, true);
    CHECK(tv2 > tv1);
  }
  SECTION("no force") {
    auto m = "OCCS"_smiles;
    REQUIRE(m);
    bool force = false;
    double tv1 = Descriptors::calcTPSA(*m, force, false);
    double tv2 = Descriptors::calcTPSA(*m, force, true);
    CHECK(tv2 > tv1);
  }
}