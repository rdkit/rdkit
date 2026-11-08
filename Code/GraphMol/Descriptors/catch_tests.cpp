//
//  Copyright (C) 2021-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Descriptors/ConnectivityDescriptors.h>
#include <GraphMol/Descriptors/OxidationNumbers.h>
#include <GraphMol/Descriptors/PMI.h>
#include <GraphMol/Descriptors/DCLV.h>
#include <GraphMol/Descriptors/BCUT.h>
#include <GraphMol/Descriptors/Property.h>
#include <cmath>
#include <fstream>
#include <memory>
#include <mutex>
#include <new>
#include <thread>
#include <vector>
#ifdef RDK_BUILD_DESCRIPTORS3D
#include <GraphMol/Descriptors/GETAWAY.h>
#include <GraphMol/Descriptors/EEM.h>
#endif

using namespace RDKit;

TEST_CASE("Kier kappa2", "[2D]") {
  SECTION("values from https://doi.org/10.1002/qsar.19860050103") {
    std::vector<std::pair<std::string, double>> data = {
        // Table 5 from the paper
        {"c1ccccc15.Cl5", 1.987},
        {"c1ccccc15.F5", 1.735},
        {"c1ccccc15.[N+]5(=O)O", 2.259},
        {"c1ccccc15.C5(=O)C", 2.444},
        /* expected values from paper (differences are due to hybridization
           mismatches)
          {"c1ccccc15.N5(C)C", 2.646},
          {"c1ccccc15.C5(=O)N", 2.416},
          {"c1ccccc15.C5(=O)O", 2.416},
          {"c1ccccc15.S5(=O)(=O)C", 2.617},
          {"c1ccccc15.O5", 1.756},
        */
        {"c1ccccc15.N5(C)C", 2.53},
        {"c1ccccc15.C5(=O)N", 2.31},
        {"c1ccccc15.C5(=O)O", 2.31},
        {"c1ccccc15.S5(=O)(=O)C", 2.42},
        {"c1ccccc15.O5", 1.65},
    };
    for (const auto &pr : data) {
      std::unique_ptr<ROMol> m(SmilesToMol(pr.first));
      REQUIRE(m);
      auto k2 = Descriptors::calcKappa2(*m);
      CHECK(k2 == Catch::Approx(pr.second).epsilon(0.01));
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
        /* expected values from paper (differences are due to hybridization
          mismatches)
          {"CCOC(=O)CC", 3.61},
          {"CCC(=O)Nc1ccc(CC)cc1", 3.65},
        */
        {"CCOC(=O)CC", 3.38},
        {"CCC(=O)Nc1ccc(CC)cc1", 3.50},
    };
    for (const auto &pr : data) {
      std::unique_ptr<ROMol> m(SmilesToMol(pr.first));
      REQUIRE(m);
      auto val = Descriptors::calcPhi(*m);
      CHECK(val == Catch::Approx(pr.second).epsilon(0.01));
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
    auto m2 = "CCCC"_smiles;
    MolOps::addHs(*m2);
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
    CHECK(Descriptors::calcNumRotatableBonds(
              *m2, Descriptors::NumRotatableBondsOptions::Strict) == 1);
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

TEST_CASE("Oxidation numbers") {
  std::string rdbase = std::getenv("RDBASE");
  SECTION("simple tests") {
    {
      std::vector<std::string> smis{"CO", "C=O", "C(=O)O", "S(=O)(=O)(O)O"};
      std::vector<std::vector<int>> expected{
          {-2, -2}, {0, -2}, {2, -2, -2}, {6, -2, -2, -2, -2}};
      for (auto i = 0u; i < smis.size(); ++i) {
        std::unique_ptr<RWMol> mol(RDKit::SmilesToMol(smis[i]));
        Descriptors::calcOxidationNumbers(*mol);
        for (const auto &a : mol->atoms()) {
          CHECK(a->getProp<int>(common_properties::OxidationNumber) ==
                expected[i][a->getIdx()]);
        }
      }
    }
  }
  SECTION("organometallics tests") {
    {
      std::string file1 =
          rdbase + "/Code/GraphMol/MolStandardize/test_data/ferrocene.mol";
      std::vector<int> expected{-2, -1, -1, -1, -1, -2, -1,
                                -1, -1, -1, 2,  0,  0};
      bool takeOwnership = true;
      SDMolSupplier mol_supplier(file1, takeOwnership);
      std::unique_ptr<ROMol> m1(mol_supplier.next());
      REQUIRE(m1);
      Descriptors::calcOxidationNumbers(*m1);
      for (const auto &a : m1->atoms()) {
        CHECK(a->getProp<int>(common_properties::OxidationNumber) ==
              expected[a->getIdx()]);
      }
    }
    {
      std::string file2 =
          rdbase + "/Code/GraphMol/MolStandardize/test_data/MOL_00002.mol";
      bool takeOwnership = true;
      SDMolSupplier mol_supplier(file2, takeOwnership);
      std::unique_ptr<ROMol> m1(mol_supplier.next());
      REQUIRE(m1);
      RWMol m2(*m1);
      RDKit::MolOps::Kekulize(m2);
      std::vector<unsigned int> ats{0, 5, 10, 13, 14, 19, 20, 21, 42, 43, 44};
      std::vector<int> expected{-2, -2, 2, 3, 3, 2, -1, -1, -1, 0, -1};
      Descriptors::calcOxidationNumbers(m2);
      for (unsigned int i = 0; i < ats.size(); ++i) {
        auto a = m2.getAtomWithIdx(ats[i]);
        CHECK(a->getProp<int>(common_properties::OxidationNumber) ==
              expected[i]);
      }
    }
    {
      std::string file3 =
          rdbase + "/Code/GraphMol/MolStandardize/test_data/MOL_00104.mol";
      bool takeOwnership = true;
      SDMolSupplier mol_supplier(file3, takeOwnership);
      std::unique_ptr<ROMol> m1(mol_supplier.next());
      REQUIRE(m1);
      RWMol m2(*m1);
      RDKit::MolOps::Kekulize(m2);
      std::vector<int> expected{-3, -1, -2, 0, 2, -3, -1, -2, 0, -1, -1, 2};
      Descriptors::calcOxidationNumbers(m2);
      for (auto &a : m2.atoms()) {
        CHECK(a->getProp<int>(common_properties::OxidationNumber) ==
              expected[a->getIdx()]);
      }
      RDKit::MolOps::hapticBondsToDative(m2);
      Descriptors::calcOxidationNumbers(m2);
      std::vector<int> expectedNoDummies{-3, -1, -2, 2, -3, -1, -2, -1, -1, 2};
      for (auto &a : m2.atoms()) {
        CHECK(a->getProp<int>(common_properties::OxidationNumber) ==
              expectedNoDummies[a->getIdx()]);
      }
    }
  }
  SECTION("Syngenta tests") {
    // These are from
    // https://github.com/syngenta/linchemin/blob/main/tests/cheminfo/test_functions.py#L385
    // and thus subject to the MIT license at
    //
    // https://github.com/syngenta/linchemin/blob/f44fda38e856eaa876483c94284ee6788d2c27f4/LICENSE
    std::vector<std::tuple<std::string, std::string, std::map<int, int>>>
        test_set{
            {"A-001",
             "O=S(=O)(O)O",
             {{0, -2}, {1, 6}, {2, -2}, {3, -2}, {4, -2}}},
            {"A-002",
             "NS(=O)(=O)O",
             {{0, -3}, {1, 6}, {2, -2}, {3, -2}, {4, -2}}},
            {"A-003",
             "NS(=O)(=O)c1ccccc1",
             {{0, -3},
              {1, 4},
              {2, -2},
              {3, -2},
              {4, 1},
              {5, -1},
              {6, -1},
              {7, -1},
              {8, -1},
              {9, -1}}},
            {"A-004",
             "O=S(=O)(O)c1ccccc1",
             {{0, -2},
              {1, 4},
              {2, -2},
              {3, -2},
              {4, 1},
              {5, -1},
              {6, -1},
              {7, -1},
              {8, -1},
              {9, -1}}},
            {"A-005",
             "O=S(=O)(Cl)c1ccccc1",
             {{0, -2},
              {1, 4},
              {2, -2},
              {3, -1},
              {4, 1},
              {5, -1},
              {6, -1},
              {7, -1},
              {8, -1},
              {9, -1}}},
            {"A-006", "S", {{0, -2}}},
            {"A-007", "CSC", {{0, -2}, {1, -2}, {2, -2}}},
            {"A-008", "CS(C)=O", {{0, -2}, {1, 0}, {2, -2}, {3, -2}}},
            {"A-009",
             "CS(C)(=O)=O",
             {{0, -2}, {1, 2}, {2, -2}, {3, -2}, {4, -2}}},
            {"A-010",
             "COS(=O)(=O)OC",
             {{0, -2}, {1, -2}, {2, 6}, {3, -2}, {4, -2}, {5, -2}, {6, -2}}},
            {"A-011",
             "COP(=O)(OC)OC",
             {{0, -2},
              {1, -2},
              {2, 5},
              {3, -2},
              {4, -2},
              {5, -2},
              {6, -2},
              {7, -2}}},
            {"A-012",
             "COP(OC)OC",
             {{0, -2}, {1, -2}, {2, 3}, {3, -2}, {4, -2}, {5, -2}, {6, -2}}},
            {"A-013",
             "COP(=O)(C#N)OC",
             {{0, -2},
              {1, -2},
              {2, 5},
              {3, -2},
              {4, 2},
              {5, -3},
              {6, -2},
              {7, -2}}},
            {"A-014",
             "CCP(=O)(CC)CC",
             {{0, -3},
              {1, -3},
              {2, 5},
              {3, -2},
              {4, -3},
              {5, -3},
              {6, -3},
              {7, -3}}},
            {"A-015",
             "CCP(CC)CC",
             {{0, -3}, {1, -3}, {2, 3}, {3, -3}, {4, -3}, {5, -3}, {6, -3}}},
            {"A-016",
             "CC[P+](CC)(CC)CC",
             {{0, -3},
              {1, -3},
              {2, 5},
              {3, -3},
              {4, -3},
              {5, -3},
              {6, -3},
              {7, -3},
              {8, -3}}},
            {"A-017",
             "c1ccncc1",
             {{0, -1}, {1, -1}, {2, 1}, {3, -3}, {4, 0}, {5, -1}}},
            {"A-018",
             "C[n+]1ccccc1",
             {{0, -2}, {1, -3}, {2, 1}, {3, -1}, {4, -1}, {5, -1}, {6, 0}}},
            {"A-019",
             "[O-][n+]1ccccc1",
             {{0, -2}, {1, -1}, {2, 1}, {3, -1}, {4, -1}, {5, -1}, {6, 0}}},
            {"A-020",
             "C[C-](C)[n+]1ccccc1",
             {{0, -3},
              {1, 0},
              {2, -3},
              {3, -3},
              {4, 1},
              {5, -1},
              {6, -1},
              {7, -1},
              {8, 0}}},
            {"A-021",
             "C1CCNCC1",
             {{0, -2}, {1, -2}, {2, -1}, {3, -3}, {4, -1}, {5, -2}}},
            {"A-022",
             "[O]N1CCCCC1",
             {{0, -1}, {1, -1}, {2, -1}, {3, -2}, {4, -2}, {5, -2}, {6, -1}}},
            {"A-023", "N", {{0, -3}}},
            {"A-024", "CN(C)C", {{0, -2}, {1, -3}, {2, -2}, {3, -2}}},
            {"A-025", "NO", {{0, -1}, {1, -2}}},
            {"A-026", "[NH4+]", {{0, -3}}},
            {"A-027",
             "C[N+](C)(C)C",
             {{0, -2}, {1, -3}, {2, -2}, {3, -2}, {4, -2}}},
            {"A-028",
             "C[N+](C)(C)[O-]",
             {{0, -2}, {1, -1}, {2, -2}, {3, -2}, {4, -2}}},
            {"A-029", "[SiH4]", {{0, 4}}},
            {"A-030",
             "C[Si](C)(C)C",
             {{0, -4}, {1, 4}, {2, -4}, {3, -4}, {4, -4}}},
            {"A-031",
             "C[Si](C)(C)Cl",
             {{0, -4}, {1, 4}, {2, -4}, {3, -4}, {4, -1}}},
            {"A-032",
             "C[Si](C)(C)O",
             {{0, -4}, {1, 4}, {2, -4}, {3, -4}, {4, -2}}},
            {"A-033", "C", {{0, -4}}},
            {"A-034", "CO", {{0, -2}, {1, -2}}},
            {"A-035", "C=O", {{0, 0}, {1, -2}}},
            {"A-036", "O=CO", {{0, -2}, {1, 2}, {2, -2}}},
            {"A-037", "O=C(O)O", {{0, -2}, {1, 4}, {2, -2}, {3, -2}}},
            {"A-038", "O=C=O", {{0, -2}, {1, 4}, {2, -2}}},
            {"A-039", "[C-]#[O+]", {{0, 2}, {1, -2}}},
            {"A-041", "CI", {{0, -2}, {1, -1}}},
            {"A-042", "ICI", {{0, -1}, {1, 0}, {2, -1}}},
            {"A-043", "IC(I)I", {{0, -1}, {1, 2}, {2, -1}, {3, -1}}},
            {"A-044",
             "IC(I)(I)I",
             {{0, -1}, {1, 4}, {2, -1}, {3, -1}, {4, -1}}},
            {"A-045",
             "FC(F)(F)I",
             {{0, -1}, {1, 4}, {2, -1}, {3, -1}, {4, -1}}},
            {"A-046", "II", {{0, 0}, {1, 0}}},
            {"A-047", "ClI", {{0, -1}, {1, 1}}},
            {"A-048",
             "[O-][I+3]([O-])([O-])[O-]",
             {{0, -2}, {1, 7}, {2, -2}, {3, -2}, {4, -2}}},
            {"A-049",
             "[O-][I+2]([O-])[O-]",
             {{0, -2}, {1, 5}, {2, -2}, {3, -2}}},
            {"A-050",
             "O=[I+]([O-])c1ccccc1",
             {{0, -2},
              {1, 3},
              {2, -2},
              {3, 1},
              {4, -1},
              {5, -1},
              {6, -1},
              {7, -1},
              {8, -1}}},
            {"A-051",
             "Ic1ccccc1",
             {{0, -1}, {1, 1}, {2, -1}, {3, -1}, {4, -1}, {5, -1}, {6, -1}}},
            {"A-052",
             "CC(=O)OI1(OC(C)=O)(OC(C)=O)OC(=O)c2ccccc21",
             {{0, -3},  {1, 3},   {2, -2},  {3, -2},  {4, 3},  {5, -2},
              {6, 3},   {7, -3},  {8, -2},  {9, -2},  {10, 3}, {11, -3},
              {12, -2}, {13, -2}, {14, 3},  {15, -2}, {16, 0}, {17, -1},
              {18, -1}, {19, -1}, {20, -1}, {21, 1}}},
            {"A-053", "[Cl-]", {{0, -1}}},
            {"A-054", "ClCl", {{0, 0}, {1, 0}}},
            {"A-055", "[O-]Cl", {{0, -2}, {1, 1}}},
            {"A-056", "[O-][Cl+][O-]", {{0, -2}, {1, 3}, {2, -2}}},
            {"A-057",
             "[O-][Cl+2]([O-])[O-]",
             {{0, -2}, {1, 5}, {2, -2}, {3, -2}}},
            {"A-58",
             "[O-][Cl+3]([O-])([O-])[O-]",
             {{0, -2}, {1, 7}, {2, -2}, {3, -2}, {4, -2}}}};

    for (const auto &test : test_set) {
      //      std::cout << "checking " << std::get<0>(test) << " : "
      //                << std::get<1>(test) << std::endl;
      std::unique_ptr<RWMol> mol(RDKit::SmilesToMol(std::get<1>(test)));
      RDKit::MolOps::Kekulize(*mol);
      Descriptors::calcOxidationNumbers(*mol);
      for (const auto &expected : std::get<2>(test)) {
        auto atom = mol->getAtomWithIdx(expected.first);
        CHECK(atom->getProp<int>(common_properties::OxidationNumber) ==
              expected.second);
      }
    }
  }
}

TEST_CASE("DCLV") {
  const PeriodicTable *tbl = PeriodicTable::getTable();
  std::string pathName = getenv("RDBASE");
  std::string pdbName =
      pathName + "/Code/GraphMol/Descriptors/test_data/1mup.pdb";
  auto m = v2::FileParsers::MolFromPDBFile(pdbName);
  std::vector<double> radii;
  for (const auto atom : m->atoms()) {
    radii.push_back(tbl->getRvdw(atom->getAtomicNum()));
  }
  REQUIRE(m);
  SECTION("defaults") {
    const bool isProtein = true;
    const bool includeLigand = false;
    const double probeRadius = 1.4;
    const int confId = -1;
    for (auto dclv :
         {Descriptors::DoubleCubicLatticeVolume(
              *m, radii, isProtein, includeLigand, probeRadius, confId),
          Descriptors::DoubleCubicLatticeVolume(*m, isProtein, includeLigand,
                                                probeRadius, confId)}) {
      CHECK(dclv.getSurfaceArea() == Catch::Approx(8306.62).epsilon(0.05));
      CHECK(dclv.getPolarSurfaceArea(false, false) ==
            Catch::Approx(4652.19).epsilon(0.05));  // N & O, no Hs
      CHECK(dclv.getPolarSurfaceArea(true, false) ==
            Catch::Approx(4673.9).epsilon(0.05));  // N, O, S & P, no Hs
      CHECK(dclv.getVolume() == Catch::Approx(29952.3).epsilon(0.05));
      CHECK(dclv.getVDWVolume() == Catch::Approx(13541.5).epsilon(0.05));
      CHECK(dclv.getPolarVolume(false, false) ==
            Catch::Approx(17096).epsilon(0.05));  // N & O, no Hs
      CHECK(dclv.getCompactness() == Catch::Approx(1.78096).epsilon(0.05));
      CHECK(dclv.getPackingDensity() == Catch::Approx(0.452103).epsilon(0.05));
    }
  }
  SECTION("non-default radius") {
    const bool isProtein = true;
    const bool includeLigand = false;
    const double probeRadius = 1.6;
    const int confId = -1;
    Descriptors::DoubleCubicLatticeVolume dclv(
        *m, radii, isProtein, includeLigand, probeRadius, confId);
    CHECK(dclv.getSurfaceArea() == Catch::Approx(8112.42).epsilon(0.05));
    CHECK(dclv.getPolarSurfaceArea(false, false) ==
          Catch::Approx(4656.8).epsilon(0.05));  // N & O, no Hs
    CHECK(dclv.getPolarSurfaceArea(true, false) ==
          Catch::Approx(4674.96).epsilon(0.05));  // N, O, S & P, no Hs
    CHECK(dclv.getVolume() == Catch::Approx(31591).epsilon(0.05));
    CHECK(dclv.getVDWVolume() == Catch::Approx(13541.5).epsilon(0.05));
    CHECK(dclv.getPolarVolume(false, false) ==
          Catch::Approx(18402.8).epsilon(0.05));  // N & O, no Hs
    CHECK(dclv.getCompactness() == Catch::Approx(1.67864).epsilon(0.05));
    CHECK(dclv.getPackingDensity() == Catch::Approx(0.428652).epsilon(0.05));
  }
  SECTION("include ligand") {
    const bool isProtein = true;
    const bool includeLigand = true;
    const double probeRadius = 1.4;
    const int confId = -1;
    Descriptors::DoubleCubicLatticeVolume dclv(
        *m, radii, isProtein, includeLigand, probeRadius, confId);
    CHECK(dclv.getSurfaceArea() == Catch::Approx(8206.1).epsilon(0.05));
    CHECK(dclv.getPolarSurfaceArea(false, false) ==
          Catch::Approx(4467.92).epsilon(0.05));  // N & O, no Hs
    CHECK(dclv.getPolarSurfaceArea(true, false) ==
          Catch::Approx(4481.19).epsilon(0.05));  // N, O, S & P, no Hs
    CHECK(dclv.getVolume() == Catch::Approx(30340.1).epsilon(0.05));
    CHECK(dclv.getVDWVolume() == Catch::Approx(13789.7).epsilon(0.05));
    CHECK(dclv.getPolarVolume(false, false) ==
          Catch::Approx(16507).epsilon(0.05));  // N & O, no Hs
    CHECK(dclv.getCompactness() == Catch::Approx(1.74438).epsilon(0.05));
    CHECK(dclv.getPackingDensity() == Catch::Approx(0.454504).epsilon(0.05));
  }
  SECTION("SDF") {
    std::string sdfName =
        pathName + "/Code/GraphMol/Descriptors/test_data/TZL_model.sdf";
    auto m = v2::FileParsers::MolFromMolFile(sdfName);
    std::vector<double> radii;
    for (const auto atom : m->atoms()) {
      radii.push_back(tbl->getRvdw(atom->getAtomicNum()));
    }
    REQUIRE(m);
    Descriptors::DoubleCubicLatticeVolume dclv(*m, radii);

    // packing density and compactness not relevant for small mol
    CHECK(dclv.getSurfaceArea() == Catch::Approx(304.239).epsilon(0.05));
    CHECK(dclv.getPolarSurfaceArea(false, false) ==
          Catch::Approx(18.7319).epsilon(0.05));  // N & O, no Hs
    CHECK(dclv.getPolarSurfaceArea(true, false) ==
          Catch::Approx(64.9764).epsilon(0.05));  // N, O, S & P, no Hs
    CHECK(dclv.getVolume() == Catch::Approx(431.35).epsilon(0.05));
    CHECK(dclv.getVDWVolume() == Catch::Approx(119.296).epsilon(0.05));
    CHECK(dclv.getPolarVolume(false, false) ==
          Catch::Approx(21.35).epsilon(0.05));  // N & O, no Hs
  }
  SECTION("Partial Test") {
    std::string sdfName =
        pathName + "/Code/GraphMol/Descriptors/test_data/ethane.sdf";
    auto m = v2::FileParsers::MolFromMolFile(sdfName);
    REQUIRE(m);

    std::vector<double> radii = {1.7, 1.7};

    const bool isProtein = false;
    const bool includeLigand = true;  // needs to be true for small mol
    const double probeRadius = 1.4;
    const int confId = -1;
    Descriptors::DoubleCubicLatticeVolume dclv(
        *m, radii, isProtein, includeLigand, probeRadius, confId);

    boost::dynamic_bitset<> partialAtoms(m->getNumAtoms());
    partialAtoms.set(0);

    boost::dynamic_bitset<> partialAtoms2(m->getNumAtoms());
    partialAtoms2.set(1);

    boost::dynamic_bitset<> partialAtoms3(m->getNumAtoms());
    partialAtoms3.set(0);
    partialAtoms3.set(1);

    CHECK(dclv.getAtomSurfaceArea(0) == dclv.getAtomSurfaceArea(1));
    CHECK(dclv.getPartialSurfaceArea(partialAtoms) ==
          dclv.getPartialSurfaceArea(partialAtoms2));
    CHECK(dclv.getPartialSurfaceArea(partialAtoms3) == dclv.getSurfaceArea());
  }

  SECTION("Surface points") {
    std::string sdfName =
        pathName + "/Code/GraphMol/Descriptors/test_data/ethane.sdf";
    auto m = v2::FileParsers::MolFromMolFile(sdfName);
    REQUIRE(m);
    Descriptors::DoubleCubicLatticeVolume dclv(*m);
    auto surfacePoints = dclv.getSurfacePoints();
    CHECK(surfacePoints.size() == m->getNumAtoms());
    auto allSurfacePoints = dclv.getSurfacePoints(true);
    CHECK(allSurfacePoints.size() == surfacePoints.size());
    for (size_t i = 0; i < surfacePoints.size(); ++i) {
      CHECK(allSurfacePoints[i].size() == 320);
      CHECK(surfacePoints[i].size() < allSurfacePoints[i].size());
    }
    // make sure calling the function again doesn't change the result
    auto surfacePoints2 = dclv.getSurfacePoints();
    CHECK(surfacePoints2.size() == m->getNumAtoms());
    for (size_t i = 0; i < surfacePoints2.size(); ++i) {
      CHECK(surfacePoints2[i].size() == surfacePoints[i].size());
    }
  }
}
#ifdef RDK_HAS_EIGEN3
TEST_CASE("Github #7364: BCUT descriptors failing for moleucles with Hs") {
  SECTION("as reported") {
    auto m = "CCOC#CCC(C(=O)c1ccc(C)cc1)N1CCCC1"_smiles;
    REQUIRE(m);
    RWMol m2 = *m;
    MolOps::addHs(m2);
    auto ref = Descriptors::BCUT2D(*m);
    auto val = Descriptors::BCUT2D(m2);
    CHECK(ref.size() == val.size());
    CHECK(ref == val);
  }
}
#endif

TEST_CASE(
    "Github #6757: numAtomStereoCenters fails if molecule is sanitized a second time") {
  SECTION("as reported") {
    auto m = "C[C@H](F)Cl"_smiles;
    REQUIRE(m);
    CHECK(Descriptors::numAtomStereoCenters(*m) == 1);
    CHECK(Descriptors::numUnspecifiedAtomStereoCenters(*m) == 0);
    MolOps::sanitizeMol(*m);
    CHECK(Descriptors::numAtomStereoCenters(*m) == 1);
  }
  SECTION("expanded") {
    auto m = "C[C@H](F)C(O)Cl"_smiles;
    REQUIRE(m);
    CHECK(Descriptors::numAtomStereoCenters(*m) == 2);
    CHECK(Descriptors::numUnspecifiedAtomStereoCenters(*m) == 1);
  }
}

#ifdef RDK_BUILD_DESCRIPTORS3D
TEST_CASE("Github #7264: GETAWAY descriptors are non-deterministic") {
  SECTION("as reported") {
    v2::SmilesParse::SmilesParserParams ps;
    ps.removeHs = false;
    auto m = v2::SmilesParse::MolFromSmiles(
        "[H]c1snnc1Br |(0.753469,1.80182,0.0404378;0.108659,0.921186,-0.0974666;-1.6304,0.84393,-0.147597;-1.60936,-1.01967,-0.408092;-0.507386,-1.21741,-0.410872;0.515991,-0.400518,-0.273344;2.36902,-0.929334,-0.305519)|",
        ps);
    REQUIRE(m);
    /*
    void GETAWAY(
    const ROMol &, std::vector<double> &res, int confId = -1,
    unsigned int precision = 2, const std::string &customAtomPropName = "");
    */
    std::vector<double> res1, res2;
    Descriptors::GETAWAY(*m, res1);
    Descriptors::GETAWAY(*m, res2);
    for (size_t i = 0; i < res1.size(); ++i) {
      INFO(i);
      CHECK(res1[i] == res2[i]);
    }
  }
}

TEST_CASE("GETAWAY on molecules with several fragments") {
  // Atoms in different fragments are a very large topological distance apart;
  // they only contribute to the totals, which therefore match those of the
  // connected molecule with the same atoms and coordinates.
  // Covers the default path (HTu, HATSu, RTu) and the custom atom property
  // path (HTc, HATSTc, RTc).
  struct Case {
    std::string propName;
    size_t numValues;
    std::vector<size_t> totals;
  };
  const auto tc = GENERATE(Case{"", 273, {13, 23, 155}},
                           Case{"customProp", 45, {13, 23, 35}});
  INFO("customAtomPropName: '" << tc.propName << "'");
  const std::string coords =
      " |(-1.2,0.1,0.3;0.2,-0.4,0.1;1.1,0.6,-0.3;2.4,0.2,0.4)|";
  auto salt = v2::SmilesParse::MolFromSmiles("CCO.C" + coords);
  auto connected = v2::SmilesParse::MolFromSmiles("CCOC" + coords);
  REQUIRE(salt);
  REQUIRE(connected);
  for (auto mol : {salt.get(), connected.get()}) {
    for (auto atom : mol->atoms()) {
      atom->setProp("customProp", static_cast<double>(atom->getAtomicNum()));
    }
  }
  std::vector<double> saltRes, connectedRes;
  Descriptors::GETAWAY(*salt, saltRes, -1, 2, tc.propName);
  Descriptors::GETAWAY(*connected, connectedRes, -1, 2, tc.propName);
  REQUIRE(saltRes.size() == tc.numValues);
  REQUIRE(connectedRes.size() == tc.numValues);
  for (const auto idx : tc.totals) {
    INFO(idx);
    CHECK(saltRes[idx] == connectedRes[idx]);
  }
}
#endif

#ifdef RDK_BUILD_DESCRIPTORS3D
TEST_CASE("EEM charges of charged molecules sum to the formal charge") {
  v2::SmilesParse::SmilesParserParams ps;
  ps.removeHs = false;
  SECTION("hydroxide") {
    auto m = v2::SmilesParse::MolFromSmiles("[O-][H] |(0,0,0;0.97,0,0)|", ps);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 2);
    std::vector<double> charges;
    Descriptors::EEM(*m, charges, -1);
    REQUIRE(charges.size() == 2);
    CHECK(charges[0] == Catch::Approx(-1.0632).margin(1e-4));
    CHECK(charges[1] == Catch::Approx(0.0632).margin(1e-4));
  }
  SECTION("sum of charges") {
    const std::vector<std::string> smis = {
        "[O-][H] |(0,0,0;0.97,0,0)|",
        "[H][O+]([H])[H] |(0.94,0,0;0,0,0;-0.31,0.89,0;-0.31,-0.44,0.77)|",
        "[H]C([H])([H])C(=O)[O-] |(-1.06,-0.98,0.31;-0.49,-0.06,0.28;-0.69,0.48,1.21;-0.83,0.55,-0.56;0.99,-0.30,0.17;1.70,0.71,0.09;1.40,-1.48,0.17)|"};
    for (const auto &smi : smis) {
      INFO(smi);
      auto m = v2::SmilesParse::MolFromSmiles(smi, ps);
      REQUIRE(m);
      std::vector<double> charges;
      Descriptors::EEM(*m, charges, -1);
      double sum = 0.0;
      for (auto q : charges) {
        sum += q;
      }
      CHECK(sum == Catch::Approx(MolOps::getFormalCharge(*m)).margin(1e-8));
    }
  }
}
#endif
TEST_CASE("DCLV ignores dummy atoms") {
  v2::SmilesParse::SmilesParserParams ps;
  ps.removeHs = false;
  auto checkSame = [&ps](const std::string &smi, const std::string &refSmi) {
    INFO(smi << " vs " << refSmi);
    auto m = v2::SmilesParse::MolFromSmiles(smi, ps);
    REQUIRE(m);
    auto ref = v2::SmilesParse::MolFromSmiles(refSmi, ps);
    REQUIRE(ref);
    Descriptors::DoubleCubicLatticeVolume dclv(*m);
    Descriptors::DoubleCubicLatticeVolume refDclv(*ref);
    CHECK(dclv.getSurfaceArea() == refDclv.getSurfaceArea());
    CHECK(dclv.getVolume() == refDclv.getVolume());
    CHECK(dclv.getVDWVolume() == refDclv.getVDWVolume());
    CHECK(dclv.getCompactness() == refDclv.getCompactness());
    CHECK(dclv.getPackingDensity() == refDclv.getPackingDensity());
  };
  SECTION("leading dummy") { checkSame("*C |(1.5,0,0;0,0,0)|", "C |(0,0,0)|"); }
  SECTION("trailing dummy") {
    checkSame("C* |(0,0,0;1.5,0,0)|", "C |(0,0,0)|");
    checkSame("CCO* |(0,0,0;1.5,0,0;2,1.3,0;8,8,8)|",
              "CCO |(0,0,0;1.5,0,0;2,1.3,0)|");
  }
}

TEST_CASE("DCLV polar contributions with H-H bonds") {
  v2::SmilesParse::SmilesParserParams ps;
  ps.removeHs = false;
  SECTION("hydrogen molecule") {
    auto m = v2::SmilesParse::MolFromSmiles("[H][H] |(0,0,0;0.74,0,0)|", ps);
    REQUIRE(m);
    Descriptors::DoubleCubicLatticeVolume dclv(*m);
    CHECK(dclv.getPolarSurfaceArea(false, true) == 0.0);
    CHECK(dclv.getPolarVolume(false, true) == 0.0);
    CHECK(dclv.getPolarSurfaceArea(true, true) == 0.0);
    CHECK(dclv.getPolarVolume(true, true) == 0.0);
  }
  SECTION("polar hydrogens are still included") {
    // methanol: the polar atoms are O (4) and the H attached to it (5)
    auto m = v2::SmilesParse::MolFromSmiles(
        "[H]C([H])([H])O[H] |(-0.36,1.03,0;0,0,0;-0.36,-0.51,0.89;-0.36,-0.51,-0.89;1.43,0,0;1.75,-0.91,0)|",
        ps);
    REQUIRE(m);
    Descriptors::DoubleCubicLatticeVolume dclv(*m);
    boost::dynamic_bitset<> polarAtoms(m->getNumAtoms());
    polarAtoms.set(4);
    polarAtoms.set(5);
    CHECK(dclv.getPolarSurfaceArea(false, true) > 0.0);
    CHECK(dclv.getPolarSurfaceArea(false, true) ==
          dclv.getPartialSurfaceArea(polarAtoms));
    CHECK(dclv.getPolarVolume(false, true) ==
          dclv.getPartialVolume(polarAtoms));
  }
}

namespace {
// lets tests use arbitrary functors without touching the global registry
class PropertiesFromFunctors : public Descriptors::Properties {
 public:
  PropertiesFromFunctors(
      std::vector<boost::shared_ptr<Descriptors::PropertyFunctor>> props) {
    m_properties = std::move(props);
  }
};

struct FailingProperty : public Descriptors::PropertyFunctor {
  FailingProperty() : PropertyFunctor("FailOnSingleAtom", "1.0.0") {}
  double operator()(const ROMol &mol) const override {
    if (mol.getNumAtoms() == 1) {
      throw ValueErrorException("single atom");
    }
    return mol.getNumAtoms();
  }
};

struct BadAllocProperty : public Descriptors::PropertyFunctor {
  BadAllocProperty() : PropertyFunctor("BadAlloc", "1.0.0") {}
  // bad_alloc is rethrown whatever the failure mode says
  Descriptors::PropertyFailureMode failureMode() const override {
    return Descriptors::PropertyFailureMode::DEFAULT_ZERO;
  }
  double operator()(const ROMol &) const override { throw std::bad_alloc(); }
};

// always throws, with THROW failure mode, so the exception must reach the
// caller
struct ThrowModeProperty : public Descriptors::PropertyFunctor {
  ThrowModeProperty() : PropertyFunctor("ThrowMode", "1.0.0") {}
  Descriptors::PropertyFailureMode failureMode() const override {
    return Descriptors::PropertyFailureMode::THROW;
  }
  double operator()(const ROMol &) const override {
    throw ValueErrorException("throw mode failure");
  }
};

// like FailingProperty, but failures give 0.0
struct ZeroModeProperty : public Descriptors::PropertyFunctor {
  ZeroModeProperty() : PropertyFunctor("ZeroMode", "1.0.0") {}
  Descriptors::PropertyFailureMode failureMode() const override {
    return Descriptors::PropertyFailureMode::DEFAULT_ZERO;
  }
  double operator()(const ROMol &mol) const override {
    if (mol.getNumAtoms() == 1) {
      throw ValueErrorException("single atom");
    }
    return mol.getNumAtoms();
  }
};

// records the id of the thread each call runs on
struct ThreadRecordingProperty : public Descriptors::PropertyFunctor {
  std::shared_ptr<std::mutex> d_mutex;
  std::shared_ptr<std::vector<std::thread::id>> d_ids;
  ThreadRecordingProperty(const std::string &name,
                          std::shared_ptr<std::mutex> mutex,
                          std::shared_ptr<std::vector<std::thread::id>> ids)
      : PropertyFunctor(name, "1.0.0"),
        d_mutex(std::move(mutex)),
        d_ids(std::move(ids)) {}
  double operator()(const ROMol &mol) const override {
    {
      std::lock_guard<std::mutex> lock(*d_mutex);
      d_ids->push_back(std::this_thread::get_id());
    }
    return mol.getNumAtoms();
  }
};

// saves the static property registry and restores it on scope exit, so tests
// that register properties cannot leak into other tests even on failure
class RegistryRestore {
  std::vector<boost::shared_ptr<Descriptors::PropertyFunctor>> d_saved;

 public:
  RegistryRestore() : d_saved(Descriptors::Properties::registry) {}
  ~RegistryRestore() { Descriptors::Properties::registry = d_saved; }
};

bool sameRow(const std::vector<double> &a, const std::vector<double> &b) {
  if (a.size() != b.size()) {
    return false;
  }
  for (size_t i = 0; i < a.size(); ++i) {
    if (std::isnan(a[i]) != std::isnan(b[i]) ||
        (!std::isnan(a[i]) && a[i] != b[i])) {
      return false;
    }
  }
  return true;
}

const std::vector<std::string> batchSmiles{
    "c1ccccc1C(=O)O",
    "CCO",
    "C1CC1CC",
    "OC(=O)CCN",
    "c1ccncc1Cl",
    "CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O",
    "C",
    "O=S(=O)(N)c1ccccc1"};

std::vector<std::unique_ptr<RWMol>> parseBatch() {
  std::vector<std::unique_ptr<RWMol>> res;
  for (const auto &smi : batchSmiles) {
    res.emplace_back(SmilesToMol(smi));
    REQUIRE(res.back());
  }
  return res;
}

std::vector<const ROMol *> asPtrs(
    const std::vector<std::unique_ptr<RWMol>> &ms) {
  std::vector<const ROMol *> res;
  for (const auto &m : ms) {
    res.push_back(m.get());
  }
  return res;
}
}  // namespace

TEST_CASE("Properties::computePropertiesForMols", "[2D][Properties]") {
  Descriptors::Properties props;  // all registered properties
  // the reference values are computed on molecules which are not reused, so
  // cached results cannot leak into the batch calculations
  std::vector<std::vector<double>> ref;
  {
    auto refMols = parseBatch();
    for (const auto &m : refMols) {
      ref.push_back(props.computeProperties(*m));
    }
  }

  SECTION("matches computeProperties") {
    for (int numThreads : {1, 2, 4, 0, -1, 100}) {
      INFO("numThreads: " << numThreads);
      auto ms = parseBatch();
      auto res = props.computePropertiesForMols(asPtrs(ms), numThreads);
      REQUIRE(res.size() == ref.size());
      for (size_t i = 0; i < res.size(); ++i) {
        CHECK(sameRow(res[i], ref[i]));
      }
    }
  }
  SECTION("column order follows the property names") {
    Descriptors::Properties p1(std::vector<std::string>{"tpsa", "exactmw"});
    Descriptors::Properties p2(std::vector<std::string>{"exactmw", "tpsa"});
    auto ms = parseBatch();
    auto r1 = p1.computePropertiesForMols(asPtrs(ms), 2);
    auto r2 = p2.computePropertiesForMols(asPtrs(ms), 2);
    for (size_t i = 0; i < ms.size(); ++i) {
      CHECK(r1[i][0] == r2[i][1]);
      CHECK(r1[i][1] == r2[i][0]);
    }
  }
  SECTION("null molecules") {
    for (int numThreads : {1, 4}) {
      auto ms = parseBatch();
      auto ptrs = asPtrs(ms);
      ptrs.insert(ptrs.begin(), nullptr);
      ptrs.insert(ptrs.begin() + 4, nullptr);
      ptrs.push_back(nullptr);
      auto res = props.computePropertiesForMols(ptrs, numThreads);
      REQUIRE(res.size() == ptrs.size());
      size_t refIdx = 0;
      for (size_t i = 0; i < ptrs.size(); ++i) {
        REQUIRE(res[i].size() == ref[0].size());
        if (!ptrs[i]) {
          for (auto v : res[i]) {
            CHECK(std::isnan(v));
          }
        } else {
          CHECK(sameRow(res[i], ref[refIdx++]));
        }
      }
    }
  }
  SECTION("empty input and no properties") {
    for (int numThreads : {1, 4, 0}) {
      CHECK(props.computePropertiesForMols({}, numThreads).empty());
    }
    Descriptors::Properties none(std::vector<std::string>{});
    auto ms = parseBatch();
    auto res = none.computePropertiesForMols(asPtrs(ms), 4);
    REQUIRE(res.size() == ms.size());
    for (const auto &row : res) {
      CHECK(row.empty());
    }
  }
  SECTION("repeated molecules") {
    // 3 threads: with a stride that doesn't divide the number of distinct
    // molecules, repeats would be handed to different threads
    for (int numThreads : {1, 3, 4}) {
      auto ms = parseBatch();
      std::vector<const ROMol *> ptrs;
      for (unsigned int rep = 0; rep < 16; ++rep) {
        for (const auto &m : ms) {
          ptrs.push_back(m.get());
        }
      }
      auto res = props.computePropertiesForMols(ptrs, numThreads);
      REQUIRE(res.size() == ptrs.size());
      for (size_t i = 0; i < res.size(); ++i) {
        CHECK(sameRow(res[i], ref[i % ms.size()]));
      }
    }
  }
  SECTION("a single molecule repeated") {
    for (int numThreads : {1, 4}) {
      auto ms = parseBatch();
      std::vector<const ROMol *> ptrs(64, ms[5].get());
      auto res = props.computePropertiesForMols(ptrs, numThreads);
      REQUIRE(res.size() == ptrs.size());
      for (const auto &row : res) {
        CHECK(sameRow(row, ref[5]));
      }
    }
  }
  SECTION("calculation failures produce NaN") {
    PropertiesFromFunctors failing(
        {boost::shared_ptr<Descriptors::PropertyFunctor>(new FailingProperty),
         Descriptors::Properties::getProperty("NumHeavyAtoms")});
    auto m1 = "C"_smiles;
    auto m2 = "CC"_smiles;
    REQUIRE(m1);
    REQUIRE(m2);
    for (int numThreads : {1, 4}) {
      auto res =
          failing.computePropertiesForMols({m1.get(), m2.get()}, numThreads);
      REQUIRE(res.size() == 2);
      CHECK(std::isnan(res[0][0]));
      CHECK(res[0][1] == 1);
      CHECK(res[1][0] == 2);
      CHECK(res[1][1] == 2);
    }
  }
  SECTION("other exceptions are propagated") {
    PropertiesFromFunctors badAlloc(
        {boost::shared_ptr<Descriptors::PropertyFunctor>(
            new BadAllocProperty)});
    auto ms = parseBatch();
    for (int numThreads : {1, 4}) {
      CHECK_THROWS_AS(badAlloc.computePropertiesForMols(asPtrs(ms), numThreads),
                      std::bad_alloc);
    }
  }
}

TEST_CASE("computePropertiesForMols failure modes", "[2D][Properties]") {
  // constructing a Properties registers the built-in descriptors
  Descriptors::Properties warm;
  auto ms = parseBatch();
  auto ptrs = asPtrs(ms);

  SECTION("THROW-mode exceptions propagate") {
    PropertiesFromFunctors p(
        {boost::shared_ptr<Descriptors::PropertyFunctor>(new ThrowModeProperty),
         Descriptors::Properties::getProperty("NumHeavyAtoms")});
    for (int numThreads : {1, 4}) {
      INFO("numThreads: " << numThreads);
      REQUIRE_THROWS_AS(p.computePropertiesForMols(ptrs, numThreads),
                        ValueErrorException);
    }
  }

  SECTION("DEFAULT_ZERO gives 0, DEFAULT_NAN gives NaN") {
    PropertiesFromFunctors p(
        {boost::shared_ptr<Descriptors::PropertyFunctor>(new ZeroModeProperty),
         boost::shared_ptr<Descriptors::PropertyFunctor>(new FailingProperty),
         Descriptors::Properties::getProperty("NumHeavyAtoms")});
    auto m1 = "C"_smiles;
    auto m2 = "CC"_smiles;
    REQUIRE(m1);
    REQUIRE(m2);
    for (int numThreads : {1, 4}) {
      INFO("numThreads: " << numThreads);
      auto res = p.computePropertiesForMols({m1.get(), m2.get()}, numThreads);
      REQUIRE(res.size() == 2);
      CHECK(res[0][0] == 0.0);
      CHECK(std::isnan(res[0][1]));
      CHECK(res[0][2] == 1);
      CHECK(res[1][0] == 2);
      CHECK(res[1][1] == 2);
      CHECK(res[1][2] == 2);
    }
  }

  SECTION("bad_alloc ignores the failure mode") {
    PropertiesFromFunctors p({boost::shared_ptr<Descriptors::PropertyFunctor>(
        new BadAllocProperty)});
    for (int numThreads : {1, 4}) {
      INFO("numThreads: " << numThreads);
      CHECK_THROWS_AS(p.computePropertiesForMols(ptrs, numThreads),
                      std::bad_alloc);
    }
  }

  SECTION("getProperties returns the instance's functors") {
    std::vector<boost::shared_ptr<Descriptors::PropertyFunctor>> functors{
        boost::shared_ptr<Descriptors::PropertyFunctor>(new ThrowModeProperty),
        Descriptors::Properties::getProperty("NumHeavyAtoms"),
        boost::shared_ptr<Descriptors::PropertyFunctor>(new ZeroModeProperty)};
    PropertiesFromFunctors p(functors);
    const auto &got = p.getProperties();
    REQUIRE(got.size() == functors.size());
    auto names = p.getPropertyNames();
    REQUIRE(names.size() == functors.size());
    for (size_t i = 0; i < got.size(); ++i) {
      CHECK(got[i] == functors[i]);
      CHECK(names[i] == got[i]->getName());
    }
  }

#ifdef RDK_TEST_MULTITHREADED
  SECTION("control: native functors do run on worker threads") {
    auto mutex = std::make_shared<std::mutex>();
    auto ids = std::make_shared<std::vector<std::thread::id>>();
    PropertiesFromFunctors p({boost::shared_ptr<Descriptors::PropertyFunctor>(
        new ThreadRecordingProperty("NativeRecord", mutex, ids))});
    p.computePropertiesForMols(ptrs, 4);
    REQUIRE(ids->size() == ptrs.size());
    bool sawOtherThread = false;
    for (const auto &id : *ids) {
      sawOtherThread = sawOtherThread || (id != std::this_thread::get_id());
    }
    CHECK(sawOtherThread);
  }
#endif

  SECTION("registry replacement does not change an existing instance") {
    RegistryRestore restore;
    boost::shared_ptr<Descriptors::PropertyFunctor> original(
        new ThrowModeProperty);
    Descriptors::Properties::registerProperty(original);
    Descriptors::Properties instance(
        std::vector<std::string>{original->getName()});
    // same name, different functor
    boost::shared_ptr<Descriptors::PropertyFunctor> replacement(
        new ThrowModeProperty);
    Descriptors::Properties::registerProperty(replacement);
    REQUIRE(instance.getProperties().size() == 1);
    CHECK(instance.getProperties()[0] == original);
    CHECK(Descriptors::Properties::getProperty(original->getName()) ==
          replacement);
  }
}

TEST_CASE("computePropertiesForMols TPSA reference values",
          "[2D][Properties]") {
  std::string fName = std::getenv("RDBASE");
  fName += "/Data/NCI/first_200.tpsa.csv";
  std::ifstream inf(fName);
  REQUIRE(inf.good());
  std::vector<std::unique_ptr<RWMol>> mols;
  std::vector<double> expected;
  std::string line;
  while (std::getline(inf, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    auto comma = line.find(',');
    REQUIRE(comma != std::string::npos);
    mols.emplace_back(SmilesToMol(line.substr(0, comma)));
    REQUIRE(mols.back());
    expected.push_back(std::stod(line.substr(comma + 1)));
  }
  REQUIRE(!mols.empty());
  Descriptors::Properties props(std::vector<std::string>{"tpsa"});
  auto res = props.computePropertiesForMols(asPtrs(mols), 4);
  REQUIRE(res.size() == expected.size());
  for (size_t i = 0; i < res.size(); ++i) {
    CHECK_THAT(res[i][0], Catch::Matchers::WithinAbs(expected[i], 1e-4));
  }
}

#ifdef RDK_TEST_MULTITHREADED
TEST_CASE("computePropertiesForMols parallel vs serial", "[2D][Properties]") {
  std::string fName = std::getenv("RDBASE");
  fName += "/Data/NCI/first_200.props.sdf";
  auto loadMols = [&fName]() {
    std::vector<std::unique_ptr<RWMol>> res;
    SDMolSupplier suppl(fName);
    while (!suppl.atEnd()) {
      std::unique_ptr<ROMol> m;
      try {
        m.reset(suppl.next());
      } catch (...) {
      }
      res.emplace_back(m ? new RWMol(*m) : nullptr);
    }
    return res;
  };
  Descriptors::Properties props;
  auto names = props.getPropertyNames();
  auto refMols = loadMols();
  auto ref = props.computePropertiesForMols(asPtrs(refMols), 1);

  // properties which are also stored as SD tags in the input file
  const std::vector<std::pair<std::string, std::string>> tags{
      {"NumHBD", "NUM_HDONORS"},
      {"NumHBA", "NUM_HACCEPTORS"},
      {"lipinskiHBD", "NUM_LIPINSKIHDONORS"},
      {"lipinskiHBA", "NUM_LIPINSKIHACCEPTORS"},
      {"NumRings", "NUM_RINGS"},
      {"NumHeteroatoms", "NUM_HETEROATOMS"}};
  for (const auto &[propName, tag] : tags) {
    auto col = std::find(names.begin(), names.end(), propName) - names.begin();
    REQUIRE(col < static_cast<long>(names.size()));
    for (size_t i = 0; i < refMols.size(); ++i) {
      if (refMols[i] && refMols[i]->hasProp(tag)) {
        CHECK(ref[i][col] == std::stoi(refMols[i]->getProp<std::string>(tag)));
      }
    }
  }

  for (int numThreads : {2, 4, 0}) {
    INFO("numThreads: " << numThreads);
    auto ms = loadMols();
    auto res = props.computePropertiesForMols(asPtrs(ms), numThreads);
    REQUIRE(res.size() == ref.size());
    for (size_t i = 0; i < res.size(); ++i) {
      CHECK(sameRow(res[i], ref[i]));
    }
    // a second run uses the values cached on the molecules
    auto again = props.computePropertiesForMols(asPtrs(ms), numThreads);
    for (size_t i = 0; i < res.size(); ++i) {
      CHECK(sameRow(again[i], res[i]));
    }
  }
}
#endif
