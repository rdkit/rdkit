//
//  Copyright (C) 2021 Greg Landrum
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
  std::string pathName = getenv("RDBASE");
  std::string pdbName =
      pathName + "/Code/GraphMol/Descriptors/test_data/1mup.pdb";
  auto m = v2::FileParsers::MolFromPDBFile(pdbName);
  REQUIRE(m);
  SECTION("defaults") {
    bool isProtein = true;
    bool includeLigand = false;
    Descriptors::DoubleCubicLatticeVolume dclv(*m, isProtein, includeLigand);
    CHECK(dclv.getSurfaceArea() == Catch::Approx(8330.59).epsilon(0.05));
    CHECK(dclv.getVolume() == Catch::Approx(31789.6).epsilon(0.05));
    CHECK(dclv.getVDWVolume() == Catch::Approx(15355.3).epsilon(0.05));
    CHECK(dclv.getCompactness() == Catch::Approx(1.7166).epsilon(0.05));
    CHECK(dclv.getPackingDensity() == Catch::Approx(0.48303).epsilon(0.05));
  }
  SECTION("depth and radius") {
    double probeRadius = 1.6;
    int depth = 6;
    bool isProtein = true;
    bool includeLigand = false;
    Descriptors::DoubleCubicLatticeVolume dclv(*m, isProtein, includeLigand,
                                               probeRadius, depth);
    CHECK(dclv.getSurfaceArea() == Catch::Approx(8186.06).epsilon(0.05));
    CHECK(dclv.getVolume() == Catch::Approx(33464.5).epsilon(0.05));
    CHECK(dclv.getVDWVolume() == Catch::Approx(15350.7).epsilon(0.05));
    CHECK(dclv.getCompactness() == Catch::Approx(1.63005).epsilon(0.05));
    CHECK(dclv.getPackingDensity() == Catch::Approx(0.458717).epsilon(0.05));
  }
  SECTION("ligand") {
    bool isProtein = true;
    bool includeLigand = true;
    Descriptors::DoubleCubicLatticeVolume dclv(*m, isProtein, includeLigand);
    CHECK(dclv.getSurfaceArea() == Catch::Approx(8010.56).epsilon(0.05));
    CHECK(dclv.getVolume() == Catch::Approx(31228.4).epsilon(0.05));
    CHECK(dclv.getVDWVolume() == Catch::Approx(15155.7).epsilon(0.05));
    CHECK(dclv.getCompactness() == Catch::Approx(1.67037).epsilon(0.05));
    CHECK(dclv.getPackingDensity() == Catch::Approx(0.48532).epsilon(0.05));
  }
  SECTION("SDF") {
    std::string sdfName =
        pathName + "/Code/GraphMol/Descriptors/test_data/TZL_model.sdf";
    auto m = v2::FileParsers::MolFromMolFile(sdfName);
    REQUIRE(m);
    bool isProtein = false;
    Descriptors::DoubleCubicLatticeVolume dclv(*m, isProtein);
    // NOTE - expected values generated from Roger's original C code
    // Original did not return surface area for Ligand only
    // so no check for Surface Area or Compactness

    CHECK(dclv.getSurfaceArea() == Catch::Approx(296.466).epsilon(0.05));
    CHECK(dclv.getVolume() == Catch::Approx(411.972).epsilon(0.05));
    CHECK(dclv.getVDWVolume() == Catch::Approx(139.97).epsilon(0.05));
  }
}

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