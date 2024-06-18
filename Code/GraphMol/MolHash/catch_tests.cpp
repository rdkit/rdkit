//
//  Copyright (C) 2019-2022 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include "MolHash.h"

#include <iostream>
#include <fstream>

using namespace RDKit;

TEST_CASE("Basic MolHash", "[molhash]") {
  SECTION("basics") {
    auto om = "C1CCCC(O)C1c1ccnc(OC)c1"_smiles;
    REQUIRE(om);
    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh =
          MolHash::MolHash(m.get(), MolHash::HashFunction::AnonymousGraph);
      CHECK(hsh == "***1****(*2*****2*)*1");
    }
    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh = MolHash::MolHash(m.get(), MolHash::HashFunction::ElementGraph);
      CHECK(hsh == "COC1CC(C2CCCCC2O)CCN1");
    }
    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh =
          MolHash::MolHash(m.get(), MolHash::HashFunction::CanonicalSmiles);
      CHECK(hsh == "COc1cc(C2CCCCC2O)ccn1");
    }
    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh =
          MolHash::MolHash(m.get(), MolHash::HashFunction::MurckoScaffold);
      CHECK(hsh == "c1cc(C2CCCCC2)ccn1");
    }
    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh =
          MolHash::MolHash(m.get(), MolHash::HashFunction::ExtendedMurcko);
      CHECK(hsh == "*c1cc(C2CCCCC2*)ccn1");
    }
    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh = MolHash::MolHash(m.get(), MolHash::HashFunction::MolFormula);
      CHECK(hsh == "C12H17NO2");
    }
    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh =
          MolHash::MolHash(m.get(), MolHash::HashFunction::AtomBondCounts);
      CHECK(hsh == "15,16");
    }
    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh = MolHash::MolHash(m.get(), MolHash::HashFunction::DegreeVector);
      CHECK(hsh == "0,4,9,2");
    }
    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh = MolHash::MolHash(m.get(), MolHash::HashFunction::Mesomer);
      CHECK(hsh == "CO[C]1[CH][C](C2CCCCC2O)[CH][CH][N]1_0");
    }
    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh = MolHash::MolHash(m.get(), MolHash::HashFunction::Regioisomer);
      CHECK(hsh == "*O.*O*.C.C1CCCCC1.c1ccncc1");
    }
    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh = MolHash::MolHash(m.get(), MolHash::HashFunction::NetCharge);
      CHECK(hsh == "0");
    }
    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh =
          MolHash::MolHash(m.get(), MolHash::HashFunction::SmallWorldIndexBR);
      CHECK(hsh == "B16R2");
    }
    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh =
          MolHash::MolHash(m.get(), MolHash::HashFunction::SmallWorldIndexBRL);
      CHECK(hsh == "B16R2L9");
    }
    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh = MolHash::MolHash(
          m.get(), MolHash::HashFunction::ArthorSubstructureOrder);
      CHECK(hsh == "000f001001000c000300005f000000");
    }
  }
  SECTION("tautomers") {
    auto om = "C(CC1=NNC=C1)C1=CNC=N1"_smiles;
    REQUIRE(om);
    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh =
          MolHash::MolHash(m.get(), MolHash::HashFunction::HetAtomTautomer);
      CHECK(hsh == "[CH]1[CH][C](CC[C]2[CH][N][CH][N]2)[N][N]1_2_0");
    }
    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh =
          MolHash::MolHash(m.get(), MolHash::HashFunction::HetAtomProtomer);
      CHECK(hsh == "[CH]1[CH][C](CC[C]2[CH][N][CH][N]2)[N][N]1_2");
    }
  }
  SECTION("tautomers 2") {
    {
      auto om = "C/C=C/C"_smiles;
      REQUIRE(om);
      auto hsh =
          MolHash::MolHash(om.get(), MolHash::HashFunction::HetAtomTautomer);
      CHECK(hsh == "C/C=C/C_0_0");
    }

    {
      auto om = "C/C=N/C"_smiles;
      REQUIRE(om);
      auto hsh =
          MolHash::MolHash(om.get(), MolHash::HashFunction::HetAtomTautomer);
      CHECK(hsh == "C[CH][N]C_0_0");
    }

    {
      auto om = "C/C=C/C=C/C"_smiles;
      REQUIRE(om);
      auto hsh =
          MolHash::MolHash(om.get(), MolHash::HashFunction::HetAtomTautomer);
      CHECK(hsh == "C[CH][CH][CH][CH]C_0_0");
    }
  }
  SECTION("tautomers bug found in testing") {
    auto m1 =
        "CCC(=Cc1sc2cc(C)c(C)cc2[n+]1CC(O)CS(=O)(=O)[O-])C=C1[Se]c2ccc(C)cc2[NH+]1CC"_smiles;
    REQUIRE(m1);
    auto m2 =
        "CCC(=Cc1[se]c2ccc(C)cc2[n+]1CC)C=C1Sc2cc(C)c(C)cc2N1CC(O)CS(=O)(=O)O"_smiles;
    REQUIRE(m2);

    std::unique_ptr<RWMol> t1(new RWMol(*m1));
    auto hsh1 =
        MolHash::MolHash(t1.get(), MolHash::HashFunction::HetAtomTautomer);
    std::unique_ptr<RWMol> t2(new RWMol(*m2));
    auto hsh2 =
        MolHash::MolHash(t2.get(), MolHash::HashFunction::HetAtomTautomer);
    CHECK(hsh1 == hsh2);
    CHECK(hsh1 ==
          "CC[C]([CH][C]1S[C]2[CH][C](C)[C](C)[CH][C]2N1CC([O])CS([O])([O])[O])"
          "[CH][C]1[Se][C]2[CH][CH][C](C)[CH][C]2N1CC_2_1");
  }
  SECTION("tautomers bug found in testing2") {
    auto m1 =
        "N/C(=N\\[N+](=O)[O-])NCCCCCCCC(=O)NC(CC(=O)OCc1ccccc1)C(=O)NCCCCN/C(N)=N/[N+](=O)[O-]"_smiles;
    REQUIRE(m1);
    auto m2 =
        "N/C(=N\\CCCCCCCC(=O)NC(CC(=O)OCc1ccccc1)C(=O)NCCCC/N=C(\\N)N[N+](=O)[O-])N[N+](=O)[O-]"_smiles;
    REQUIRE(m2);

    std::unique_ptr<RWMol> t1(new RWMol(*m1));
    auto hsh1 =
        MolHash::MolHash(t1.get(), MolHash::HashFunction::HetAtomTautomer);
    std::unique_ptr<RWMol> t2(new RWMol(*m2));
    auto hsh2 =
        MolHash::MolHash(t2.get(), MolHash::HashFunction::HetAtomTautomer);
    CHECK(hsh1 == hsh2);
    CHECK(hsh1 ==
          "[N][C]([N]CCCCCCC[C]([O])[N]C(C[C]([O])OC[C]1[CH][CH][CH][CH][CH]1)["
          "C]([O])[N]CCCC[N][C]([N])[N]N([O])[O])[N]N([O])[O]_8_0");
  }
}

TEST_CASE("Tautomers and chirality", "[molhash]") {
  SECTION("basics") {
    auto om = "C[C@H](C(=O)O)C(=O)[O-]"_smiles;
    REQUIRE(om);

    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh =
          MolHash::MolHash(m.get(), MolHash::HashFunction::CanonicalSmiles);
      CHECK(hsh == "C[C@@H](C(=O)[O-])C(=O)O");
    }
    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh =
          MolHash::MolHash(m.get(), MolHash::HashFunction::HetAtomTautomer);
      CHECK(hsh == "CC([C]([O])[O])[C]([O])[O]_1_-1");
    }
    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh =
          MolHash::MolHash(m.get(), MolHash::HashFunction::HetAtomProtomer);
      CHECK(hsh == "CC([C]([O])[O])[C]([O])[O]_2");
    }
  }
}

TEST_CASE("Molecular formula with fragments", "[molhash]") {
  SECTION("basics") {
    auto om = "CC(=O)[O-].C[N+](C)(C)C"_smiles;
    REQUIRE(om);
    auto hsh = MolHash::MolHash(om.get(), MolHash::HashFunction::MolFormula);
    CHECK(hsh == "C6H15NO2");
  }
}

TEST_CASE("Github issues", "[molhash]") {
  SECTION("Issue #4222: MolHash fails on non-standard valences") {
    SmilesParserParams p;
    p.sanitize = false;
    std::unique_ptr<RWMol> mol(SmilesToMol("C[Cl]C", p));
    REQUIRE(mol);
    auto hsh = MolHash::MolHash(mol.get(), MolHash::HashFunction::MolFormula);
    CHECK(hsh == "C2H6Cl");
  }
}

TEST_CASE("MolHash with CX extensions", "[molhash]") {
  SECTION("Tautomer") {
    auto mol =
        "C[C@@H](O)[C@@H](C)[C@@H](C)C[C@H](C1=CN=CN1)C1=CNC=N1 |o1:8,5,&1:1,3,r,c:11,18,t:9,15|"_smiles;
    REQUIRE(mol);

    {
      RWMol cp(*mol);
      auto hsh = MolHash::MolHash(&cp, MolHash::HashFunction::HetAtomTautomer);
      CHECK(
          hsh ==
          "C[C@H]([C@@H](C)[O])[C@@H](C)CC([C]1[CH][N][CH][N]1)[C]1[CH][N][CH][N]1_3_0");
    }
    {
      RWMol cp(*mol);

      auto hsh =
          MolHash::MolHash(&cp, MolHash::HashFunction::HetAtomTautomer, true);
      CHECK(
          hsh ==
          "C[C@H](CC([C]1[CH][N][CH][N]1)[C]1[CH][N][CH][N]1)[C@@H](C)[C@H](C)[O]_3_0 |o1:1,&1:14,16|");
    }
  }
  SECTION("no coordinates please") {
    auto mol = R"CTAB(
  Mrv2108 03032205502D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 15 16 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C -0.4657 -3.589 0 0 CFG=1
M  V30 2 C -0.4657 -2.049 0 0
M  V30 3 C 0.8679 -4.359 0 0
M  V30 4 C 2.2016 -3.589 0 0 CFG=2
M  V30 5 C 3.5353 -4.359 0 0
M  V30 6 C 4.9422 -3.7327 0 0
M  V30 7 N 5.9726 -4.8771 0 0
M  V30 8 C 5.2026 -6.2108 0 0
M  V30 9 N 3.6963 -5.8906 0 0
M  V30 10 C 2.2016 -2.049 0 0
M  V30 11 C 0.9557 -1.1438 0 0
M  V30 12 N 1.4316 0.3208 0 0
M  V30 13 C 2.9716 0.3208 0 0
M  V30 14 N 3.4475 -1.1438 0 0
M  V30 15 F -1.7994 -4.359 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 CFG=1
M  V30 2 1 1 3
M  V30 3 1 4 3 CFG=1
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 7
M  V30 7 2 7 8
M  V30 8 1 8 9
M  V30 9 1 5 9
M  V30 10 1 4 10
M  V30 11 2 10 11
M  V30 12 1 11 12
M  V30 13 1 12 13
M  V30 14 2 13 14
M  V30 15 1 10 14
M  V30 16 1 1 15
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEREL1 ATOMS=(1 1)
M  V30 MDLV30/STERAC1 ATOMS=(1 4)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);

    {
      RWMol cp(*mol);
      auto hsh =
          MolHash::MolHash(&cp, MolHash::HashFunction::HetAtomTautomer, true);
      CHECK(hsh ==
            "C[C@H](F)CC([C]1[CH][N][CH][N]1)[C]1[CH][N][CH][N]1_2_0 |o1:1|");
    }
  }

  SECTION("Mesomer") {
    auto mol = "C[C@H](F)C[C@@](C([NH-])=O)C([O-])=N |o1:1,&1:4|"_smiles;
    REQUIRE(mol);

    {
      RWMol cp(*mol);
      auto hsh = MolHash::MolHash(&cp, MolHash::HashFunction::Mesomer);
      CHECK(hsh == "C[C@H](F)C[C]([C]([NH])[O])[C]([NH])[O]_-2");
    }
    {
      RWMol cp(*mol);

      auto hsh = MolHash::MolHash(&cp, MolHash::HashFunction::Mesomer, true);
      CHECK(hsh == "C[C@H](F)C[C]([C]([NH])[O])[C]([NH])[O]_-2 |o1:1|");
    }
  }
  SECTION("Extended Murcko") {
    auto mol =
        "CC1=CC=CC=C1[C@@H](C[C@@H](C1CC1)C1CCC1)C1=CC=CC=C1O |o1:9,&1:7|"_smiles;
    REQUIRE(mol);

    {
      RWMol cp(*mol);
      auto hsh = MolHash::MolHash(&cp, MolHash::HashFunction::ExtendedMurcko);
      CHECK(hsh == "*c1ccccc1C(C[C@H](C1CCC1)C1CC1)c1ccccc1*");
    }
    {
      RWMol cp(*mol);

      auto hsh =
          MolHash::MolHash(&cp, MolHash::HashFunction::ExtendedMurcko, true);
      CHECK(hsh == "*c1ccccc1C(C[C@H](C1CCC1)C1CC1)c1ccccc1* |o1:9|");
    }
  }
  SECTION("Murcko") {
    auto mol =
        "CC1=CC=CC=C1[C@@H](C[C@@H](C1CC1)C1CCC1)C1=CC=CC=C1O |o1:9,&1:7|"_smiles;
    REQUIRE(mol);

    {
      RWMol cp(*mol);
      auto hsh = MolHash::MolHash(&cp, MolHash::HashFunction::MurckoScaffold);
      CHECK(hsh == "c1ccc(C(C[C@H](C2CCC2)C2CC2)c2ccccc2)cc1");
    }
    {
      RWMol cp(*mol);

      auto hsh =
          MolHash::MolHash(&cp, MolHash::HashFunction::MurckoScaffold, true);
      CHECK(hsh == "c1ccc(C(C[C@H](C2CCC2)C2CC2)c2ccccc2)cc1 |o1:6|");
    }
  }
  SECTION("Element") {
    auto mol =
        "C([C@@H](C1CC1)C1CCC1)[C@@H](C1CCCCC1)C1=CC=CC=C1 |o1:1,&1:9,c:21,23,t:19|"_smiles;
    REQUIRE(mol);

    {
      RWMol cp(*mol);
      auto hsh = MolHash::MolHash(&cp, MolHash::HashFunction::ElementGraph);
      CHECK(hsh == "C1CCC(C(C[C@H](C2CCC2)C2CC2)C2CCCCC2)CC1");
    }
    {
      RWMol cp(*mol);

      auto hsh =
          MolHash::MolHash(&cp, MolHash::HashFunction::ElementGraph, true);
      CHECK(hsh == "C1CCC(C(C[C@H](C2CCC2)C2CC2)C2CCCCC2)CC1 |o1:6|");
    }
  }
  SECTION("Anonymous") {
    auto mol =
        "C([C@@H](C1CC1)C1CCC1)[C@@H](C1CCCCC1)C1=CC=CC=N1 |o1:1,&1:9,c:21,23,t:19|"_smiles;
    REQUIRE(mol);

    {
      RWMol cp(*mol);
      auto hsh = MolHash::MolHash(&cp, MolHash::HashFunction::AnonymousGraph);
      CHECK(hsh == "*1***(*(**(*2***2)*2**2)*2*****2)**1");
    }
    {
      RWMol cp(*mol);

      auto hsh =
          MolHash::MolHash(&cp, MolHash::HashFunction::AnonymousGraph, true);
      CHECK(hsh == "*1***(*(**(*2***2)*2**2)*2*****2)**1");
    }
  }
}

TEST_CASE("tautomer v2") {
  SECTION("matches") {
    // pairs of {molecules with the same hash} {molecules with different hashes
    // from those}
    std::vector<std::pair<std::vector<std::string>, std::vector<std::string>>>
        data = {
            {{"CC=O", "C=CO"}, {}},
            {{"CCC=O", "CC=CO"}, {"C=CCO"}},
            {{"CC(=O)CC(=O)C", "C=C(O)CC(=O)C", "CC(=O)C=C(O)C",
              "C=C(O)C=C(O)C", "C=C(O)CC(O)=C"},
             {}},
            {{"CN=CCF", "CNC=CF"}, {"C=NCCF", "CNCCF"}},
            {{"CN=C(C)F", "CNC(=C)F"}, {"C=NC(C)F"}},
            {{"Cc1n[nH]cc1", "Cc1[nH][n]cc1", "CC1=NN=CC1", "CC1N=NCC=1"}, {}},
            {{"O=C1C=CC(=O)C=C1"}, {"Oc1ccc(O)cc1", "O=C1C=CC(O)C=C1"}},
            {{"CC(=O)CCC(=O)C", "CC(=O)CCC(O)=C", "C=C(O)CCC(O)=C"},
             {"CC(O)C=CC(=O)C", "CC(=O)C=CC(O)C"}},
            {{"c1ccccc1/C=C/c1ccccc1"},
             {"c1ccccc1/C=C\\c1ccccc1", "c1ccccc1C=Cc1ccccc1"}},
            // imine stereochemistry is lost:
            {{"CC/C=N/C", "CC/C=N\\C", "CCC=NC", "C/C=C/NC"}, {}},
            // but only when tautomers can happen:
            {{"FC(F)(F)/C(F)=N/C(F)(F)F"},
             {"FC(F)(F)/C(F)=N\\C(F)(F)F", "FC(F)(F)C(F)=NC(F)(F)F"}},
            {{"NC(=N)CC(=O)C", "NC(N)=CC(=O)C", "NC(=N)CC(O)=C",
              "NC(=N)C=C(O)C"},
             {}},
            {{"CC(=O)C=CC", "C=C(O)C=CC"}, {"CC(=O)CC=C"}},
            {{"N=C1N=CN(C)C2N=CNC=21", "NC1N=CN(C)C2=NC=NC2=1"}, {}},
            {{"S=C1N=CN=C2NC=NC12", "S=C1NC=NC2N=CNC1=2"}, {}},
            {{"CC1=CN=CN1", "CC1CN=CN=1"}, {}},
            {{"N1C(=O)NC(=O)C2C=NNC=21", "N1C(=O)NC(=O)C2=CNNC2=1",
              "N1C(=O)NC(=O)C2=CNNC2=1", "N1C(=O)NC(=O)C2CN=NC2=1"},
             {
                 "N1C(=O)NC(=O)C2CN=NC=21",
             }},
            // ---------------------------
            // more stereochemistry
            // ---------------------------
            {{"C[C@H](F)C=O", "C[C@@H](F)C=O", "CC(F)C=O"}, {}},
            {{"C[C@H](F)CC=O"}, {"C[C@@H](F)CC=O", "CC(F)CC=O"}},
            {{"C/C=C/O", "C/C=C\\O", "CC=CO"}, {}},
            {{"C/C=C/C=O", "C/C=C\\C=O", "CC=CC=O"}, {}},
            {{"C/C=C/CC=O"}, {"C/C=C\\CC=O", "CC=CCC=O"}},
            {{"C1C=CC=C2C=C(O)NC=12", "C1C=CC=C2CC(=O)NC=12"}, {}},
            // ---------------------------
            // some areas for potential improvement
            // ---------------------------
            // these two are tautomers, but the current algorithm does not catch
            // them
            // problematic because pyridine is recognized as tautomeric
            {{"c1ccccc1/C=C/c1ncccc1", "c1ccccc1/C=C\\c1ncccc1",
              "c1ccccc1C=Cc1ncccc1"},
             {}},
        };
    for (const auto &[same, diff] : data) {
      std::unique_ptr<RWMol> m{SmilesToMol(same[0])};
      REQUIRE(m);
      RWMol cp(*m);
      auto ref =
          MolHash::MolHash(&cp, MolHash::HashFunction::HetAtomTautomerv2);
      for (auto i = 1u; i < same.size(); ++i) {
        INFO(same[0] + "." + same[i]);
        std::unique_ptr<RWMol> m2{SmilesToMol(same[i])};
        REQUIRE(m2);
        RWMol cp(*m2);
        auto hsh =
            MolHash::MolHash(&cp, MolHash::HashFunction::HetAtomTautomerv2);
        CHECK(hsh == ref);
      }
      for (auto i = 0u; i < diff.size(); ++i) {
        INFO(same[0] + "." + diff[i]);
        std::unique_ptr<RWMol> m2{SmilesToMol(diff[i])};
        REQUIRE(m2);
        RWMol cp(*m2);
        auto hsh =
            MolHash::MolHash(&cp, MolHash::HashFunction::HetAtomTautomerv2);
        CHECK(hsh != ref);
      }
    }
  }

  SECTION("basics") {
    std::vector<std::tuple<std::string, std::string, std::string>> data = {
        {"C=O", "[CH2][O]_0_0", "[CH2]=[O]_0_0"},
        {"CC=O", "C[CH][O]_0_0", "[C]:[C]:[O]_4_0"},
        {"C=CO", "[CH2][CH][O]_1_0", "[C]:[C]:[O]_4_0"},
        {"c1ccccc1", "[CH]1[CH][CH][CH][CH][CH]1_0_0",
         "[cH]1:[cH]:[cH]:[cH]:[cH]:[cH]:1_0_0"},
        {"n1ccccc1", "[CH]1[CH][CH][N][CH][CH]1_0_0",
         "[C]1:[C]:[C]:[N]:[C]:[C]:1_5_0"},
        {"Nc1ccccc1", "[N][C]1[CH][CH][CH][CH][CH]1_2_0",
         "[N]:[C]1:[C]:[C]:[C]:[C]:[C]:1_7_0"},
        {"C=COC", "[CH2][CH]OC_0_0", "[CH2]=[CH]-[O]-[CH3]_0_0"},
        {"CC(C)(C)C=O", "CC(C)(C)[CH][O]_0_0",
         "[CH3]-[C](-[CH3])(-[CH3])-[CH]=[O]_0_0"},
        {"CC(C)=CO", "C[C](C)[CH][O]_1_0", "[CH3]-[C](-[CH3]):[C]:[O]_2_0"},
        {"COC=O", "CO[CH][O]_0_0", "[CH3]-[O]:[C]:[O]_1_0"},
        {"CNC=O", "C[N][CH][O]_1_0", "[CH3]-[N]:[C]:[O]_2_0"},
        {"CN(C)C=O", "CN(C)[CH][O]_0_0", "[CH3]-[N](-[CH3]):[C]:[O]_1_0"},
        {"CC(C)(C)NC=O", "CC(C)(C)[N][CH][O]_1_0",
         "[CH3]-[C](-[CH3])(-[CH3])-[N]:[C]:[O]_2_0"},
        {"CC(C)=O", "C[C](C)[O]_0_0", "[C]:[C](:[C]):[O]_6_0"},
        {"C=C(C)O", "[CH2][C](C)[O]_1_0", "[C]:[C](:[C]):[O]_6_0"},
        {"N1CCC1", "C1C[N]C1_1_0", "[CH2]1-[CH2]-[NH]-[CH2]-1_0_0"},
        {"CC=CC(=O)C", "C[CH][CH][C](C)[O]_0_0",
         "[C]:[C](:[O]):[C]:[C]-[CH3]_5_0"},
        {"N1C=CCC(F)C1", "FC1C[CH][CH][N]C1_1_0",
         "[F]-[CH]1-[CH2]-[C]:[C]:[N]-[CH2]-1_3_0"},
        {"CCC=C(O)C", "CC[CH][C](C)[O]_1_0",
         "[C]:[C](:[O]):[C]-[CH2]-[CH3]_5_0"},
        {"CCCC(=O)C", "CCC[C](C)[O]_0_0", "[C]:[C](:[O]):[C]-[CH2]-[CH3]_5_0"},
        {"CCCC(O)=C", "[CH2][C]([O])CCC_1_0",
         "[C]:[C](:[O]):[C]-[CH2]-[CH3]_5_0"},
        {"C=CCC(O)C", "C=CCC(C)[O]_1_0",
         "[CH2]=[CH]-[CH2]-[CH](-[CH3])-[OH]_0_0"},
        {"C=NC(=O)C", "[CH2][N][C](C)[O]_0_0", "[C]:[N]:[C](:[C]):[O]_5_0"},
        {"C=NC(O)=C", "[CH2][N][C]([CH2])[O]_1_0", "[C]:[N]:[C](:[C]):[O]_5_0"},
        {"CC(=O)CC(=O)C", "C[C]([O])C[C](C)[O]_0_0",
         "[C]:[C](:[O]):[C]:[C](:[C]):[O]_8_0"},
        {"CC(=O)C=C(O)C", "C[C]([O])[CH][C](C)[O]_1_0",
         "[C]:[C](:[O]):[C]:[C](:[C]):[O]_8_0"},
        {"C=C(O)C=C(O)C", "[CH2][C]([O])[CH][C](C)[O]_2_0",
         "[C]:[C](:[O]):[C]:[C](:[C]):[O]_8_0"},
        {"C=C(O)CC(O)=C", "[CH2][C]([O])C[C]([CH2])[O]_2_0",
         "[C]:[C](:[O]):[C]:[C](:[C]):[O]_8_0"},
        {"CC(=O)CCC(=O)C", "C[C]([O])CC[C](C)[O]_0_0",
         "[C]:[C](:[O]):[C]-[C]:[C](:[C]):[O]_10_0"},
        {"CC(=O)C=CC(=O)C", "C[C]([O])[CH][CH][C](C)[O]_0_0",
         "[C]:[C](:[O]):[C]:[C]:[C](:[C]):[O]_8_0"},
    };
    for (const auto &tpl : data) {
      INFO(std::get<0>(tpl));
      std::unique_ptr<RWMol> m{SmilesToMol(std::get<0>(tpl))};
      REQUIRE(m);
      {
        RWMol cp(*m);
        auto hsh1 =
            MolHash::MolHash(&cp, MolHash::HashFunction::HetAtomTautomer);
        CHECK(hsh1 == std::get<1>(tpl));
      }
      {
        RWMol cp(*m);
        auto hsh2 =
            MolHash::MolHash(&cp, MolHash::HashFunction::HetAtomTautomerv2);
        CHECK(hsh2 == std::get<2>(tpl));
      }
    }
  }
}

TEST_CASE("tautomer hash problem cases") {
#if 1
  SECTION("sulfur problem") {
    auto m = R"CTAB(
     RDKit          2D

 22 24  0  0  0  0  0  0  0  0999 V2000
   -1.3203  -10.1153    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9127  -10.8290    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0927  -10.8317    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3205  -10.1216    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0921   -9.4072    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9107   -9.4080    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1477  -10.1129    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5593   -9.3951    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5637  -10.8282    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3911  -10.8258    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8786  -10.1547    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6663  -10.4082    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6688  -11.2356    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8825  -11.4935    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.1480  -10.1229    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5606  -10.8402    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5629   -9.4070    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.3903   -9.4084    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8743  -10.0739    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.6616   -9.8194    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.6630   -8.9920    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8765   -8.7351    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0
  5  6  2  0
  6  1  1  0
  1  2  2  0
 11 12  1  0
 12 13  2  0
 13 14  1  0
 14 10  2  0
  1  7  1  0
  4 15  1  0
  3  4  2  0
 15 16  2  0
  7  8  2  0
 15 17  1  0
 17 18  1  0
 18 19  2  0
  7  9  1  0
  4  5  1  0
  9 10  1  0
 10 11  1  0
 19 20  1  0
 20 21  2  0
 21 22  1  0
 22 18  1  0
M  END
)CTAB"_ctab;
    REQUIRE(m);
    auto hsh =
        MolHash::MolHash(m.get(), MolHash::HashFunction::HetAtomTautomerv2);
    CHECK(hsh.find("[s]") == std::string::npos);
  }
  SECTION("atom order") {
    auto m = R"CTAB(
     RDKit          2D

 17 18  0  0  0  0  0  0  0  0999 V2000
   12.9442  -15.7431    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   10.7279  -14.6717    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.7266  -15.4976    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.4381  -15.9096    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.4361  -14.2599    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.1526  -14.6678    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.1578  -15.4930    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   13.4251  -15.0725    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.9359  -14.4079    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0130  -14.2612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0123  -13.4380    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.3004  -14.6734    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   14.2501  -15.0676    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   14.6661  -15.7779    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   14.6573  -14.3521    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   15.4893  -15.7728    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   15.9055  -16.4831    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  8  1  0
  8  9  2  0
  9  6  1  0
  6  5  2  0
  5  2  1  0
  6  7  1  0
 10 11  1  0
 10 12  2  0
  2 10  1  0
  3  4  1  0
  4  7  2  0
 13 14  1  0
 13 15  2  0
  8 13  1  0
  2  3  2  0
 14 16  1  0
  7  1  1  0
 16 17  1  0
M  END

 > <chembl_id>
 CHEMBL503643

 > <chembl_pref_name>
 None
 )CTAB"_ctab;
    REQUIRE(m);
    std::vector<std::string> row = {"CCOC(=O)c1cc2cc(C(=O)O)ccc2[nH]1",
                                    "CCOC(=O)c1cc2cc(ccc2[nH]1)C(O)=O",
                                    "O(C(=O)c1cc2c(ccc(c2)C(=O)O)[nH]1)CC"};
    auto hsh =
        MolHash::MolHash(m.get(), MolHash::HashFunction::HetAtomTautomerv2);
    for (const auto &smi : row) {
      std::unique_ptr<RWMol> mi{SmilesToMol(smi)};
      REQUIRE(mi);
      auto hshi =
          MolHash::MolHash(mi.get(), MolHash::HashFunction::HetAtomTautomerv2);
      CHECK(hsh == hshi);
      break;
    }
  }
#endif
  SECTION("atom order 2") {
    // {molecules with the same hash}
    std::vector<std::vector<std::string>> data = {
        {"Ic1ccc(Cn2cc[n+](Cc3ccc(I)cc3)c2)cc1",
         "c1[n+](Cc2ccc(cc2)I)cn(Cc2ccc(cc2)I)c1"},
        {"CN1C(=O)c2cccnc2NC1c1cccnc1", "c1ncccc1C1N(C)C(c2cccnc2N1)=O",
         "c1nc2c(cc1)C(N(C)C(N2)c1cccnc1)=O"},
        {"CC(=O)OCC1=C(C(=O)[O-])N2C(=O)C(=C(Br)Br)[C@H]2S(=O)(=O)C1",
         "BrC(=C1C(=O)N2C(=C(COC(=O)C)CS([C@@H]21)(=O)=O)C(=O)[O-])Br"}};
    for (const auto &same : data) {
      std::unique_ptr<RWMol> m{SmilesToMol(same[0])};
      REQUIRE(m);
      RWMol cp(*m);
      auto ref =
          MolHash::MolHash(&cp, MolHash::HashFunction::HetAtomTautomerv2);
      for (auto i = 1u; i < same.size(); ++i) {
        INFO(same[0] + "->" + same[i]);
        std::unique_ptr<RWMol> m2{SmilesToMol(same[i])};
        REQUIRE(m2);
        RWMol cp(*m2);
        auto hsh =
            MolHash::MolHash(&cp, MolHash::HashFunction::HetAtomTautomerv2);
        CHECK(ref == hsh);
      }
    }
  }
}

TEST_CASE("GitHub Issue #6505") {
  const auto m = "CCCCCC[NH3+] |SgD:6:lambda max:230:=:nm::|"_smiles;
  REQUIRE(m);
  REQUIRE(getSubstanceGroups(*m).size() == 1);

  const auto use_cx_smiles = true;

  SECTION("Do not skip any CX flags") {
    const auto cx_to_skip = SmilesWrite::CXSmilesFields::CX_NONE;
    const auto hsh1 =
        MolHash::MolHash(m.get(), MolHash::HashFunction::HetAtomTautomerv2,
                         use_cx_smiles, cx_to_skip);
    CHECK(
        hsh1 ==
        "[CH3]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[NH3+]_0_0 |SgD:6:lambda max:230:=:nm::|");
  }

  SECTION("Strip all CX flags") {
    const auto cx_to_skip = SmilesWrite::CXSmilesFields::CX_ALL;
    const auto hsh2 =
        MolHash::MolHash(m.get(), MolHash::HashFunction::HetAtomTautomerv2,
                         use_cx_smiles, cx_to_skip);
    CHECK(hsh2 == "[CH3]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[NH3+]_0_0");
  }
}

TEST_CASE("Github Issue #6855 MakeScaffoldGeneric isotope removal") {
  SECTION("Extended Murcko") {
    auto mol = "[235U]C1CC1"_smiles;
    REQUIRE(mol);
    {
      RWMol cp(*mol);
      auto hsh = MolHash::MolHash(&cp, MolHash::HashFunction::ExtendedMurcko);
      CHECK(hsh == "*C1CC1");
    }
  }
  SECTION("Anonymous") {
    auto mol = "[235U]1CC1"_smiles;
    REQUIRE(mol);

    {
      RWMol cp(*mol);
      auto hsh = MolHash::MolHash(&cp, MolHash::HashFunction::AnonymousGraph);
      CHECK(hsh == "*1**1");
    }
  }
}

TEST_CASE("Github Issue #6472 non-matching element and anononymous graph") {
  SECTION("Element graph test1") {
    auto mol = "C1COC(C1)C1=NC=NC=C1"_smiles;
    REQUIRE(mol);
    {
      RWMol cp(*mol);
      auto hsh = MolHash::MolHash(&cp, MolHash::HashFunction::ElementGraph);
      CHECK(hsh == "C1COC(C2CCNCN2)C1");
    }
  }
  SECTION("Element graph test2") {
    auto mol = "C1CC(N=CN1)C1=CC=CO1"_smiles;
    REQUIRE(mol);
    {
      RWMol cp(*mol);
      auto hsh = MolHash::MolHash(&cp, MolHash::HashFunction::ElementGraph);
      CHECK(hsh == "C1COC(C2CCNCN2)C1");
    }
  }
  SECTION("Anonymous graph test 1") {
    auto mol = "C1COC(C1)C1=NC=NC=C1"_smiles;
    REQUIRE(mol);

    {
      RWMol cp(*mol);
      auto hsh = MolHash::MolHash(&cp, MolHash::HashFunction::AnonymousGraph);
      CHECK(hsh == "*1***(*2****2)**1");
    }
  }
  SECTION("Anonymous graph test 2") {
    auto mol = "C1CC(N=CN1)C1=CC=CO1"_smiles;
    REQUIRE(mol);

    {
      RWMol cp(*mol);
      auto hsh = MolHash::MolHash(&cp, MolHash::HashFunction::AnonymousGraph);
      CHECK(hsh == "*1***(*2****2)**1");
    }
  }
}

TEST_CASE("tautomer overreach") {
  SECTION("as reported") {
    auto mol = "C1=CN(C[C@H]2CNCCO2)N=C1"_smiles;
    REQUIRE(mol);
    {
      RWMol cp(*mol);
      auto hsh =
          MolHash::MolHash(&cp, MolHash::HashFunction::HetAtomTautomerv2);
      CHECK(
          hsh ==
          "[C]1:[C]:[N]:[N](-[CH2]-[C@H]2-[CH2]-[NH]-[CH2]-[CH2]-[O]-2):[C]:1_3_0");
    }
  }
  SECTION("dbw example") {
    auto mol = "c1cccn1C[C@H](C)COC"_smiles;
    REQUIRE(mol);
    {
      RWMol cp(*mol);
      auto hsh =
          MolHash::MolHash(&cp, MolHash::HashFunction::HetAtomTautomerv2);
      CHECK(
          hsh ==
          "[CH3]-[O]-[CH2]-[C@@H](-[CH3])-[CH2]-[n]1:[cH]:[cH]:[cH]:[cH]:1_0_0");
    }
  }
}

TEST_CASE("HetAtomProtomerv2") {
  SECTION("matches") {
    // pairs of {molecules with the same hash} {molecules with different hashes
    // from those}
    std::vector<std::pair<std::vector<std::string>, std::vector<std::string>>>
        data = {
            // example from the NextMove documentation
            {{"Cc1c[nH]cn1", "Cc1cnc[nH]1", "Cc1c[nH]c[nH+]1"}, {}},
            {{"CC=CO", "CCC=O"}, {"C=CCO"}},

        };
    for (const auto &[same, diff] : data) {
      std::unique_ptr<RWMol> m{SmilesToMol(same[0])};
      REQUIRE(m);
      RWMol cp(*m);
      auto ref =
          MolHash::MolHash(&cp, MolHash::HashFunction::HetAtomProtomerv2);
      for (auto i = 1u; i < same.size(); ++i) {
        INFO(same[0] + "->" + same[i]);
        std::unique_ptr<RWMol> m2{SmilesToMol(same[i])};
        REQUIRE(m2);
        RWMol cp(*m2);
        auto hsh =
            MolHash::MolHash(&cp, MolHash::HashFunction::HetAtomProtomerv2);
        CHECK(hsh == ref);
      }
      for (auto i = 0u; i < diff.size(); ++i) {
        INFO(same[0] + "->" + diff[i]);
        std::unique_ptr<RWMol> m2{SmilesToMol(diff[i])};
        REQUIRE(m2);
        RWMol cp(*m2);
        auto hsh =
            MolHash::MolHash(&cp, MolHash::HashFunction::HetAtomProtomerv2);
        CHECK(hsh != ref);
      }
    }
  }
}

TEST_CASE("overreach with v2 tautomer hashes and imines") {
  SECTION("basics") {
    std::vector<std::string> smileses = {"C[C@H](F)NC1=CCCCC1",
                                         "C[C@H](F)N=C1CCCCC1"};
    for (const auto &smiles : smileses) {
      auto m = v2::SmilesParse::MolFromSmiles(smiles);
      REQUIRE(m);
      auto hsh =
          MolHash::MolHash(m.get(), MolHash::HashFunction::HetAtomTautomerv2);
      CHECK(hsh ==
            "[CH3]-[C@H](-[F])-[N]:[C]1:[C]-[CH2]-[CH2]-[CH2]-[C]:1_4_0");
    }
  }
}