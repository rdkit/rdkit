//
//  Copyright (C) 2019-2022 Greg Landrum
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
        "C[C@@H](O)[C@@H](C)[C@@H](C)C[C@H](C1=CN=CN1)C1=CNC=N1 "
        "|o1:8,5,&1:1,3,r,c:11,18,t:9,15|"_smiles;
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
          "C[C@H]([C@@H](C)[O])[C@@H](C)CC([C]1[CH][N][CH][N]1)[C]1[CH][N][CH][N]1_3_0 |o1:5,&1:1,2|");
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
