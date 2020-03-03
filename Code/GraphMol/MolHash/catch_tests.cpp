//
//  Copyright (C) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
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
