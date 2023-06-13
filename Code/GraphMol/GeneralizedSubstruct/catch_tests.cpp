//
//  Copyright (c) 2023, Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Tests of the generalized substructure searching code
//

#include "catch.hpp"

#include <tuple>
#include <utility>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/GeneralizedSubstruct/XQMol.h>
#include <GraphMol/TautomerQuery/TautomerQuery.h>
#include <GraphMol/MolEnumerator/MolEnumerator.h>

using namespace RDKit;
using namespace RDKit::GeneralizedSubstruct;

TEST_CASE("molecule basics") {
  auto mol = "Cc1n[nH]c(F)c1"_smarts;
  REQUIRE(mol);
  ExtendedQueryMol xqm = std::make_unique<RWMol>(*mol);
  SECTION("substructure matching and serialization") {
    ExtendedQueryMol xqm2(xqm.toBinary());
    ExtendedQueryMol xqm3(xqm.toJSON(), true);
    for (const auto xq : {&xqm, &xqm2, &xqm3}) {
      CHECK(SubstructMatch(*"CCc1n[nH]c(F)c1"_smiles, *xq).size() == 1);
      CHECK(SubstructMatch(*"CCc1[nH]nc(F)c1"_smiles, *xq).empty());
      CHECK(hasSubstructMatch(*"CCc1n[nH]c(F)c1"_smiles, *xq));
      CHECK(!hasSubstructMatch(*"CCc1[nH]nc(F)c1"_smiles, *xq));
    }
  }
}

TEST_CASE("enumeration basics") {
  auto mol = "COCC |LN:1:1.3|"_smiles;
  REQUIRE(mol);
  ExtendedQueryMol xqm =
      std::make_unique<MolBundle>(MolEnumerator::enumerate(*mol));
  SECTION("substructure matching and serialization") {
    ExtendedQueryMol xqm2(xqm.toBinary());
    ExtendedQueryMol xqm3(xqm.toJSON(), true);

    for (const auto xq : {&xqm, &xqm2, &xqm3}) {
      CHECK(SubstructMatch(*"COCC"_smiles, *xq).size() == 1);
      CHECK(SubstructMatch(*"COOCC"_smiles, *xq).size() == 1);
      CHECK(SubstructMatch(*"COOOCC"_smiles, *xq).size() == 1);
      CHECK(SubstructMatch(*"COOOOCC"_smiles, *xq).empty());
    }
  }
}

TEST_CASE("result counts") {
  auto mol = "COC |LN:1:1.3|"_smiles;
  REQUIRE(mol);
  ExtendedQueryMol xqm =
      std::make_unique<MolBundle>(MolEnumerator::enumerate(*mol));
  SECTION("substructure matching and serialization") {
    ExtendedQueryMol xqm2(xqm.toBinary());
    ExtendedQueryMol xqm3(xqm.toJSON(), true);
    SubstructMatchParameters ps;
    ps.uniquify = false;
    for (const auto xq : {&xqm, &xqm2, &xqm3}) {
      CHECK(SubstructMatch(*"COCC"_smiles, *xq, ps).size() == 2);
      CHECK(SubstructMatch(*"COOCC"_smiles, *xq, ps).size() == 2);
      CHECK(SubstructMatch(*"COOOCC"_smiles, *xq, ps).size() == 2);
      CHECK(SubstructMatch(*"COOOOCC"_smiles, *xq, ps).empty());
    }
  }
}

TEST_CASE("tautomer basics") {
  auto mol = "Cc1n[nH]c(F)c1"_smiles;
  REQUIRE(mol);
  ExtendedQueryMol xqm =
      std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*mol));
  const TautomerQuery &otq =
      *std::get<ExtendedQueryMol::TautomerQuery_T>(xqm.xqmol);

  SECTION("substructure matching and serialization") {
    ExtendedQueryMol xqm2(xqm.toBinary());
    ExtendedQueryMol xqm3(xqm.toJSON(), true);
    for (const auto xq : {&xqm, &xqm2, &xqm3}) {
      CHECK(SubstructMatch(*"CCc1n[nH]c(F)c1"_smiles, *xq).size() == 1);
      CHECK(SubstructMatch(*"CCc1[nH]nc(F)c1"_smiles, *xq).size() == 1);
      CHECK(SubstructMatch(*"CCc1[nH]ncc1"_smiles, *xq).empty());
    }
  }
}

TEST_CASE("tautomer bundle basics") {
  auto mol1 = "Cc1n[nH]c(F)c1"_smiles;
  REQUIRE(mol1);
  auto mol2 = "Cc1n[nH]cc1F"_smiles;
  REQUIRE(mol2);
  std::vector<std::unique_ptr<TautomerQuery>> tbndl;
  tbndl.emplace_back(
      std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*mol1)));
  tbndl.emplace_back(
      std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*mol2)));
  ExtendedQueryMol xqm =
      std::make_unique<std::vector<std::unique_ptr<TautomerQuery>>>(
          std::move(tbndl));
  SECTION("substructure matching and serialization") {
    ExtendedQueryMol xqm2(xqm.toBinary());
    ExtendedQueryMol xqm3(xqm.toJSON(), true);
    for (const auto xq : {&xqm, &xqm2, &xqm3}) {
      CHECK(SubstructMatch(*"CCc1n[nH]c(F)c1"_smiles, *xq).size() == 1);
      CHECK(SubstructMatch(*"CCc1[nH]nc(F)c1"_smiles, *xq).size() == 1);
      CHECK(SubstructMatch(*"CCc1[nH]ncc1F"_smiles, *xq).size() == 1);
      CHECK(SubstructMatch(*"CCc1n[nH]cc1F"_smiles, *xq).size() == 1);
      CHECK(SubstructMatch(*"CCc1[nH]ncc1"_smiles, *xq).empty());
    }
  }
}

TEST_CASE("createExtendedQueryMol") {
  SECTION("RWMol") {
    auto mol = "COCC"_smiles;
    REQUIRE(mol);
    auto xqm = createExtendedQueryMol(*mol);
    CHECK(std::holds_alternative<ExtendedQueryMol::RWMol_T>(xqm.xqmol));
    CHECK(SubstructMatch(*"COCC"_smiles, xqm).size() == 1);
    CHECK(SubstructMatch(*"COOCC"_smiles, xqm).empty());
  }
  SECTION("MolBundle") {
    auto mol = "COCC |LN:1:1.3|"_smiles;
    REQUIRE(mol);
    auto xqm = createExtendedQueryMol(*mol);
    CHECK(std::holds_alternative<ExtendedQueryMol::MolBundle_T>(xqm.xqmol));
    CHECK(SubstructMatch(*"COCC"_smiles, xqm).size() == 1);
    CHECK(SubstructMatch(*"COOCC"_smiles, xqm).size() == 1);
    CHECK(SubstructMatch(*"COOOCC"_smiles, xqm).size() == 1);
    CHECK(SubstructMatch(*"COOOOCC"_smiles, xqm).empty());
  }
  SECTION("TautomerQuery") {
    auto mol1 = "Cc1n[nH]c(F)c1"_smiles;
    REQUIRE(mol1);
    auto xqm = createExtendedQueryMol(*mol1);
    CHECK(std::holds_alternative<ExtendedQueryMol::TautomerQuery_T>(xqm.xqmol));
    CHECK(SubstructMatch(*"CCc1n[nH]c(F)c1"_smiles, xqm).size() == 1);
    CHECK(SubstructMatch(*"CCc1[nH]nc(F)c1"_smiles, xqm).size() == 1);
    CHECK(SubstructMatch(*"CCc1[nH]ncc1"_smiles, xqm).empty());
  }
  SECTION("TautomerBundle") {
    auto mol1 = "COCc1n[nH]c(F)c1 |LN:1:1.3|"_smiles;
    REQUIRE(mol1);
    auto xqm = createExtendedQueryMol(*mol1);
    CHECK(
        std::holds_alternative<ExtendedQueryMol::TautomerBundle_T>(xqm.xqmol));
    CHECK(SubstructMatch(*"COCc1n[nH]c(F)c1"_smiles, xqm).size() == 1);
    CHECK(SubstructMatch(*"COOCc1n[nH]c(F)c1"_smiles, xqm).size() == 1);
    CHECK(SubstructMatch(*"COCc1[nH]nc(F)c1"_smiles, xqm).size() == 1);
    CHECK(SubstructMatch(*"COOCc1[nH]nc(F)c1"_smiles, xqm).size() == 1);
  }
}

TEST_CASE("test SRUs") {
  SECTION("basics") {
    auto mol = "FCN(CC)-* |Sg:n:2,5:1-2:ht|"_smiles;
    REQUIRE(mol);
    auto xqm = createExtendedQueryMol(*mol);
    CHECK(std::holds_alternative<ExtendedQueryMol::MolBundle_T>(xqm.xqmol));
    // as of v2023.03.1 the SRU enumerator ignores repeat counts, so
    // we won't test limits here.
    CHECK(SubstructMatch(*"FCN(C)CC"_smiles, xqm).size() == 1);
    CHECK(SubstructMatch(*"FCN(O)N(C)CC"_smiles, xqm).size() == 1);
  }
}

// there's some redundancy in testing with what's above, but duplicating tests
// isn't a terrible thing
TEST_CASE("adjustQueryProperties") {
  SECTION("RWMol") {
    auto mol = "COC1CC1"_smiles;
    REQUIRE(mol);
    auto xqm1 = createExtendedQueryMol(*mol);
    auto xqm2 = createExtendedQueryMol(*mol, true);
    MolOps::AdjustQueryParameters aqps =
        MolOps::AdjustQueryParameters::noAdjustments();
    aqps.makeAtomsGeneric = true;
    auto xqm3 = createExtendedQueryMol(*mol, true, aqps);
    CHECK(std::holds_alternative<ExtendedQueryMol::RWMol_T>(xqm1.xqmol));
    CHECK(std::holds_alternative<ExtendedQueryMol::RWMol_T>(xqm2.xqmol));
    CHECK(std::holds_alternative<ExtendedQueryMol::RWMol_T>(xqm3.xqmol));
    CHECK(SubstructMatch(*"COC1CC1"_smiles, xqm1).size() == 1);
    CHECK(SubstructMatch(*"COC1CC1"_smiles, xqm2).size() == 1);
    CHECK(SubstructMatch(*"COC1CC1"_smiles, xqm3).size() == 1);
    CHECK(SubstructMatch(*"COC1C(C)C1"_smiles, xqm1).size() == 1);
    CHECK(SubstructMatch(*"COC1C(C)C1"_smiles, xqm2).empty());
    CHECK(SubstructMatch(*"COC1C(C)C1"_smiles, xqm3).size() == 1);
    CHECK(SubstructMatch(*"COC1OC1"_smiles, xqm1).empty());
    CHECK(SubstructMatch(*"COC1OC1"_smiles, xqm2).empty());
    CHECK(SubstructMatch(*"COC1OC1"_smiles, xqm3).size() == 1);
  }
  SECTION("MolBundle") {
    auto mol = "COCC |LN:1:1.3|"_smiles;
    REQUIRE(mol);
    auto xqm1 = createExtendedQueryMol(*mol);
    MolOps::AdjustQueryParameters aqps =
        MolOps::AdjustQueryParameters::noAdjustments();
    aqps.makeBondsGeneric = true;
    auto xqm2 = createExtendedQueryMol(*mol, true, aqps);
    CHECK(std::holds_alternative<ExtendedQueryMol::MolBundle_T>(xqm1.xqmol));
    CHECK(std::holds_alternative<ExtendedQueryMol::MolBundle_T>(xqm2.xqmol));
    CHECK(SubstructMatch(*"COC=C"_smiles, xqm1).empty());
    CHECK(SubstructMatch(*"COOC=C"_smiles, xqm1).empty());
    CHECK(SubstructMatch(*"COC=C"_smiles, xqm2).size() == 1);
    CHECK(SubstructMatch(*"COOC=C"_smiles, xqm2).size() == 1);
  }
  SECTION("TautomerQuery") {
    auto mol1 = "Cc1n[nH]cc1"_smiles;
    REQUIRE(mol1);
    auto xqm1 = createExtendedQueryMol(*mol1);
    auto xqm2 = createExtendedQueryMol(*mol1, true);

    CHECK(
        std::holds_alternative<ExtendedQueryMol::TautomerQuery_T>(xqm1.xqmol));
    CHECK(
        std::holds_alternative<ExtendedQueryMol::TautomerQuery_T>(xqm2.xqmol));
    CHECK(SubstructMatch(*"CCc1[nH]ncc1"_smiles, xqm1).size() == 1);
    CHECK(SubstructMatch(*"CCc1[nH]ncc1"_smiles, xqm2).size() == 1);
    CHECK(SubstructMatch(*"CCc1n[nH]c(F)c1"_smiles, xqm1).size() == 1);
    CHECK(SubstructMatch(*"CCc1[nH]nc(F)c1"_smiles, xqm1).size() == 1);
    CHECK(SubstructMatch(*"CCc1n[nH]c(F)c1"_smiles, xqm2).empty());
    CHECK(SubstructMatch(*"CCc1[nH]nc(F)c1"_smiles, xqm2).empty());
  }
  SECTION("TautomerBundle") {
    auto mol1 = "COCc1n[nH]cc1 |LN:1:1.3|"_smiles;
    REQUIRE(mol1);
    auto xqm1 = createExtendedQueryMol(*mol1);
    auto xqm2 = createExtendedQueryMol(*mol1, true);
    CHECK(
        std::holds_alternative<ExtendedQueryMol::TautomerBundle_T>(xqm1.xqmol));
    CHECK(SubstructMatch(*"COCc1n[nH]cc1"_smiles, xqm1).size() == 1);
    CHECK(SubstructMatch(*"COOCc1[nH]ncc1"_smiles, xqm1).size() == 1);
    CHECK(SubstructMatch(*"COCc1n[nH]cc1"_smiles, xqm2).size() == 1);
    CHECK(SubstructMatch(*"COOCc1[nH]ncc1"_smiles, xqm2).size() == 1);
    CHECK(SubstructMatch(*"COCc1n[nH]c(F)c1"_smiles, xqm1).size() == 1);
    CHECK(SubstructMatch(*"COOCc1[nH]nc(F)c1"_smiles, xqm1).size() == 1);
    CHECK(SubstructMatch(*"COCc1n[nH]c(F)c1"_smiles, xqm2).empty());
    CHECK(SubstructMatch(*"COOCc1[nH]nc(F)c1"_smiles, xqm2).empty());
  }
}
