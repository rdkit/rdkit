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

#include <catch2/catch_all.hpp>

#include <tuple>
#include <utility>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/GeneralizedSubstruct/XQMol.h>
#include <GraphMol/TautomerQuery/TautomerQuery.h>
#include <GraphMol/MolEnumerator/MolEnumerator.h>
#include <GraphMol/GenericGroups/GenericGroups.h>

using namespace RDKit;
using namespace RDKit::GeneralizedSubstruct;

bool fingerprintsMatch(const ROMol& target, const ExtendedQueryMol& xqm) {
  const auto queryFingerprint = xqm.patternFingerprintQuery();
  const auto targetFingerprint = patternFingerprintTargetMol(target);
  CHECK(queryFingerprint->getNumOnBits() > 0);
  CHECK(targetFingerprint->getNumOnBits() > 0);
  const auto match = AllProbeBitsMatch(*queryFingerprint, *targetFingerprint);
  return match;
}

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
      CHECK(fingerprintsMatch(*"CCc1n[nH]c(F)c1"_smiles, *xq));
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
      CHECK(fingerprintsMatch(*"COCC"_smiles, *xq));
      CHECK(fingerprintsMatch(*"COOCC"_smiles, *xq));
      CHECK(fingerprintsMatch(*"COOOCC"_smiles, *xq));
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
      CHECK(fingerprintsMatch(*"COCC"_smiles, *xq));
      CHECK(fingerprintsMatch(*"COOCC"_smiles, *xq));
      CHECK(fingerprintsMatch(*"COOOCC"_smiles, *xq));
    }
  }
}

TEST_CASE("tautomer basics") {
  auto mol = "Cc1n[nH]c(F)c1"_smiles;
  REQUIRE(mol);
  ExtendedQueryMol xqm =
      std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*mol));

  SECTION("substructure matching and serialization") {
    ExtendedQueryMol xqm2(xqm.toBinary());
    ExtendedQueryMol xqm3(xqm.toJSON(), true);
    for (const auto xq : {&xqm, &xqm2, &xqm3}) {
      CHECK(SubstructMatch(*"CCc1n[nH]c(F)c1"_smiles, *xq).size() == 1);
      CHECK(SubstructMatch(*"CCc1[nH]nc(F)c1"_smiles, *xq).size() == 1);
      CHECK(SubstructMatch(*"CCc1[nH]ncc1"_smiles, *xq).empty());
      CHECK(fingerprintsMatch(*"CCc1n[nH]c(F)c1"_smiles, *xq));
      CHECK(fingerprintsMatch(*"CCc1[nH]nc(F)c1"_smiles, *xq));
      CHECK(!fingerprintsMatch(*"CCc1[nH]ncc1"_smiles, *xq));
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
      CHECK(fingerprintsMatch(*"CCc1n[nH]c(F)c1"_smiles, *xq));
      CHECK(fingerprintsMatch(*"CCc1[nH]nc(F)c1"_smiles, *xq));
      CHECK(fingerprintsMatch(*"CCc1[nH]ncc1F"_smiles, *xq));
      CHECK(fingerprintsMatch(*"CCc1n[nH]cc1F"_smiles, *xq));
      CHECK(!fingerprintsMatch(*"CCc1[nH]ncc1"_smiles, *xq));
    }
  }
}

TEST_CASE("createExtendedQueryMol and copy ctors") {
  SECTION("RWMol") {
    auto mol = "COCC"_smiles;
    REQUIRE(mol);
    auto txqm = createExtendedQueryMol(*mol);
    ExtendedQueryMol xqm1(txqm);
    ExtendedQueryMol xqm2(std::make_unique<RWMol>(*mol));
    xqm2 = txqm;

    for (const auto &xqm : {txqm, xqm1, xqm2}) {
      CHECK(std::holds_alternative<ExtendedQueryMol::RWMol_T>(xqm.xqmol));
      CHECK(SubstructMatch(*"COCC"_smiles, xqm).size() == 1);
      CHECK(SubstructMatch(*"COOCC"_smiles, xqm).empty());
      CHECK(fingerprintsMatch(*"COCC"_smiles, xqm));
      CHECK(!fingerprintsMatch(*"COOCC"_smiles, xqm));
    }
  }
  SECTION("MolBundle") {
    auto mol = "COCC |LN:1:1.3|"_smiles;
    REQUIRE(mol);
    auto txqm = createExtendedQueryMol(*mol);
    ExtendedQueryMol xqm1(txqm);
    ExtendedQueryMol xqm2(std::make_unique<RWMol>(*mol));
    xqm2 = txqm;

    for (const auto &xqm : {txqm, xqm1, xqm2}) {
      CHECK(std::holds_alternative<ExtendedQueryMol::MolBundle_T>(xqm.xqmol));
      CHECK(SubstructMatch(*"COCC"_smiles, xqm).size() == 1);
      CHECK(SubstructMatch(*"COOCC"_smiles, xqm).size() == 1);
      CHECK(SubstructMatch(*"COOOCC"_smiles, xqm).size() == 1);
      CHECK(SubstructMatch(*"COOOOCC"_smiles, xqm).empty());
      CHECK(fingerprintsMatch(*"COCC"_smiles, xqm));
      CHECK(fingerprintsMatch(*"COOCC"_smiles, xqm));
      CHECK(fingerprintsMatch(*"COOOCC"_smiles, xqm));
    }
  }
  SECTION("TautomerQuery") {
    auto mol1 = "CC1OC(N)=N1"_smiles;
    REQUIRE(mol1);
    auto txqm = createExtendedQueryMol(*mol1);
    ExtendedQueryMol xqm1(txqm);
    ExtendedQueryMol xqm2(std::make_unique<RWMol>(*mol1));
    xqm2 = txqm;

    for (const auto &xqm : {txqm, xqm1, xqm2}) {
      CHECK(
          std::holds_alternative<ExtendedQueryMol::TautomerQuery_T>(xqm.xqmol));
      CHECK(SubstructMatch(*"CCC1OC(N)=N1"_smiles, xqm).size() == 1);
      CHECK(SubstructMatch(*"CCC1OC(=N)N1"_smiles, *mol1).empty());
      CHECK(SubstructMatch(*"CCC1OC(=N)N1"_smiles, xqm).size() == 1);
      CHECK(SubstructMatch(*"c1[nH]ncc1"_smiles, xqm).empty());
      CHECK(fingerprintsMatch(*"CCC1OC(N)=N1"_smiles, xqm));
      CHECK(fingerprintsMatch(*"CCC1OC(=N)N1"_smiles, xqm));
      CHECK(!fingerprintsMatch(*"c1[nH]ncc1"_smiles, xqm));
    }
  }
  SECTION("TautomerBundle") {
    auto mol1 = "COCC1OC(N)=N1 |LN:1:1.3|"_smiles;
    REQUIRE(mol1);
    auto txqm = createExtendedQueryMol(*mol1);
    ExtendedQueryMol xqm1(txqm);
    ExtendedQueryMol xqm2(std::make_unique<RWMol>(*mol1));
    xqm2 = txqm;

    for (const auto &xqm : {txqm, xqm1, xqm2}) {
      CHECK(std::holds_alternative<ExtendedQueryMol::TautomerBundle_T>(
          xqm.xqmol));
      CHECK(SubstructMatch(*"COCC1(F)OC(N)=N1"_smiles, xqm).size() == 1);
      CHECK(SubstructMatch(*"COOCC1(F)OC(=N)N1"_smiles, xqm).size() == 1);
      CHECK(SubstructMatch(*"COCC1OC(N)=N1"_smiles, xqm).size() == 1);
      CHECK(SubstructMatch(*"COOOOCC1OC(=N)N1"_smiles, xqm).empty());
      CHECK(fingerprintsMatch(*"COCC1(F)OC(N)=N1"_smiles, xqm));
      CHECK(fingerprintsMatch(*"COOCC1(F)OC(=N)N1"_smiles, xqm));
      CHECK(fingerprintsMatch(*"COCC1OC(N)=N1"_smiles, xqm));
    }
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
    CHECK(fingerprintsMatch(*"FCN(C)CC"_smiles, xqm));
    CHECK(fingerprintsMatch(*"FCN(O)N(C)CC"_smiles, xqm));
  }
}

// there's some redundancy in testing with what's above, but duplicating tests
// isn't a terrible thing
TEST_CASE("adjustQueryProperties") {
  bool doEnumeration = true;
  bool doTautomers = true;
  bool adjustQueryProperties = true;
  SECTION("RWMol") {
    auto mol = "COC1CC1"_smiles;
    REQUIRE(mol);
    auto xqm1 = createExtendedQueryMol(*mol, doEnumeration, doTautomers);
    auto xqm2 = createExtendedQueryMol(*mol, doEnumeration, doTautomers,
                                       adjustQueryProperties);
    MolOps::AdjustQueryParameters aqps =
        MolOps::AdjustQueryParameters::noAdjustments();
    aqps.makeAtomsGeneric = true;
    auto xqm3 = createExtendedQueryMol(*mol, doEnumeration, doTautomers,
                                       adjustQueryProperties, aqps);
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
    CHECK(fingerprintsMatch(*"COC1CC1"_smiles, xqm1));
    CHECK(fingerprintsMatch(*"COC1CC1"_smiles, xqm2));
    CHECK(fingerprintsMatch(*"COC1CC1"_smiles, xqm3));
    CHECK(fingerprintsMatch(*"COC1C(C)C1"_smiles, xqm1));
    CHECK(fingerprintsMatch(*"COC1C(C)C1"_smiles, xqm3));
    CHECK(!fingerprintsMatch(*"COC1OC1"_smiles, xqm1));
    CHECK(!fingerprintsMatch(*"COC1OC1"_smiles, xqm2));
    CHECK(fingerprintsMatch(*"COC1OC1"_smiles, xqm3));
  }
  SECTION("MolBundle") {
    auto mol = "COCC |LN:1:1.3|"_smiles;
    REQUIRE(mol);
    auto xqm1 = createExtendedQueryMol(*mol, doEnumeration, doTautomers);
    MolOps::AdjustQueryParameters aqps =
        MolOps::AdjustQueryParameters::noAdjustments();
    aqps.makeBondsGeneric = true;
    auto xqm2 = createExtendedQueryMol(*mol, doEnumeration, doTautomers,
                                       adjustQueryProperties, aqps);
    CHECK(std::holds_alternative<ExtendedQueryMol::MolBundle_T>(xqm1.xqmol));
    CHECK(std::holds_alternative<ExtendedQueryMol::MolBundle_T>(xqm2.xqmol));
    CHECK(SubstructMatch(*"COC=C"_smiles, xqm1).empty());
    CHECK(SubstructMatch(*"COOC=C"_smiles, xqm1).empty());
    CHECK(SubstructMatch(*"COC=C"_smiles, xqm2).size() == 1);
    CHECK(SubstructMatch(*"COOC=C"_smiles, xqm2).size() == 1);
    CHECK(!fingerprintsMatch(*"COC=C"_smiles, xqm1));
    CHECK(!fingerprintsMatch(*"COOC=C"_smiles, xqm1));
    CHECK(fingerprintsMatch(*"COC=C"_smiles, xqm2));
    CHECK(fingerprintsMatch(*"COOC=C"_smiles, xqm2));
  }
  SECTION("TautomerQuery") {
    auto mol1 = "CC1OC(N)=N1"_smiles;
    REQUIRE(mol1);
    auto xqm1 = createExtendedQueryMol(*mol1, doEnumeration, doTautomers);
    auto xqm2 = createExtendedQueryMol(*mol1, doEnumeration, doTautomers,
                                       adjustQueryProperties);

    CHECK(
        std::holds_alternative<ExtendedQueryMol::TautomerQuery_T>(xqm1.xqmol));
    CHECK(
        std::holds_alternative<ExtendedQueryMol::TautomerQuery_T>(xqm2.xqmol));
    CHECK(SubstructMatch(*"CC1OC(N)=N1"_smiles, xqm1).size() == 1);
    CHECK(SubstructMatch(*"CC1OC(N)=N1"_smiles, xqm2).size() == 1);
    CHECK(SubstructMatch(*"CC1(F)OC(N)=N1"_smiles, xqm1).size() == 1);
    CHECK(SubstructMatch(*"CC1(F)OC(=N)N1"_smiles, xqm1).size() == 1);
    CHECK(SubstructMatch(*"CC1(F)OC(N)=N1"_smiles, xqm2).empty());
    CHECK(SubstructMatch(*"CC1(F)OC(=N)N1"_smiles, xqm2).empty());
    CHECK(fingerprintsMatch(*"CC1OC(N)=N1"_smiles, xqm1));
    CHECK(fingerprintsMatch(*"CC1OC(N)=N1"_smiles, xqm2));
    CHECK(fingerprintsMatch(*"CC1(F)OC(N)=N1"_smiles, xqm1));
    CHECK(fingerprintsMatch(*"CC1(F)OC(=N)N1"_smiles, xqm1));
  }
  SECTION("TautomerBundle") {
    auto mol1 = "COCC1OC(N)=N1 |LN:1:1.3|"_smiles;
    REQUIRE(mol1);
    auto xqm1 = createExtendedQueryMol(*mol1, doEnumeration, doTautomers);
    auto xqm2 = createExtendedQueryMol(*mol1, doEnumeration, doTautomers,
                                       adjustQueryProperties);
    CHECK(
        std::holds_alternative<ExtendedQueryMol::TautomerBundle_T>(xqm1.xqmol));
    CHECK(SubstructMatch(*"COCC1OC(N)=N1"_smiles, xqm1).size() == 1);
    CHECK(SubstructMatch(*"COOCC1OC(=N)N1"_smiles, xqm1).size() == 1);
    CHECK(SubstructMatch(*"COCC1OC(N)=N1"_smiles, xqm2).size() == 1);
    CHECK(SubstructMatch(*"COOCC1OC(=N)N1"_smiles, xqm2).size() == 1);
    CHECK(SubstructMatch(*"COCC1(F)OC(N)=N1"_smiles, xqm1).size() == 1);
    CHECK(SubstructMatch(*"COOCC1(F)OC(=N)N1"_smiles, xqm1).size() == 1);
    CHECK(SubstructMatch(*"COCC1(F)OC(N)=N1"_smiles, xqm2).empty());
    CHECK(SubstructMatch(*"COOCC1(F)OC(=N)N1"_smiles, xqm2).empty());
    CHECK(fingerprintsMatch(*"COCC1OC(N)=N1"_smiles, xqm1));
    CHECK(fingerprintsMatch(*"COOCC1OC(=N)N1"_smiles, xqm1));
    CHECK(fingerprintsMatch(*"COCC1OC(N)=N1"_smiles, xqm2));
    CHECK(fingerprintsMatch(*"COOCC1OC(=N)N1"_smiles, xqm2));
    CHECK(fingerprintsMatch(*"COCC1(F)OC(N)=N1"_smiles, xqm1));
    CHECK(fingerprintsMatch(*"COOCC1(F)OC(=N)N1"_smiles, xqm1));
  }
}

TEST_CASE("controlling which steps are applied") {
  auto mol1 = "COCC1OC(N)=N1 |LN:1:1.3|"_smiles;
  bool doEnumeration = true;
  bool doTautomers = true;
  REQUIRE(mol1);
  {
    auto xqm = createExtendedQueryMol(*mol1, doEnumeration, doTautomers);
    CHECK(
        std::holds_alternative<ExtendedQueryMol::TautomerBundle_T>(xqm.xqmol));
  }
  {
    auto xqm = createExtendedQueryMol(*mol1, doEnumeration, false);
    CHECK(std::holds_alternative<ExtendedQueryMol::MolBundle_T>(xqm.xqmol));
  }
  {
    auto xqm = createExtendedQueryMol(*mol1, false, doTautomers);
    CHECK(std::holds_alternative<ExtendedQueryMol::TautomerQuery_T>(xqm.xqmol));
  }
  {
    auto xqm = createExtendedQueryMol(*mol1, false, false);
    CHECK(std::holds_alternative<ExtendedQueryMol::RWMol_T>(xqm.xqmol));
  }
}

TEST_CASE("interaction with generic groups") {
  // setup a molecule with a generic query:
  auto baseQ = "COC1=NNC(*)=C1 |$;;;;;;AEL_p;$,LN:1:1.3|"_smiles;
  REQUIRE(baseQ);
  auto aqps = MolOps::AdjustQueryParameters::noAdjustments();
  aqps.makeDummiesQueries = true;
  MolOps::adjustQueryProperties(*baseQ, &aqps);
  GenericGroups::setGenericQueriesFromProperties(*baseQ);
  CHECK(baseQ->getAtomWithIdx(6)->hasProp(
      common_properties::_QueryAtomGenericLabel));

  auto xqm = createExtendedQueryMol(*baseQ);
  CHECK(std::holds_alternative<ExtendedQueryMol::TautomerBundle_T>(xqm.xqmol));
  const auto &tquery = std::get<ExtendedQueryMol::TautomerBundle_T>(xqm.xqmol);
  CHECK(tquery->at(0)->getTautomers()[0]->getAtomWithIdx(7)->hasProp(
      common_properties::_QueryAtomGenericLabel));
  CHECK(tquery->at(0)->getTautomers()[1]->getAtomWithIdx(7)->hasProp(
      common_properties::_QueryAtomGenericLabel));
  SECTION("serialization") {
    ExtendedQueryMol xqm2(xqm.toBinary());
    CHECK(
        std::holds_alternative<ExtendedQueryMol::TautomerBundle_T>(xqm2.xqmol));
    const auto &tquery =
        std::get<ExtendedQueryMol::TautomerBundle_T>(xqm2.xqmol);
    CHECK(tquery->at(0)->getTautomers()[0]->getAtomWithIdx(7)->hasProp(
        common_properties::_QueryAtomGenericLabel));
    CHECK(tquery->at(0)->getTautomers()[1]->getAtomWithIdx(7)->hasProp(
        common_properties::_QueryAtomGenericLabel));
  }
}