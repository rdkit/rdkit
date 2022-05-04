//
// Copyright 2020-2022 Schrodinger, Inc and other RDKit contributors
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/TautomerQuery/TautomerQuery.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <DataStructs/BitOps.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

// #define VERBOSE 1

using namespace RDKit;

TEST_CASE("TEMPLATE_ERROR") {
  // for this guy the template needs to account for bonds modified when
  // tautomers are sanitized
  auto mol = "Cc1nc2ccccc2[nH]1"_smiles;
  REQUIRE(mol);
  auto target = "CN1C2=C(C(=O)Nc3ccccc3)C(=O)CCN2c2ccccc21"_smiles;
  REQUIRE(target);

  auto tautomerQuery =
      std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*mol));
  auto match = false;
  MatchVectType matchVect;
  for (auto taut : tautomerQuery->getTautomers()) {
    auto test = SubstructMatch(*target, *taut, matchVect);
#ifdef VERBOSE
    std::cout << "Tautomer " << MolToSmiles(*taut) << " match " << test
              << std::endl;
#endif
    if (test) {
      match = true;
    }
  }
  CHECK(match);

  SubstructMatchParameters params;
  std::vector<ROMOL_SPTR> matchingTautomers;
  auto matches =
      tautomerQuery->substructOf(*target, params, &matchingTautomers);
  CHECK(matches.size() == 1);
  REQUIRE(matchingTautomers.size() == 1);
}

TEST_CASE("TEST_UNIQUIFY") {
  auto mol = "O=C1CCCCC1"_smiles;
  REQUIRE(mol);
  auto target = "O=C1CCCC(CC)C1"_smiles;
  REQUIRE(target);

  auto tautomerQuery =
      std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*mol));
  auto tautomers = tautomerQuery->getTautomers();
  SubstructMatchParameters params;
  params.maxMatches = 1000;
  std::vector<ROMOL_SPTR> matchingTautomers;

  params.uniquify = true;
  auto matches =
      tautomerQuery->substructOf(*target, params, &matchingTautomers);
  CHECK(matches.size() == 1);
  REQUIRE(matchingTautomers.size() == 1);

  params.uniquify = false;
  matches = tautomerQuery->substructOf(*target, params, &matchingTautomers);
  CHECK(matches.size() == 2);
  REQUIRE(matchingTautomers.size() == 2);
  CHECK(matchingTautomers[0] == matchingTautomers[1]);
}

TEST_CASE("DIFFERENT_TO_ENUMERATED") {
  // test shows we need to set uniquify = false when matching template
  auto mol = "NC(N)=O"_smiles;
  auto tautomerQuery =
      std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*mol));
  auto tautomers = tautomerQuery->getTautomers();
  // auto target =
  // "NC1=NC2(CO1)c1cc(-c3cccnc3F)ccc1Oc1cnc(C3=CCCOC3)cc12"_smiles;
  auto target = "NC1=NCCO1"_smiles;
  auto enumMatch = false;
  MatchVectType matchVect;
  for (auto t : tautomers) {
    if (SubstructMatch(*target, *t, matchVect)) {
      enumMatch = true;
      break;
    }
  }
  CHECK(enumMatch);

  auto match = tautomerQuery->isSubstructOf(*target);
  CHECK(match);
}

TEST_CASE("SIMPLE_ERROR") {
  auto mol = "CC=O"_smiles;
  REQUIRE(mol);
  auto tautomerQuery =
      std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*mol));
  auto target = "OC(C)=O"_smiles;
  auto matches = SubstructMatch(*target, *mol);
  CHECK(matches.size() == 1);
  auto match = tautomerQuery->isSubstructOf(*target);
  CHECK(match);
  target = "CCN(CC)C(=O)COP(=O)(O)COCCn1cnc2c(N)ncnc21"_smiles;
  match = tautomerQuery->isSubstructOf(*target);
  CHECK(match);
}

TEST_CASE("TEST_ENOL") {
  auto mol = "O=C1CCCCC1"_smiles;

  REQUIRE(mol);
  auto tautomerQuery =
      std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*mol));
  auto tautomers = tautomerQuery->getTautomers();
  CHECK(tautomers.size() == 2);
  auto modifiedAtoms = tautomerQuery->getModifiedAtoms();
  CHECK(modifiedAtoms.size() == 3);
  auto modifiedBonds = tautomerQuery->getModifiedBonds();
  CHECK(modifiedBonds.size() == 3);

  auto target1 = "OC1=CCCC(CC)C1"_smiles;
  REQUIRE(target1);
  SubstructMatchParameters params;
  std::vector<ROMOL_SPTR> matchingTautomers;
  auto matches =
      tautomerQuery->substructOf(*target1, params, &matchingTautomers);
  CHECK(matches.size() == 1);
  auto match = tautomerQuery->isSubstructOf(*target1, params);
  CHECK(match);

  REQUIRE(matchingTautomers.size() == 1);
  auto tautomerSmiles = MolToSmiles(*matchingTautomers[0]);
  CHECK(tautomerSmiles == "OC1=CCCCC1");

  auto target2 = "O=C1CCCC(CC)C1"_smiles;
  REQUIRE(target2);
  matches = tautomerQuery->substructOf(*target2, params, &matchingTautomers);

  CHECK(matches.size() == 1);
  REQUIRE(matchingTautomers.size() == 1);
  tautomerSmiles = MolToSmiles(*matchingTautomers[0]);
  CHECK(tautomerSmiles == "O=C1CCCCC1");

  MatchVectType matchVect;
  auto hasMatch = SubstructMatch(*target1, *tautomerQuery, matchVect);
  CHECK(hasMatch);

  auto templateFingerpint = tautomerQuery->patternFingerprintTemplate();
  REQUIRE(templateFingerpint);

  auto target1Fingerprint = TautomerQuery::patternFingerprintTarget(*target1);
  REQUIRE(target1Fingerprint);
  CHECK(AllProbeBitsMatch(*templateFingerpint, *target1Fingerprint));
  delete target1Fingerprint;

  auto target2Fingerprint = TautomerQuery::patternFingerprintTarget(*target2);
  REQUIRE(target2Fingerprint);
  CHECK(AllProbeBitsMatch(*templateFingerpint, *target2Fingerprint));
  delete target2Fingerprint;
  delete templateFingerpint;
}

TEST_CASE("TEST_COMPLEX") {
  auto mol = "Nc1nc(=O)c2nc[nH]c2[nH]1"_smiles;
  REQUIRE(mol);
  auto tautomerQuery =
      std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*mol));
  CHECK(15 == tautomerQuery->getTautomers().size());

  auto queryFingerprint = tautomerQuery->patternFingerprintTemplate();
  REQUIRE(queryFingerprint);
  std::vector<std::string> targetSmis{"CCc1nc2[nH]c(=N)nc(O)c2[nH]1",
                                      "CN1C2=NC=NC2=C(O)N=C1N"};
  for (auto targetSmiles : targetSmis) {
    auto target = SmilesToMol(targetSmiles);
    REQUIRE(target);
    CHECK(tautomerQuery->isSubstructOf(*target));
    auto targetFingerprint = TautomerQuery::patternFingerprintTarget(*target);
    REQUIRE(targetFingerprint);
    CHECK(AllProbeBitsMatch(*queryFingerprint, *targetFingerprint));
    delete targetFingerprint;
    delete target;
  }
  delete queryFingerprint;
}

TEST_CASE("TEST_PICKLE") {
  auto mol = "O=C1CCCCC1"_smiles;
  REQUIRE(mol);
  auto tautomerQuery =
      std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*mol));
  auto templateMol = tautomerQuery->getTemplateMolecule();

  std::string pickle;
  MolPickler::pickleMol(templateMol, pickle);
  ROMol pickleMol;
  MolPickler::molFromPickle(pickle, pickleMol);

  for (auto modifiedBondIdx : tautomerQuery->getModifiedBonds()) {
    auto modifiedBond = pickleMol.getBondWithIdx(modifiedBondIdx);
    REQUIRE(modifiedBond->hasQuery());
    CHECK(modifiedBond->getQuery()->getDescription() ==
          "SingleOrDoubleOrAromaticBond");
  }
}

TEST_CASE("TEST_FINGERPRINT") {
  auto mol = "O=C1CCCCC1"_smiles;
  REQUIRE(mol);
  auto tautomerQuery =
      std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*mol));
  auto templateMol = tautomerQuery->getTemplateMolecule();

  // this test molecule has complex query bonds where the template has query
  // bonds, but they are not identified as tautomer bonds.
  RWMol molWithoutTautomerBonds(*mol);
  std::vector<std::pair<int, int>> atomIndexes;
  for (auto modifiedBondIdx : tautomerQuery->getModifiedBonds()) {
    auto queryBond = new QueryBond();
    queryBond->setQuery(makeBondOrderEqualsQuery(Bond::BondType::SINGLE));
    queryBond->expandQuery(makeBondOrderEqualsQuery(Bond::BondType::AROMATIC),
                           Queries::COMPOSITE_OR);
    molWithoutTautomerBonds.replaceBond(modifiedBondIdx, queryBond, true);
    delete queryBond;
  }

  // The molecule without tautomer bonds has the same regular fingerprint as the
  // template
#ifdef VERBOSE
  std::cout << std::endl << "fingerprinting template" << std::endl;
#endif
  auto templateQueryFingerprint = PatternFingerprintMol(templateMol);
#ifdef VERBOSE
  std::cout << "fingerprinting mol without bonds" << std::endl;
#endif
  auto molWithoutTautomerBondsFingerprint =
      PatternFingerprintMol(molWithoutTautomerBonds);

  REQUIRE(templateQueryFingerprint);
  REQUIRE(molWithoutTautomerBondsFingerprint);
  CHECK(AllProbeBitsMatch(*templateQueryFingerprint,
                          *molWithoutTautomerBondsFingerprint));
  CHECK(AllProbeBitsMatch(*molWithoutTautomerBondsFingerprint,
                          *templateQueryFingerprint));
  delete templateQueryFingerprint;
  delete molWithoutTautomerBondsFingerprint;

  // The tautomer fingerprint for the molecule without tautomer bonds has a
  // subset of the template's tautomeric fingerprint.
  templateQueryFingerprint = tautomerQuery->patternFingerprintTemplate();
  molWithoutTautomerBondsFingerprint =
      TautomerQuery::patternFingerprintTarget(molWithoutTautomerBonds);
  REQUIRE(templateQueryFingerprint);
  REQUIRE(molWithoutTautomerBondsFingerprint);
#ifdef VERBOSE
  std::cout << std::endl << "molWithoutTautomerBonds" << std::endl;
  std::vector<int> onBits;
  molWithoutTautomerBondsFingerprint->getOnBits(onBits);
  for (auto bit : onBits) {
    std::cout << bit << " ";
  }
  std::cout << std::endl;
  std::cout << std::endl << "template" << std::endl;
  templateQueryFingerprint->getOnBits(onBits);
  for (auto bit : onBits) {
    std::cout << bit << " ";
  }
  std::cout << std::endl;
#endif
  CHECK(AllProbeBitsMatch(*molWithoutTautomerBondsFingerprint,
                          *templateQueryFingerprint));
  CHECK(!AllProbeBitsMatch(*templateQueryFingerprint,
                           *molWithoutTautomerBondsFingerprint));

  // The tautomer fingerprint for the template is a subset of the tautomer
  // fingerprint for the original molecule.
  auto molFingerprint = TautomerQuery::patternFingerprintTarget(*mol);
  CHECK(AllProbeBitsMatch(*templateQueryFingerprint, *molFingerprint));
  // This expected bit count for the template tautomeric fingerprint applies for
  // all queries if there are no bit clashes
  auto expectedTemplateBitsCount =
      molWithoutTautomerBondsFingerprint->getNumOnBits() +
      (molFingerprint->getNumOnBits() -
       molWithoutTautomerBondsFingerprint->getNumOnBits()) /
          2;
  CHECK(expectedTemplateBitsCount == templateQueryFingerprint->getNumOnBits());

  delete templateQueryFingerprint;
  delete molWithoutTautomerBondsFingerprint;
  delete molFingerprint;
}

TEST_CASE("TEST_NOT_TAUTOMER") {
  auto mol = "c1ccccc1"_smiles;
  REQUIRE(mol);
  auto tautomerQuery =
      std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*mol));
  CHECK(1 == tautomerQuery->getTautomers().size());
  CHECK(0 == tautomerQuery->getModifiedAtoms().size());
  CHECK(0 == tautomerQuery->getModifiedBonds().size());
  auto target = "CC1=NC2=CC=CC=C2O1"_smiles;
  REQUIRE(target);
  CHECK(tautomerQuery->isSubstructOf(*target));
}

TEST_CASE("github #3821 TAUTOMERQUERY_COPY_CONSTRUCTOR") {
  auto mol = "c1ccccc1"_smiles;
  auto tautomerQuery =
      std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*mol));
  auto tautomerQueryCopyConstructed =
      std::unique_ptr<TautomerQuery>(new TautomerQuery(*tautomerQuery));
  CHECK(&(tautomerQuery->getTemplateMolecule()) !=
        &tautomerQueryCopyConstructed->getTemplateMolecule());
}

TEST_CASE("github #3821 check TAUTOMERQUERY_OPERATOR= does a deep copy") {
  auto mol = "c1ccccc1"_smiles;
  auto tautomerQuery =
      std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*mol));
  auto tautomerQueryAssigned = *tautomerQuery;
  CHECK(&(tautomerQuery->getTemplateMolecule()) !=
        &tautomerQueryAssigned.getTemplateMolecule());
}

TEST_CASE("Serialization") {
#ifdef RDK_USE_BOOST_SERIALIZATION
  SECTION("basics") {
    auto mol = "Nc1nc(=O)c2nc[nH]c2[nH]1"_smiles;
    REQUIRE(mol);
    auto tautomerQuery =
        std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*mol));
    CHECK(15 == tautomerQuery->getTautomers().size());

    std::string pickle = tautomerQuery->serialize();
    TautomerQuery serialized(pickle);
    CHECK(serialized.getTautomers().size() ==
          tautomerQuery->getTautomers().size());

    auto queryFingerprint = serialized.patternFingerprintTemplate();
    REQUIRE(queryFingerprint);
    std::vector<std::string> targetSmis{"CCc1nc2[nH]c(=N)nc(O)c2[nH]1",
                                        "CN1C2=NC=NC2=C(O)N=C1N"};
    for (auto targetSmiles : targetSmis) {
      auto target = SmilesToMol(targetSmiles);
      REQUIRE(target);
      CHECK(serialized.isSubstructOf(*target));
      auto targetFingerprint = TautomerQuery::patternFingerprintTarget(*target);
      REQUIRE(targetFingerprint);
      CHECK(AllProbeBitsMatch(*queryFingerprint, *targetFingerprint));
      delete targetFingerprint;
      delete target;
    }
    delete queryFingerprint;
  }
#endif
}
