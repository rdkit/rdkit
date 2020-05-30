
#define CATCH_CONFIG_MAIN

#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/TautomerQuery/TautomerQuery.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <DataStructs/BitOps.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/QueryOps.h>

// #define VERBOSE 1

using namespace RDKit;

TEST_CASE("TEST_ENOL") {
  auto mol = "O=C1CCCCC1"_smiles;

  REQUIRE(mol);
  auto tautomerQuery = TautomerQuery::fromMol(*mol);
  auto tautomers = tautomerQuery->getTautomers();
  CHECK(tautomers.size() == 2);

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
  auto nMatches = SubstructMatch(*target1, *tautomerQuery, matchVect);
  CHECK(nMatches == 1);

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
  auto tautomerQuery = TautomerQuery::fromMol(*mol);
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
  auto tautomerQuery = TautomerQuery::fromMol(*mol);
  auto templateMol = tautomerQuery->getTemplateMolecule();

  std::string pickle;
  MolPickler::pickleMol(*templateMol, pickle);
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
  auto tautomerQuery = TautomerQuery::fromMol(*mol);
  auto templateMol = tautomerQuery->getTemplateMolecule();

  // this test molecule has complex query bonds where the template has query
  // bonds, but they are not tautomer bonds.
  RWMol molWithoutTautomerBonds(*mol);
  std::vector<std::pair<int, int>> atomIndexes;
  for (auto modifiedBondIdx : tautomerQuery->getModifiedBonds()) {
    auto bondQuery = makeSingleOrAromaticBondQuery();
    auto queryBond = new QueryBond();
    queryBond->setQuery(bondQuery);
    molWithoutTautomerBonds.replaceBond(modifiedBondIdx, queryBond, true);
    delete queryBond;
  }

  // The molecule without tautomer bonds has the same regular fingerprint as the
  // template
#ifdef VERBOSE
  std::cout << std::endl << "fingerprinting template" << std::endl;
#endif
  auto templateQueryFingerprint = PatternFingerprintMol(*templateMol);
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
    auto tautomerQuery = TautomerQuery::fromMol(*mol);
    CHECK(1 == tautomerQuery->getTautomers().size());
    CHECK(0 == tautomerQuery->getModifiedAtoms().size());
    CHECK(0 == tautomerQuery->getModifiedBonds().size());
    auto target = "CC1=NC2=CC=CC=C2O1"_smiles;
    REQUIRE(target);
    CHECK(tautomerQuery->isSubstructOf(*target));
}

