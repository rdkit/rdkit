
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
      tautomerQuery->SubstructMatch(*target1, params, &matchingTautomers);

  CHECK(matches.size() == 1);
  REQUIRE(matchingTautomers.size() == 1);
  auto tautomerSmiles = MolToSmiles(*matchingTautomers[0]);
  CHECK(tautomerSmiles == "OC1=CCCCC1");

  auto target2 = "O=C1CCCC(CC)C1"_smiles;
  REQUIRE(target2);
  matches = tautomerQuery->SubstructMatch(*target2, params, &matchingTautomers);

  CHECK(matches.size() == 1);
  REQUIRE(matchingTautomers.size() == 1);
  tautomerSmiles = MolToSmiles(*matchingTautomers[0]);
  CHECK(tautomerSmiles == "O=C1CCCCC1");

  MatchVectType matchVect;
  auto nMatches = SubstructMatch(*target1, *tautomerQuery, matchVect);
  CHECK(nMatches == 1);

  // I know the fingerprinter is setting bits using the tautomer query
  // bonds (from verbose output), but need a test to prove
  auto templateFingerpint = tautomerQuery->patternFingerprintTemplate();
  REQUIRE(templateFingerpint);

  auto target1Fingerprint =
      PatternFingerprintMol(*target1, 2048U, nullptr, nullptr, true);
  REQUIRE(target1Fingerprint);
  CHECK(AllProbeBitsMatch(*templateFingerpint, *target1Fingerprint));
  delete target1Fingerprint;

  auto target2Fingerprint =
      PatternFingerprintMol(*target2, 2048U, nullptr, nullptr, true);
  REQUIRE(target2Fingerprint);
  CHECK(AllProbeBitsMatch(*templateFingerpint, *target2Fingerprint));
  delete target2Fingerprint;

  delete templateFingerpint;
}
