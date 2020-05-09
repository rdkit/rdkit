
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/TautomerQuery/TautomerQuery.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

TEST_CASE("TEST_ENOL") {
    auto mol = "O=C1CCCCC1"_smiles;
    
    REQUIRE(mol);
    auto tautomerQuery = TautomerQuery::fromMol(*mol);
    auto target = "OC1=CCCC(CC)C1"_smiles;
    SubstructMatchParameters params;
    std::vector<ROMOL_SPTR> matchingTautomers;
    auto matches = tautomerQuery->SubstructMatch(*target, params, matchingTautomers);

    CHECK(matches.size() ==1);
    REQUIRE(matchingTautomers.size() == 1);
    auto tautomerSmiles = MolToSmiles(*matchingTautomers[0]);
    CHECK(tautomerSmiles == "OC1=CCCCC1");
}
