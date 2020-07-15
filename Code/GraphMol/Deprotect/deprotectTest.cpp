#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Deprotect/Deprotect.h>

using namespace RDKit;

TEST_CASE("Standard deprotections", "[deprotect]") {
  SECTION("simple deprotections") {
    auto m = "N(C(=O)OC(C)(C)C)Cc1ccccc1NC(=O)OC(C)(C)C"_smiles;
    auto res = deprotect(*m);
    REQUIRE(res);
    CHECK(MolToSmiles(*res) == "NCc1ccccc1N");
    CHECK(res->getProp<int>("DEPROTECTION_COUNT") == 2);
    std::vector<std::string> expected {"Boc", "Boc"};
    CHECK(res->getProp<std::vector<std::string>>("DEPROTECTIONS") ==
	  expected);
  }
}
