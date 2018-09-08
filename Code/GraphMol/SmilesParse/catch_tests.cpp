#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

TEST_CASE("Github #1972", "[SMILES,bug]") {
  SECTION("basics") {
    std::vector<std::vector<std::string>> smiles = {
        {"[C@@]1(Cl)(F)(I).Br1", "[C@@](Br)(Cl)(F)(I)"},
        {"[C@@](Cl)(F)(I)1.Br1", "[C@@](Cl)(F)(I)Br"},
        {"[C@@](Cl)1(F)(I).Br1", "[C@@](Cl)(Br)(F)(I)"},
        {"[C@@](Cl)(F)1(I).Br1", "[C@@](Cl)(F)(Br)(I)"}};
    for (const auto &pr : smiles) {
      std::unique_ptr<ROMol> m1(SmilesToMol(pr[0]));
      std::unique_ptr<ROMol> m2(SmilesToMol(pr[1]));
      REQUIRE(m1);
      REQUIRE(m2);
      auto csmi1 = MolToSmiles(*m1);
      auto csmi2 = MolToSmiles(*m2);
      CHECK(csmi1 == csmi2);
    }
  }
  SECTION("further examples") {
    std::vector<std::vector<std::string>> smiles = {
        {"[C@@]1(Cl)2(I).Br1.F2", "[C@@](Br)(Cl)(F)(I)"},
        {"[C@@](Cl)2(I)1.Br1.F2", "[C@@](Cl)(F)(I)Br"},
        {"[C@@]12(Cl)(I).Br1.F2", "[C@@](Br)(F)(Cl)(I)"},
        {"[C@@]21(Cl)(I).Br1.F2", "[C@@](F)(Br)(Cl)(I)"},
        {"[C@@](Cl)12(I).Br1.F2", "[C@@](Cl)(Br)(F)(I)"},
        {"[C@@](Cl)21(I).Br1.F2", "[C@@](Cl)(F)(Br)(I)"},
        {"[C@@](Cl)(I)21.Br1.F2", "[C@@](Cl)(I)(F)(Br)"},
        {"[C@@](Cl)(I)12.Br1.F2", "[C@@](Cl)(I)(Br)(F)"}};
    for (const auto &pr : smiles) {
      std::unique_ptr<ROMol> m1(SmilesToMol(pr[0]));
      std::unique_ptr<ROMol> m2(SmilesToMol(pr[1]));
      REQUIRE(m1);
      REQUIRE(m2);
      auto csmi1 = MolToSmiles(*m1);
      auto csmi2 = MolToSmiles(*m2);
      CHECK(csmi1 == csmi2);
    }
  }
}

TEST_CASE("Github #2029", "[SMILES,bug]") {
  SECTION("basics") {
    std::unique_ptr<ROMol> m1(SmilesToMol("CN[C@H](Cl)C(=O)O"));
    REQUIRE(m1);
    m1->getBondWithIdx(1)->setBondDir(Bond::BEGINWEDGE);
    bool doKekule = false, allBondsExplicit = false;
    auto bsmi = SmilesWrite::GetBondSmiles(m1->getBondWithIdx(1), -1, doKekule,
                                           allBondsExplicit);
    CHECK(bsmi == "");
    allBondsExplicit = true;
    bsmi = SmilesWrite::GetBondSmiles(m1->getBondWithIdx(1), -1, doKekule,
                                      allBondsExplicit);
    CHECK(bsmi == "-");
  }
}
