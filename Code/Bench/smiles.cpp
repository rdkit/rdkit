#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

TEST_CASE("SmilesToMol", "[smiles]") {
  for (auto smiles : bench_common::CASES) {
    BENCHMARK("SmilesToMol: " + std::string(smiles)) {
      std::unique_ptr<ROMol> mol{SmilesToMol(smiles)};
      REQUIRE(mol);
      return mol;
    };
  }
}

TEST_CASE("MolToSmiles", "[smiles]") {
  for (auto smiles : bench_common::CASES) {
    std::unique_ptr<ROMol> mol{SmilesToMol(smiles)};
    REQUIRE(mol);
    BENCHMARK("MolToSmiles: " + std::string(smiles)) {
      return MolToSmiles(*mol);
    };
  }
}
