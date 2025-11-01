#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

TEST_CASE("SmilesToMol", "[smiles]") {
  BENCHMARK("SmilesToMol") {
    auto total_atoms = 0;
    for (auto smiles : bench_common::SAMPLES) {
      auto mol = v2::SmilesParse::MolFromSmiles(smiles);
      REQUIRE(mol);
	  total_atoms += mol->getNumAtoms();
    }
    return total_atoms;
  };
}

TEST_CASE("MolToSmiles", "[smiles]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("MolToSmiles") {
    auto total_length = 0;
    for (auto &mol : samples) {
      auto smiles = MolToSmiles(mol);
      total_length += smiles.size();
    }
    return total_length;
  };
}
