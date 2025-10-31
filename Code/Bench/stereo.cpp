#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

TEST_CASE("Chirality::findPotentialStereo", "[stereo]") {
  for (auto smiles : bench_common::CASES) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);

    BENCHMARK("Chirality::findPotentialStereo: " + std::string(smiles)) {
      return Chirality::findPotentialStereo(*mol);
    };
  }
}

TEST_CASE("CIPLabeler::CIPLabeler", "[stereo]") {
  for (auto smiles : bench_common::CASES) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);

    BENCHMARK("CIPLabeler::assignCIPLabels: " + std::string(smiles)) {
      return Chirality::findPotentialStereo(*mol);
    };
  }
}
