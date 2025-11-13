#include <catch2/catch_all.hpp>
#include <string>
#include <sstream>

#include "bench_common.hpp"

#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

TEST_CASE("MolPickler::pickleMol", "[pickle]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("MolPickler::pickleMol") {
    std::stringstream buf;
    for (auto &mol : samples) {
      MolPickler::pickleMol(mol, buf);
    }
    return buf.str().size();
  };
}

TEST_CASE("MolPickler::molFromPickle", "[pickle]") {
  auto samples = bench_common::load_samples();
  std::vector<std::string> pickles;
  pickles.reserve(samples.size());
  for (auto &mol : samples) {
    std::string pickled;
    MolPickler::pickleMol(mol, pickled);
    pickles.push_back(std::move(pickled));
  }
  BENCHMARK("MolPickler::molFromPickle") {
    auto total_atoms = 0;
    for (auto &pickled : pickles) {
      ROMol res;
      MolPickler::molFromPickle(pickled, res);
      total_atoms += res.getNumAtoms();
    }
    REQUIRE(total_atoms > 0);
    return total_atoms;
  };
}
