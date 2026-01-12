#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <INCHI-API/inchi.h>

using namespace RDKit;

TEST_CASE("MolToInchi", "[inchi]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("MolToInchi") {
    std::vector<std::string> inchis;
    for (auto &mol : samples) {
      ExtraInchiReturnValues rv;
      inchis.push_back(MolToInchi(mol, rv));
    }
    return inchis;
  };
}

TEST_CASE("InchiToInchiKey", "[inchi]") {
  auto samples = bench_common::load_samples();
  std::vector<std::string> inchis;
  for (auto &mol : samples) {
    ExtraInchiReturnValues rv;
    inchis.push_back(MolToInchi(mol, rv));
  }
  BENCHMARK("InchiToInchiKey") {
    std::vector<std::string> inchikeys;
    for (auto &inchi : inchis) {
      inchikeys.push_back(InchiToInchiKey(inchi));
    }
    return inchikeys;
  };
}

TEST_CASE("InchiToMol", "[inchi]") {
  auto samples = bench_common::load_samples();
  std::vector<std::string> inchis;
  for (auto &mol : samples) {
    ExtraInchiReturnValues rv;
    inchis.push_back(MolToInchi(mol, rv));
  }
  BENCHMARK("InchiToMol") {
    std::vector<std::unique_ptr<ROMol>> mols;
    for (auto &inchi : inchis) {
      ExtraInchiReturnValues rv_inner;
      mols.emplace_back(InchiToMol(inchi, rv_inner));
    }
    return mols;
  };
}
