#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Descriptors/Lipinski.h>

using namespace RDKit;

TEST_CASE("Descriptors::calcNumSpiroAtoms", "[descriptors]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("Descriptors::calcNumSpiroAtoms") {
    auto sum = 0;
    for (auto &mol : samples) {
      sum += Descriptors::calcNumSpiroAtoms(mol);
    }
    return sum;
  };
}

TEST_CASE("Descriptors::calcNumBridgeheadAtoms", "[descriptors]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("Descriptors::calcNumBridgeheadAtoms") {
    auto sum = 0;
    for (auto &mol : samples) {
      sum += Descriptors::calcNumBridgeheadAtoms(mol);
    }
    return sum;
  };
}
