#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>

using namespace RDKit;

TEST_CASE("MorganFingerprints::getFingerprint", "[fingerprint]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("MorganFingerprints::getFingerprint") {
    auto sum = 0;
    for (auto &mol : samples) {
      std::unique_ptr<SparseIntVect<uint32_t>> fp(
          MorganFingerprints::getFingerprint(mol, 2, nullptr, nullptr, false,
                                             true, true, false));
      sum += fp->getTotalVal();
    }
    return sum;
  };
}
