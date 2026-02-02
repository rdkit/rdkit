#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>

using namespace RDKit;

TEST_CASE("MorganFingerprints::getFingerprint", "[fingerprint]") {
  auto samples = bench_common::load_samples();
  const auto radius = 2;
  std::unique_ptr<FingerprintGenerator<uint64_t>> gen(
      MorganFingerprint::getMorganGenerator<uint64_t>(radius));

  BENCHMARK("MorganFingerprints::getFingerprint") {
    auto sum = 0;
    for (auto &mol : samples) {
      std::unique_ptr<ExplicitBitVect> fp(gen->getFingerprint(mol));
      sum += fp->getNumOnBits();
    }
    return sum;
  };
}
