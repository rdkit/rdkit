#include <catch2/catch_all.hpp>
#include <memory>
#include <string>
#include <vector>

#include "bench_common.hpp"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <DataStructs/BitOps.h>

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

TEST_CASE("Bitmap similarity", "[fingerprint][bitmap]") {
  // Morgan fingerprints laid out as raw bitmaps, the form the FPB reader
  // and the PostgreSQL cartridge compare
  constexpr unsigned int nBits = 2048;
  constexpr unsigned int nBytes = nBits / 8;
  auto samples = bench_common::load_samples();
  std::unique_ptr<FingerprintGenerator<uint64_t>> gen(
      MorganFingerprint::getMorganGenerator<uint64_t>(2));  // 2048 bits
  std::vector<std::vector<unsigned char>> bitmaps;
  for (auto &mol : samples) {
    std::unique_ptr<ExplicitBitVect> fp(gen->getFingerprint(mol));
    std::vector<unsigned char> bytes(nBytes, 0);
    std::vector<int> onBits;
    fp->getOnBits(onBits);
    for (auto bit : onBits) {
      bytes[bit / 8] |= 1 << (bit % 8);
    }
    bitmaps.push_back(std::move(bytes));
  }

  BENCHMARK("CalcBitmapPopcount") {
    unsigned int total = 0;
    for (const auto &bm : bitmaps) {
      total += CalcBitmapPopcount(bm.data(), nBytes);
    }
    return total;
  };

  BENCHMARK("CalcBitmapTanimoto all pairs") {
    double total = 0.0;
    for (size_t i = 0; i < bitmaps.size(); ++i) {
      for (size_t j = 0; j < bitmaps.size(); ++j) {
        total += CalcBitmapTanimoto(bitmaps[i].data(), bitmaps[j].data(),
                                    nBytes);
      }
    }
    return total;
  };
}
