//
//  Copyright (C) 2026 Clay Moore and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include "bench_common.hpp"

#include <DataStructs/BitOps.h>
#include <DataStructs/BulkSimilarity/BulkSimilarity.h>
#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>

#include <cstdint>
#include <memory>
#include <vector>

using namespace RDKit;

namespace {

// Build a pool of fingerprints by replicating the sample molecules and
// xor-perturbing each copy with a deterministic random mask. This keeps the
// content varied enough that the Tanimoto values are non-trivial without
// depending on a large on-disk dataset.
std::vector<std::unique_ptr<ExplicitBitVect>> buildPopulation(
    std::size_t count) {
  auto samples = bench_common::load_samples();
  std::unique_ptr<FingerprintGenerator<std::uint64_t>> gen(
      MorganFingerprint::getMorganGenerator<std::uint64_t>(/*radius=*/2));

  std::vector<std::unique_ptr<ExplicitBitVect>> base;
  base.reserve(samples.size());
  for (auto &mol : samples) {
    base.emplace_back(gen->getFingerprint(mol));
  }

  std::vector<std::unique_ptr<ExplicitBitVect>> out;
  out.reserve(count);
  const unsigned int nBits = base.front()->getNumBits();
  for (std::size_t i = 0; i < count; ++i) {
    auto fp = std::make_unique<ExplicitBitVect>(*base[i % base.size()]);
    const std::uint64_t r = bench_common::nth_random(i);
    fp->setBit(static_cast<unsigned int>(r % nBits));
    fp->setBit(static_cast<unsigned int>((r >> 32) % nBits));
    out.emplace_back(std::move(fp));
  }
  return out;
}

std::vector<const ExplicitBitVect *> viewsOf(
    const std::vector<std::unique_ptr<ExplicitBitVect>> &fps) {
  std::vector<const ExplicitBitVect *> views;
  views.reserve(fps.size());
  for (const auto &p : fps) {
    views.push_back(p.get());
  }
  return views;
}

constexpr std::size_t kProbeCount = 64;
constexpr std::size_t kTargetCount = 1024;

}  // namespace

TEST_CASE("Similarity matrix: naive nested TanimotoSimilarity loop",
          "[similarity-matrix]") {
  auto probes = buildPopulation(kProbeCount);
  auto targets = buildPopulation(kTargetCount);

  BENCHMARK("nested TanimotoSimilarity (64 x 1024)") {
    double sum = 0.0;
    for (const auto &p : probes) {
      for (const auto &t : targets) {
        sum += TanimotoSimilarity(*p, *t);
      }
    }
    return sum;
  };
}

TEST_CASE("Similarity matrix: BulkSimilarity::tanimotoMatrix",
          "[similarity-matrix]") {
  auto probes = buildPopulation(kProbeCount);
  auto targets = buildPopulation(kTargetCount);
  auto probeViews = viewsOf(probes);
  auto targetViews = viewsOf(targets);

  // Pre-pack so the benchmark times the inner loop, not the packing.
  std::size_t probeBits = 0;
  std::size_t targetBits = 0;
  auto packedProbes =
      BulkSimilarity::packFingerprints(probeViews, probeBits);
  auto packedTargets =
      BulkSimilarity::packFingerprints(targetViews, targetBits);
  const std::size_t words = BulkSimilarity::wordsForBits(probeBits);
  std::vector<double> out(probes.size() * targets.size(), 0.0);

  BENCHMARK("BulkSimilarity::tanimotoMatrix (64 x 1024)") {
    BulkSimilarity::tanimotoMatrix(packedProbes.data(), probes.size(),
                                   packedTargets.data(), targets.size(),
                                   words, out.data());
    return out[0];
  };
}
