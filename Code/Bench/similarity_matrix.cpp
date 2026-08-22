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
#include <iostream>
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

// A "small" matrix: below the kernel's multi-threading work threshold, so
// tanimotoMatrix runs single-threaded here. The small-matrix benchmarks
// therefore isolate the SIMD / loop-structure win from the threading win.
constexpr std::size_t kProbeCountSmall = 64;
constexpr std::size_t kTargetCountSmall = 1024;

// A "large" matrix: above the threshold, so tanimotoMatrix splits the work
// across several worker threads. Comparing the large naive loop against the
// large bulk call shows the combined SIMD + threading win.
constexpr std::size_t kProbeCountLarge = 256;
constexpr std::size_t kTargetCountLarge = 4096;

// Runs a pre-packed tanimotoMatrix benchmark for the given dimensions.
void benchBulk(const char *name, std::size_t numProbes,
               std::size_t numTargets) {
  auto probes = buildPopulation(numProbes);
  auto targets = buildPopulation(numTargets);
  auto probeViews = viewsOf(probes);
  auto targetViews = viewsOf(targets);

  // Pre-pack so the benchmark times the inner loop, not the packing.
  std::size_t probeBits = 0;
  std::size_t targetBits = 0;
  auto packedProbes = BulkSimilarity::packFingerprints(probeViews, probeBits);
  auto packedTargets = BulkSimilarity::packFingerprints(targetViews, targetBits);
  const std::size_t words = BulkSimilarity::wordsForBits(probeBits);
  std::vector<double> out(probes.size() * targets.size(), 0.0);

  BENCHMARK(name) {
    BulkSimilarity::tanimotoMatrix(packedProbes.data(), probes.size(),
                                   packedTargets.data(), targets.size(), words,
                                   out.data());
    return out[0];
  };
}

// Runs the naive nested-loop reference for the given dimensions.
void benchNaive(const char *name, std::size_t numProbes,
                std::size_t numTargets) {
  auto probes = buildPopulation(numProbes);
  auto targets = buildPopulation(numTargets);

  BENCHMARK(name) {
    double sum = 0.0;
    for (const auto &p : probes) {
      for (const auto &t : targets) {
        sum += TanimotoSimilarity(*p, *t);
      }
    }
    return sum;
  };
}

}  // namespace

TEST_CASE("Similarity matrix: small (single-threaded, isolates SIMD)",
          "[similarity-matrix]") {
  // Report the popcount kernel the runtime dispatcher actually selected
  // (authoritative; CPUID-based, unlike grepping /proc/cpuinfo which some
  // virtualized hosts mask).
  std::cout << "[bulksim] active popcount kernel = "
            << BulkSimilarity::activeKernel() << std::endl;

  benchNaive("naive nested TanimotoSimilarity (64 x 1024)", kProbeCountSmall,
             kTargetCountSmall);
  benchBulk("BulkSimilarity::tanimotoMatrix (64 x 1024)", kProbeCountSmall,
            kTargetCountSmall);
}

TEST_CASE("Similarity matrix: large (multi-threaded, SIMD + threads)",
          "[similarity-matrix]") {
  benchNaive("naive nested TanimotoSimilarity (256 x 4096)", kProbeCountLarge,
             kTargetCountLarge);
  benchBulk("BulkSimilarity::tanimotoMatrix (256 x 4096)", kProbeCountLarge,
            kTargetCountLarge);
}
