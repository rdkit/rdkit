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

#include "BulkSimilarity.h"

#include <DataStructs/BitOps.h>
#include <DataStructs/ExplicitBitVect.h>

#include <memory>
#include <random>
#include <vector>

using RDKit::BulkSimilarity::activeKernel;
using RDKit::BulkSimilarity::packFingerprints;
using RDKit::BulkSimilarity::tanimotoMatrix;
using RDKit::BulkSimilarity::wordsForBits;

namespace {

std::unique_ptr<ExplicitBitVect> makeRandomFp(std::size_t nBits,
                                              std::mt19937 &rng,
                                              double density) {
  auto fp = std::make_unique<ExplicitBitVect>(static_cast<unsigned int>(nBits));
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  for (std::size_t i = 0; i < nBits; ++i) {
    if (dist(rng) < density) {
      fp->setBit(static_cast<unsigned int>(i));
    }
  }
  return fp;
}

std::vector<std::unique_ptr<ExplicitBitVect>> makeRandomFps(
    std::size_t count, std::size_t nBits, std::mt19937 &rng,
    double density = 0.1) {
  std::vector<std::unique_ptr<ExplicitBitVect>> result;
  result.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    result.push_back(makeRandomFp(nBits, rng, density));
  }
  return result;
}

std::vector<const ExplicitBitVect *> asViews(
    const std::vector<std::unique_ptr<ExplicitBitVect>> &fps) {
  std::vector<const ExplicitBitVect *> views;
  views.reserve(fps.size());
  for (const auto &p : fps) {
    views.push_back(p.get());
  }
  return views;
}

}  // namespace

TEST_CASE("BulkSimilarity: wordsForBits validates input",
          "[BulkSimilarity]") {
  CHECK(wordsForBits(64) == 1);
  CHECK(wordsForBits(2048) == 32);
  CHECK_THROWS(wordsForBits(0));
  CHECK_THROWS(wordsForBits(100));
  CHECK_THROWS(wordsForBits(1023));
}

TEST_CASE("BulkSimilarity: packFingerprints round-trips bit values",
          "[BulkSimilarity]") {
  constexpr unsigned int nBits = 256;
  ExplicitBitVect fp(nBits);
  fp.setBit(0);
  fp.setBit(63);
  fp.setBit(64);
  fp.setBit(130);
  fp.setBit(255);

  std::vector<const ExplicitBitVect *> view{&fp};
  std::size_t actualBits = 0;
  auto packed = packFingerprints(view, actualBits);
  REQUIRE(actualBits == nBits);
  REQUIRE(packed.size() == 4);

  CHECK(((packed[0] >> 0) & 1ULL) == 1ULL);
  CHECK(((packed[0] >> 63) & 1ULL) == 1ULL);
  CHECK(((packed[1] >> 0) & 1ULL) == 1ULL);
  CHECK(((packed[2] >> (130 - 128)) & 1ULL) == 1ULL);
  CHECK(((packed[3] >> (255 - 192)) & 1ULL) == 1ULL);
  CHECK(((packed[0] >> 5) & 1ULL) == 0ULL);
}

TEST_CASE("BulkSimilarity: matches per-pair TanimotoSimilarity reference",
          "[BulkSimilarity]") {
  std::mt19937 rng(0xC0FFEE);
  constexpr std::size_t nBits = 512;
  auto probes = makeRandomFps(7, nBits, rng);
  auto targets = makeRandomFps(11, nBits, rng);
  auto probeViews = asViews(probes);
  auto targetViews = asViews(targets);

  std::vector<double> matrix;
  tanimotoMatrix(probeViews, targetViews, matrix);
  REQUIRE(matrix.size() == probes.size() * targets.size());

  for (std::size_t i = 0; i < probes.size(); ++i) {
    for (std::size_t j = 0; j < targets.size(); ++j) {
      const double expected = TanimotoSimilarity(*probes[i], *targets[j]);
      const double actual = matrix[i * targets.size() + j];
      CHECK(actual == Catch::Approx(expected).margin(1e-12));
    }
  }
}

TEST_CASE("BulkSimilarity: handles fingerprint sizes that aren't 8-word "
          "multiples (exercises AVX-512 tail mask)",
          "[BulkSimilarity]") {
  std::mt19937 rng(0xBADCAFE);
  // 320 bits = 5 uint64_t words; the AVX-512 loop processes 8 words per
  // iteration so this size forces the mask-load tail.
  constexpr std::size_t nBits = 320;
  auto probes = makeRandomFps(4, nBits, rng);
  auto targets = makeRandomFps(6, nBits, rng);
  auto probeViews = asViews(probes);
  auto targetViews = asViews(targets);

  std::vector<double> matrix;
  tanimotoMatrix(probeViews, targetViews, matrix);
  for (std::size_t i = 0; i < probes.size(); ++i) {
    for (std::size_t j = 0; j < targets.size(); ++j) {
      const double expected = TanimotoSimilarity(*probes[i], *targets[j]);
      const double actual = matrix[i * targets.size() + j];
      CHECK(actual == Catch::Approx(expected).margin(1e-12));
    }
  }
}

TEST_CASE("BulkSimilarity: zero fingerprints produce zero similarity",
          "[BulkSimilarity]") {
  constexpr std::size_t nBits = 128;
  ExplicitBitVect a(nBits);
  ExplicitBitVect b(nBits);
  std::vector<const ExplicitBitVect *> probes{&a};
  std::vector<const ExplicitBitVect *> targets{&b};
  std::vector<double> out;
  tanimotoMatrix(probes, targets, out);
  REQUIRE(out.size() == 1);
  CHECK(out[0] == 0.0);
}

TEST_CASE("BulkSimilarity: identical fingerprints score 1.0",
          "[BulkSimilarity]") {
  constexpr std::size_t nBits = 128;
  ExplicitBitVect a(nBits);
  a.setBit(1);
  a.setBit(7);
  a.setBit(99);
  std::vector<const ExplicitBitVect *> probes{&a};
  std::vector<const ExplicitBitVect *> targets{&a};
  std::vector<double> out;
  tanimotoMatrix(probes, targets, out);
  REQUIRE(out.size() == 1);
  CHECK(out[0] == Catch::Approx(1.0));
}

TEST_CASE("BulkSimilarity: empty inputs leave the result vector empty",
          "[BulkSimilarity]") {
  std::vector<const ExplicitBitVect *> probes;
  std::vector<const ExplicitBitVect *> targets;
  std::vector<double> out{1.0, 2.0, 3.0};
  tanimotoMatrix(probes, targets, out);
  CHECK(out.empty());
}

TEST_CASE("BulkSimilarity: mismatched bit sizes are rejected",
          "[BulkSimilarity]") {
  ExplicitBitVect a(128);
  ExplicitBitVect b(256);
  std::vector<const ExplicitBitVect *> probes{&a};
  std::vector<const ExplicitBitVect *> targets{&b};
  std::vector<double> out;
  CHECK_THROWS(tanimotoMatrix(probes, targets, out));
}

TEST_CASE("BulkSimilarity: active kernel name is one of the known values",
          "[BulkSimilarity]") {
  const std::string k = activeKernel();
  CHECK((k == "scalar" || k == "avx512vpopcntdq"));
}
