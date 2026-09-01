//
//  Copyright (C) 2026 Clay Moore and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  CPU implementation of the bulk Tanimoto similarity kernel.
//
//  The fingerprints are kept in a row-major packed layout of uint64_t words
//  so the inner loop is a sequence of popcount-AND and popcount operations.
//  The popcount step dispatches at startup to one of two kernels:
//
//    * an AVX-512 implementation that uses VPOPCNTQ to popcount eight
//      uint64_t lanes per instruction, defined in BulkSimilarityCpuAvx512.cpp
//      and compiled with the appropriate -m / /arch flags;
//    * a scalar implementation defined in this file that uses the standard
//      64-bit popcount intrinsic.

#include "BulkSimilarity.h"

#include <DataStructs/ExplicitBitVect.h>
#include <RDGeneral/Exceptions.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <thread>
#include <vector>

#if defined(_MSC_VER)
#include <intrin.h>
#endif
#if defined(__GNUC__) || defined(__clang__)
#include <cpuid.h>
#endif

namespace RDKit {
namespace BulkSimilarity {

#ifdef RDK_BULKSIM_HAS_AVX512_VPOPCNTQ
namespace detail {
unsigned int popcountWordsAvx512(const std::uint64_t *data, std::size_t words);
void intersectAndTargetPopAvx512(const std::uint64_t *probe,
                                 const std::uint64_t *target,
                                 std::size_t words, unsigned int &interOut,
                                 unsigned int &targetOut);
}  // namespace detail
#endif

namespace {

inline unsigned int popcount64(std::uint64_t x) {
#if defined(_MSC_VER)
  return static_cast<unsigned int>(__popcnt64(x));
#else
  return static_cast<unsigned int>(__builtin_popcountll(x));
#endif
}

unsigned int popcountWordsScalar(const std::uint64_t *data,
                                 std::size_t words) {
  unsigned int total = 0;
  for (std::size_t w = 0; w < words; ++w) {
    total += popcount64(data[w]);
  }
  return total;
}

void intersectAndTargetPopScalar(const std::uint64_t *probe,
                                 const std::uint64_t *target,
                                 std::size_t words, unsigned int &interOut,
                                 unsigned int &targetOut) {
  unsigned int inter = 0;
  unsigned int targ = 0;
  for (std::size_t w = 0; w < words; ++w) {
    const std::uint64_t a = probe[w];
    const std::uint64_t b = target[w];
    targ += popcount64(b);
    inter += popcount64(a & b);
  }
  interOut = inter;
  targetOut = targ;
}

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || \
    defined(_M_IX86)
#define RDK_BULKSIM_HAVE_CPUID 1
#endif

#ifdef RDK_BULKSIM_HAVE_CPUID
inline void cpuidCount(int leaf, int subleaf, int regs[4]) {
#if defined(_MSC_VER)
  __cpuidex(regs, leaf, subleaf);
#else
  unsigned int eax = 0, ebx = 0, ecx = 0, edx = 0;
  __cpuid_count(leaf, subleaf, eax, ebx, ecx, edx);
  regs[0] = static_cast<int>(eax);
  regs[1] = static_cast<int>(ebx);
  regs[2] = static_cast<int>(ecx);
  regs[3] = static_cast<int>(edx);
#endif
}

inline unsigned long long readXcr0() {
#if defined(_MSC_VER)
  return _xgetbv(0);
#else
  unsigned int eax = 0, edx = 0;
  __asm__ volatile("xgetbv" : "=a"(eax), "=d"(edx) : "c"(0));
  return (static_cast<unsigned long long>(edx) << 32) | eax;
#endif
}
#endif  // RDK_BULKSIM_HAVE_CPUID

bool cpuSupportsAvx512VpopcntDq() {
#ifdef RDK_BULKSIM_HAVE_CPUID
  int regs[4];
  cpuidCount(0, 0, regs);
  if (regs[0] < 7) {
    return false;
  }
  cpuidCount(1, 0, regs);
  const bool osxsave = (regs[2] & (1 << 27)) != 0;
  if (!osxsave) {
    return false;
  }
  const unsigned long long xcr0 = readXcr0();
  // AVX-512 needs OS save/restore of opmask, ZMM_Hi256 and Hi16_ZMM.
  if ((xcr0 & 0xE0ULL) != 0xE0ULL) {
    return false;
  }
  cpuidCount(7, 0, regs);
  const bool avx512f = (regs[1] & (1 << 16)) != 0;       // EBX bit 16
  const bool vpopcntdq = (regs[2] & (1 << 14)) != 0;     // ECX bit 14
  return avx512f && vpopcntdq;
#else
  return false;
#endif
}

using PopcountFn = unsigned int (*)(const std::uint64_t *, std::size_t);
using IntersectFn = void (*)(const std::uint64_t *, const std::uint64_t *,
                             std::size_t, unsigned int &, unsigned int &);

struct KernelTable {
  PopcountFn popcount;
  IntersectFn intersect;
  const char *name;
};

KernelTable pickKernel() {
#ifdef RDK_BULKSIM_HAS_AVX512_VPOPCNTQ
  if (cpuSupportsAvx512VpopcntDq()) {
    return {&detail::popcountWordsAvx512, &detail::intersectAndTargetPopAvx512,
            "avx512vpopcntdq"};
  }
#endif
  return {&popcountWordsScalar, &intersectAndTargetPopScalar, "scalar"};
}

const KernelTable &kernel() {
  static const KernelTable k = pickKernel();
  return k;
}

}  // namespace

std::string activeKernel() { return kernel().name; }

std::size_t wordsForBits(std::size_t numBits) {
  if (numBits == 0 || (numBits % 64) != 0) {
    throw ValueErrorException(
        "BulkSimilarity requires fingerprint bit length to be a positive "
        "multiple of 64");
  }
  return numBits / 64;
}

std::vector<std::uint64_t> packFingerprints(
    const std::vector<const ExplicitBitVect *> &fps, std::size_t &numBits) {
  if (fps.empty()) {
    numBits = 0;
    return {};
  }
  if (fps.front() == nullptr || fps.front()->dp_bits == nullptr) {
    throw ValueErrorException("null fingerprint in BulkSimilarity input");
  }
  numBits = fps.front()->getNumBits();
  const std::size_t words = wordsForBits(numBits);

  std::vector<std::uint64_t> packed(fps.size() * words, 0);
  using block_t = boost::dynamic_bitset<>::block_type;
  constexpr std::size_t bitsPerBlock = sizeof(block_t) * 8;
  static_assert(bitsPerBlock == 32 || bitsPerBlock == 64,
                "Unexpected dynamic_bitset block size");

  std::vector<block_t> scratch;
  for (std::size_t k = 0; k < fps.size(); ++k) {
    const ExplicitBitVect *fp = fps[k];
    if (fp == nullptr || fp->dp_bits == nullptr) {
      throw ValueErrorException("null fingerprint in BulkSimilarity input");
    }
    if (fp->getNumBits() != numBits) {
      throw ValueErrorException(
          "all fingerprints must share the same bit length");
    }
    const std::size_t numBlocks = fp->dp_bits->num_blocks();
    if (scratch.size() < numBlocks) {
      scratch.resize(numBlocks);
    }
    boost::to_block_range(*fp->dp_bits, scratch.data());

    std::uint64_t *dest = packed.data() + k * words;
    if constexpr (bitsPerBlock == 64) {
      std::memcpy(dest, scratch.data(), words * sizeof(std::uint64_t));
    } else {
      for (std::size_t w = 0; w < words; ++w) {
        const std::uint64_t lo =
            (2 * w < numBlocks) ? static_cast<std::uint64_t>(scratch[2 * w])
                                : 0ULL;
        const std::uint64_t hi =
            (2 * w + 1 < numBlocks)
                ? static_cast<std::uint64_t>(scratch[2 * w + 1])
                : 0ULL;
        dest[w] = lo | (hi << 32);
      }
    }
  }
  return packed;
}

void tanimotoMatrix(const std::uint64_t *probes, std::size_t numProbes,
                    const std::uint64_t *targets, std::size_t numTargets,
                    std::size_t words, double *out) {
  if (numProbes == 0 || numTargets == 0) {
    return;
  }
  if (probes == nullptr || targets == nullptr || out == nullptr) {
    throw ValueErrorException("null buffer passed to tanimotoMatrix");
  }
  if (words == 0) {
    throw ValueErrorException("words must be > 0");
  }

  const KernelTable &k = kernel();
  const PopcountFn popcountFn = k.popcount;
  const IntersectFn intersectFn = k.intersect;

  unsigned int hwThreads = std::thread::hardware_concurrency();
  if (hwThreads == 0) {
    hwThreads = 1;
  }
  // Decide how many worker threads to use. Spawning a std::thread costs on
  // the order of tens of microseconds, so threading only pays off once each
  // worker has enough word-level popcount work to amortize that cost. Cap the
  // worker count so every worker gets at least kMinWordOpsPerThread words of
  // work; for small matrices this collapses to a single thread and avoids the
  // spawn overhead entirely, which would otherwise dominate and can make a
  // naive one-thread-per-row split slower than a plain serial loop.
  constexpr std::size_t kMinWordOpsPerThread = std::size_t{1} << 21;  // ~2M
  const std::size_t totalWordOps = numProbes * numTargets * words;
  std::size_t workersByWork = totalWordOps / kMinWordOpsPerThread;
  if (workersByWork < 1) {
    workersByWork = 1;
  }
  const std::size_t numWorkers =
      std::min({static_cast<std::size_t>(hwThreads), numProbes, workersByWork});

  // Precompute each probe's popcount once. It is reused for every target and
  // every target tile, so there is no need to recompute it in the inner loop.
  std::vector<unsigned int> probePops(numProbes);
  for (std::size_t i = 0; i < numProbes; ++i) {
    probePops[i] = popcountFn(probes + i * words, words);
  }

  // Cache-block the targets. Streaming the whole target set once per probe row
  // makes wide matrices memory-bandwidth bound: the target data no longer fits
  // in cache and is refetched from memory for every probe row. Instead, walk
  // the targets in tiles small enough to stay resident in L2 and reuse each
  // tile across all of a worker's probe rows, which cuts target memory traffic
  // by roughly the number of probe rows each worker handles.
  constexpr std::size_t kTargetTileBytes = std::size_t{128} << 10;  // 128 KiB
  const std::size_t bytesPerTarget = words * sizeof(std::uint64_t);
  const std::size_t targetTile =
      std::max<std::size_t>(1, kTargetTileBytes / bytesPerTarget);

  auto worker = [&](std::size_t probeStart, std::size_t probeEnd) {
    for (std::size_t jTile = 0; jTile < numTargets; jTile += targetTile) {
      const std::size_t jEnd = std::min(jTile + targetTile, numTargets);
      for (std::size_t i = probeStart; i < probeEnd; ++i) {
        const std::uint64_t *pRow = probes + i * words;
        const unsigned int probePop = probePops[i];
        double *outRow = out + i * numTargets;
        for (std::size_t j = jTile; j < jEnd; ++j) {
          const std::uint64_t *tRow = targets + j * words;
          unsigned int interPop = 0;
          unsigned int targetPop = 0;
          intersectFn(pRow, tRow, words, interPop, targetPop);
          const unsigned int unionPop = probePop + targetPop - interPop;
          outRow[j] =
              unionPop == 0 ? 0.0 : static_cast<double>(interPop) / unionPop;
        }
      }
    }
  };

  if (numWorkers == 1) {
    worker(0, numProbes);
    return;
  }

  std::vector<std::thread> threads;
  threads.reserve(numWorkers - 1);
  const std::size_t chunk = (numProbes + numWorkers - 1) / numWorkers;
  for (std::size_t t = 0; t + 1 < numWorkers; ++t) {
    const std::size_t start = t * chunk;
    const std::size_t end = std::min(start + chunk, numProbes);
    if (start >= end) {
      break;
    }
    threads.emplace_back(worker, start, end);
  }
  const std::size_t lastStart = (numWorkers - 1) * chunk;
  if (lastStart < numProbes) {
    worker(lastStart, numProbes);
  }
  for (auto &th : threads) {
    th.join();
  }
}

void tanimotoMatrix(const std::vector<const ExplicitBitVect *> &probes,
                    const std::vector<const ExplicitBitVect *> &targets,
                    std::vector<double> &out) {
  out.assign(probes.size() * targets.size(), 0.0);
  if (probes.empty() || targets.empty()) {
    return;
  }
  std::size_t probeBits = 0;
  std::size_t targetBits = 0;
  std::vector<std::uint64_t> packedProbes = packFingerprints(probes, probeBits);
  std::vector<std::uint64_t> packedTargets =
      packFingerprints(targets, targetBits);
  if (probeBits != targetBits) {
    throw ValueErrorException(
        "probe and target fingerprints must have the same bit length");
  }
  const std::size_t words = wordsForBits(probeBits);
  tanimotoMatrix(packedProbes.data(), probes.size(), packedTargets.data(),
                 targets.size(), words, out.data());
}

}  // namespace BulkSimilarity
}  // namespace RDKit
