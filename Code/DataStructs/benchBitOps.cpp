// benchBitOps.cpp - Performance benchmark for SIMD-accelerated bitmap operations
//
// Build as part of RDKit (add to CMakeLists.txt) or standalone:
//   g++ -O2 -std=c++17 -mavx2 -DRDK_OPTIMIZE_POPCNT -mpopcnt -I../../ \
//       benchBitOps.cpp BitOps.cpp -o benchBitOps
//
// Reports throughput (million pairs/sec) for CalcBitmapTanimoto, Dice, Tversky,
// Popcount, and AllProbeBitsMatch at standard fingerprint sizes.

#include <chrono>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

#include "BitOps.h"

// Detect SIMD level for reporting
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
#define BENCH_X86
#ifdef _MSC_VER
#include <intrin.h>
#else
#include <immintrin.h>
#endif
#endif

static const char *detect_simd_level() {
#ifdef BENCH_X86
#ifdef _MSC_VER
  int info[4];
  __cpuid(info, 1);
  bool osxsave = (info[2] >> 27) & 1;
  if (osxsave) {
    unsigned long long xcr0 = _xgetbv(0);
    if ((xcr0 & 0x6) == 0x6) {
      __cpuidex(info, 7, 0);
      if ((info[1] >> 5) & 1) return "AVX2";
    }
  }
  if ((info[2] >> 9) & 1) return "SSSE3";
  if ((info[2] >> 23) & 1) return "POPCNT-only";
#elif defined(__GNUC__) || defined(__clang__)
  __builtin_cpu_init();
  if (__builtin_cpu_supports("avx2")) return "AVX2";
  if (__builtin_cpu_supports("ssse3")) return "SSSE3";
  if (__builtin_cpu_supports("popcnt")) return "POPCNT-only";
#endif
#endif
  return "scalar (no SIMD)";
}

struct BenchResult {
  double best_ms;
  double mpairs_sec;
  double checksum;
};

template <typename Func>
static BenchResult bench(Func fn, int trials, long long total_ops) {
  double best = 1e30;
  double checksum = 0;
  for (int t = 0; t < trials; t++) {
    auto t0 = std::chrono::high_resolution_clock::now();
    checksum = fn();
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    if (ms < best) best = ms;
  }
  return {best, total_ops / best / 1e3, checksum};
}

int main() {
  const int NUM_FPS = 10000;
  const int REPS = 100;
  const int TRIALS = 5;

  printf("=== RDKit Bitmap Similarity Benchmark ===\n");
  printf("Active SIMD level: %s\n", detect_simd_level());
  printf("Database: %d fingerprints, %d repetitions, %d trials (best-of)\n\n",
         NUM_FPS, REPS, TRIALS);

  // Test at common fingerprint sizes
  const int fp_bits[] = {1024, 2048, 4096};
  const int num_sizes = 3;

  std::mt19937 rng(42);

  for (int s = 0; s < num_sizes; s++) {
    const int fp_bytes = fp_bits[s] / 8;
    printf("--- %d-bit fingerprints (%d bytes) ---\n", fp_bits[s], fp_bytes);

    // Generate random fingerprints (~50% density, typical for Morgan FPs)
    std::vector<std::vector<unsigned char>> fps(NUM_FPS);
    for (auto &fp : fps) {
      fp.resize(fp_bytes);
      for (int j = 0; j < fp_bytes; j++) fp[j] = rng() & 0xFF;
    }

    const unsigned char *query = fps[0].data();
    long long total_pairs = (long long)REPS * NUM_FPS;

    // Tanimoto
    auto r = bench(
        [&]() {
          double sum = 0;
          for (int rep = 0; rep < REPS; rep++)
            for (int i = 0; i < NUM_FPS; i++)
              sum += CalcBitmapTanimoto(query, fps[i].data(), fp_bytes);
          return sum;
        },
        TRIALS, total_pairs);
    printf("  Tanimoto:     %8.1f ms  %8.1f M pairs/sec  [checksum=%.2f]\n",
           r.best_ms, r.mpairs_sec, r.checksum);

    // Dice
    r = bench(
        [&]() {
          double sum = 0;
          for (int rep = 0; rep < REPS; rep++)
            for (int i = 0; i < NUM_FPS; i++)
              sum += CalcBitmapDice(query, fps[i].data(), fp_bytes);
          return sum;
        },
        TRIALS, total_pairs);
    printf("  Dice:         %8.1f ms  %8.1f M pairs/sec  [checksum=%.2f]\n",
           r.best_ms, r.mpairs_sec, r.checksum);

    // Tversky (a=0.7, b=0.3)
    r = bench(
        [&]() {
          double sum = 0;
          for (int rep = 0; rep < REPS; rep++)
            for (int i = 0; i < NUM_FPS; i++)
              sum += CalcBitmapTversky(query, fps[i].data(), fp_bytes, 0.7, 0.3);
          return sum;
        },
        TRIALS, total_pairs);
    printf("  Tversky:      %8.1f ms  %8.1f M pairs/sec  [checksum=%.2f]\n",
           r.best_ms, r.mpairs_sec, r.checksum);

    // Popcount
    r = bench(
        [&]() {
          double sum = 0;
          for (int rep = 0; rep < REPS; rep++)
            for (int i = 0; i < NUM_FPS; i++)
              sum += CalcBitmapPopcount(fps[i].data(), fp_bytes);
          return sum;
        },
        TRIALS, total_pairs);
    printf("  Popcount:     %8.1f ms  %8.1f M ops/sec    [checksum=%.0f]\n",
           r.best_ms, r.mpairs_sec, r.checksum);

    // AllProbeBitsMatch
    r = bench(
        [&]() {
          double sum = 0;
          for (int rep = 0; rep < REPS; rep++)
            for (int i = 0; i < NUM_FPS; i++)
              sum += CalcBitmapAllProbeBitsMatch(query, fps[i].data(), fp_bytes);
          return sum;
        },
        TRIALS, total_pairs);
    printf("  ProbeMatch:   %8.1f ms  %8.1f M ops/sec    [checksum=%.0f]\n",
           r.best_ms, r.mpairs_sec, r.checksum);

    printf("\n");
  }

  return 0;
}
