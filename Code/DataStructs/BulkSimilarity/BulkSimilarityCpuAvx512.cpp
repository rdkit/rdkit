//
//  Copyright (C) 2026 Clay Moore and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  AVX-512 VPOPCNTQ popcount kernels for the BulkSimilarity library.
//
//  This translation unit is compiled with AVX-512 flags (see the module's
//  CMakeLists.txt) and is only included in the library when the host
//  compiler supports those flags. The symbols defined here are only ever
//  called by the runtime dispatcher in BulkSimilarityCpu.cpp after a
//  CPUID + XGETBV check confirms the CPU and the OS support
//  AVX-512F + AVX512_VPOPCNTDQ. Keep this file small so the compiler does
//  not emit AVX-512 instructions outside of these functions.

#include "BulkSimilarity.h"

#include <immintrin.h>

#include <cstddef>
#include <cstdint>

namespace RDKit {
namespace BulkSimilarity {
namespace detail {

unsigned int popcountWordsAvx512(const std::uint64_t *data,
                                 std::size_t words) {
  std::size_t w = 0;
  __m512i acc = _mm512_setzero_si512();
  for (; w + 8 <= words; w += 8) {
    const __m512i v = _mm512_loadu_si512(
        reinterpret_cast<const __m512i *>(data + w));
    acc = _mm512_add_epi64(acc, _mm512_popcnt_epi64(v));
  }
  unsigned int total =
      static_cast<unsigned int>(_mm512_reduce_add_epi64(acc));

  // Mask-load the tail (0-7 trailing words) so we stay in AVX-512 for the
  // whole loop and avoid leaking a scalar dependency back into the caller.
  const std::size_t tail = words - w;
  if (tail != 0) {
    const __mmask8 mask = static_cast<__mmask8>((1u << tail) - 1u);
    const __m512i v = _mm512_maskz_loadu_epi64(
        mask, reinterpret_cast<const __m512i *>(data + w));
    total += static_cast<unsigned int>(
        _mm512_reduce_add_epi64(_mm512_popcnt_epi64(v)));
  }
  return total;
}

void intersectAndTargetPopAvx512(const std::uint64_t *probe,
                                 const std::uint64_t *target,
                                 std::size_t words, unsigned int &interOut,
                                 unsigned int &targetOut) {
  std::size_t w = 0;
  __m512i interAcc = _mm512_setzero_si512();
  __m512i targAcc = _mm512_setzero_si512();
  for (; w + 8 <= words; w += 8) {
    const __m512i a = _mm512_loadu_si512(
        reinterpret_cast<const __m512i *>(probe + w));
    const __m512i b = _mm512_loadu_si512(
        reinterpret_cast<const __m512i *>(target + w));
    targAcc = _mm512_add_epi64(targAcc, _mm512_popcnt_epi64(b));
    interAcc =
        _mm512_add_epi64(interAcc,
                         _mm512_popcnt_epi64(_mm512_and_si512(a, b)));
  }
  unsigned int inter =
      static_cast<unsigned int>(_mm512_reduce_add_epi64(interAcc));
  unsigned int targ =
      static_cast<unsigned int>(_mm512_reduce_add_epi64(targAcc));

  const std::size_t tail = words - w;
  if (tail != 0) {
    const __mmask8 mask = static_cast<__mmask8>((1u << tail) - 1u);
    const __m512i a = _mm512_maskz_loadu_epi64(
        mask, reinterpret_cast<const __m512i *>(probe + w));
    const __m512i b = _mm512_maskz_loadu_epi64(
        mask, reinterpret_cast<const __m512i *>(target + w));
    targ += static_cast<unsigned int>(
        _mm512_reduce_add_epi64(_mm512_popcnt_epi64(b)));
    inter += static_cast<unsigned int>(
        _mm512_reduce_add_epi64(_mm512_popcnt_epi64(_mm512_and_si512(a, b))));
  }
  interOut = inter;
  targetOut = targ;
}

}  // namespace detail
}  // namespace BulkSimilarity
}  // namespace RDKit
