//
//  Copyright (C) 2026 Clay Moore and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file BulkSimilarity.h

  \brief Bulk Tanimoto similarity over many fingerprints.

  Computes an M x N matrix of Tanimoto coefficients between two sets of
  ExplicitBitVect fingerprints. The kernel is multi-threaded across probe
  rows and dispatches at runtime to an AVX-512 popcount path on CPUs that
  support \c VPOPCNTDQ (Intel Ice Lake and later, AMD Zen 4 and later);
  otherwise it uses the scalar 64-bit popcount intrinsic.

  All fingerprints in a single call must share the same bit length, and
  that length must be a positive multiple of 64.
*/

#include <RDGeneral/export.h>
#ifndef RDKIT_BULK_SIMILARITY_H
#define RDKIT_BULK_SIMILARITY_H

#include <cstdint>
#include <string>
#include <vector>

class ExplicitBitVect;

namespace RDKit {
namespace BulkSimilarity {

//! Number of \c uint64_t words required to store a fingerprint of
//! \c numBits bits. Throws if \c numBits is not a positive multiple of 64.
RDKIT_BULKSIMILARITY_EXPORT std::size_t wordsForBits(std::size_t numBits);

//! Pack a vector of fingerprints into a contiguous, row-major
//! \c uint64_t buffer.
/*!
  All fingerprints must share the same bit length, and that length must be a
  positive multiple of 64. Bit \c i of fingerprint \c k is stored in bit
  <tt>i % 64</tt> of word <tt>k * words + i / 64</tt>, where
  <tt>words = numBits / 64</tt>.

  \param fps        the fingerprints to pack; pointers must be non-null
  \param numBits    output: the bit length shared by all fingerprints
                    (set to 0 when \c fps is empty)
  \return the packed buffer, of length <tt>fps.size() * (numBits / 64)</tt>
*/
RDKIT_BULKSIMILARITY_EXPORT std::vector<std::uint64_t> packFingerprints(
    const std::vector<const ExplicitBitVect *> &fps, std::size_t &numBits);

//! Compute an M x N Tanimoto similarity matrix from pre-packed buffers.
/*!
  \param probes      packed probe fingerprints, of length
                     <tt>numProbes * words</tt>
  \param numProbes   number of probe fingerprints (M)
  \param targets     packed target fingerprints, of length
                     <tt>numTargets * words</tt>
  \param numTargets  number of target fingerprints (N)
  \param words       number of \c uint64_t words per fingerprint
  \param out         row-major M x N output buffer, of length
                     <tt>numProbes * numTargets</tt>, allocated by the caller.
                     Entry <tt>[i, j]</tt> receives the Tanimoto coefficient
                     between probe \c i and target \c j. If both
                     fingerprints are all-zero, the result is 0.0.

  When \c numProbes or \c numTargets is zero the call is a no-op.
*/
RDKIT_BULKSIMILARITY_EXPORT void tanimotoMatrix(const std::uint64_t *probes,
                                                std::size_t numProbes,
                                                const std::uint64_t *targets,
                                                std::size_t numTargets,
                                                std::size_t words,
                                                double *out);

//! \overload Convenience wrapper that packs the inputs and resizes \c out.
RDKIT_BULKSIMILARITY_EXPORT void tanimotoMatrix(
    const std::vector<const ExplicitBitVect *> &probes,
    const std::vector<const ExplicitBitVect *> &targets,
    std::vector<double> &out);

//! Name of the CPU popcount kernel selected at startup: \c "avx512vpopcntdq"
//! when the AVX-512 path is active, otherwise \c "scalar".
RDKIT_BULKSIMILARITY_EXPORT std::string activeKernel();

}  // namespace BulkSimilarity
}  // namespace RDKit

#endif
