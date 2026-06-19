# BulkSimilarity

Bulk Tanimoto similarity over many `ExplicitBitVect` fingerprints with a
multi-threaded, optionally AVX-512-accelerated, CPU kernel.

## What it does

`tanimotoMatrix` computes an `M x N` matrix of Tanimoto coefficients between
two sets of `ExplicitBitVect` fingerprints. The result is written into a
caller-allocated row-major `double` buffer of length `M * N`, where entry
`[i, j]` is the Tanimoto coefficient between probe `i` and target `j`.

All fingerprints in a single call must share the same bit length, and that
length must be a positive multiple of 64.

The packed-buffer overload lets callers pre-pack their fingerprints once
and reuse them across many calls; the convenience overload accepts
`std::vector<const ExplicitBitVect*>` directly.

## C++ usage

```cpp
#include <DataStructs/BulkSimilarity/BulkSimilarity.h>
#include <DataStructs/ExplicitBitVect.h>

std::vector<const ExplicitBitVect*> probes  = ...;  // M fingerprints
std::vector<const ExplicitBitVect*> targets = ...;  // N fingerprints

std::vector<double> matrix;  // length M * N on return, row-major
RDKit::BulkSimilarity::tanimotoMatrix(probes, targets, matrix);

// matrix[i * targets.size() + j] is Tanimoto(probes[i], targets[j])
```

For a hot path that screens the same targets against many probe batches,
pre-pack once:

```cpp
std::size_t bits = 0;
auto packedTargets = RDKit::BulkSimilarity::packFingerprints(targets, bits);
const std::size_t words = RDKit::BulkSimilarity::wordsForBits(bits);

// ... later, possibly many times ...
RDKit::BulkSimilarity::tanimotoMatrix(
    packedProbes.data(), probes.size(),
    packedTargets.data(), targets.size(),
    words, out.data());
```

## Python usage

```python
from rdkit.DataStructs import BulkTanimotoMatrix

# probes, targets: lists of ExplicitBitVect (all same bit length)
matrix = BulkTanimotoMatrix(probes, targets)   # numpy float64 (M, N)
```

## Kernel design

* The fingerprints are kept in a row-major packed layout of `uint64_t`
  words.
* The inner loop is parallelised across probe rows. The worker count is
  chosen from the matrix size: the kernel only spawns threads once there is
  enough work to amortize thread-creation cost (~2M word-ops per worker),
  and caps the count so each worker gets a coarse, contiguous block of probe
  rows. Small matrices run single-threaded. The probe-row popcount is hoisted
  out of the target loop so each `(probe, target)` pair costs two popcounts
  instead of three.
* Two implementations of the popcount step are provided:
  * a **scalar** path using `__builtin_popcountll` / `__popcnt64`;
  * an **AVX-512** path using `_mm512_popcnt_epi64`, which processes eight
    `uint64_t` words per instruction (Intel Ice Lake and later, AMD Zen 4
    and later).
* Selection is **runtime**: at first call, `CPUID` + `XGETBV` confirm both
  hardware and OS support for AVX-512F + `AVX512_VPOPCNTDQ`, and the
  dispatcher swaps two function pointers. The compile-time check in the
  module's `CMakeLists.txt` decides whether the AVX-512 translation unit
  is included at build time, so prebuilt distributions (conda-forge,
  manylinux) ship a binary that runs everywhere and uses AVX-512 only on
  hardware that has it.
* The tail of the AVX-512 loop uses `_mm512_maskz_loadu_epi64` so any
  remaining 0-7 words still go through the SIMD path.

`activeKernel()` returns `"scalar"` or `"avx512vpopcntdq"` for diagnostics.

## Build

The module is built by default. To force the scalar-only build (for
testing or to keep the binary smaller), pass
`-DRDK_BULKSIM_COMPILER_HAS_AVX512_VPOPCNTQ=OFF` to CMake; otherwise the
build system probes the compiler and includes the AVX-512 translation
unit when it accepts the necessary flags
(`-mavx512f -mavx512vpopcntdq` on GCC/Clang, `/arch:AVX512` on MSVC).

## Future work

* GPU backend (CUDA `__popcll`-based kernel for an M x N matrix). A
  prototype lives on the `feature/cuda-tanimoto-backend` branch.
* AVX2 Harley-Seal fallback for the gap between scalar and AVX-512.
* `SparseBitVect` overload using a different packing strategy.
* `SparseIntVect`-based "count" Tanimoto matrix variant.
