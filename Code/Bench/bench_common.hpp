#pragma once

#include <cstdint>
#include <cstddef>

namespace bench_common {
constexpr const char *CASES[] = {
    "COC1/C=C/OC2(C)Oc3c(C)c(O)c4c(O)c(c(/C=N/OC(c5ccccc5)c5ccccc5)cc4c3C2=O)NC(=O)/C(C)=C\\C=C\\C(C)C(O)C(C)C(O)C(C)C(OC(C)=O)C1C",
    "Cn1cnc2n(C)c(=O)n(C)c(=O)c12",
    "c1ccc2c(c1)c3ccccc3c4ccccc24",
    "F[C@@H](Cl)[C@H](Br)[C@](I)(O)[C@@](N)(C#N)C(=O)O",
    "C[S+](C)(C)[O-]",
};

inline uint64_t nth_random(uint64_t block) noexcept {
  // https://en.wikipedia.org/wiki/Feistel_cipher

  constexpr uint32_t K[4] = {
      0xa1509c6bu,
      0x92a7d026u,
      0xf6486c80u,
      0x4ca7238bu,
  };
  constexpr size_t ROUNDS = 24;

  auto rotl = [](uint32_t x, unsigned r) noexcept {
    r &= 31u;
    return (x << r) | (x >> (32 - r));
  };
  auto mix = [](uint32_t v) noexcept {
    v ^= v >> 16;
    v *= 0x7feb352du;
    v ^= v >> 15;
    v *= 0x846ca68bu;
    v ^= v >> 16;
    return v;
  };
  auto subkey = [&](size_t i) noexcept {
    uint32_t base = K[i & 3] + uint32_t(i) * 0x9E3779B9u;
    return rotl(base, uint32_t(i));
  };

  auto L = uint32_t(block);
  auto R = uint32_t(block >> 32);

  for (size_t i = 0; i < ROUNDS; ++i) {
    uint32_t f = mix(R ^ subkey(i));
    uint32_t T = L ^ f;
    L = R;
    R = T;
  }

  return (uint64_t(R) << 32) | L;
}
}  // namespace bench_common
