//
//  Copyright (C) 2026 RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "CharacteristicPolynomial.h"

#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>

#if defined(_MSC_VER) && !defined(__clang__) && \
    (defined(_M_X64) || defined(_M_ARM64))
#include <intrin.h>
#endif

namespace RDKit {
namespace Descriptors {
namespace detail {
namespace {

using BigInt = boost::multiprecision::cpp_int;

// The 64 largest primes below 2^62. Residues modulo all of them determine any
// integer of up to about 3,960 bits.
constexpr std::array<std::uint64_t, 64> largePrimes = {
    0x3fffffffffffffc7ull, 0x3fffffffffffffa9ull, 0x3fffffffffffff8bull,
    0x3fffffffffffff71ull, 0x3fffffffffffff67ull, 0x3fffffffffffff59ull,
    0x3fffffffffffff55ull, 0x3fffffffffffff3dull, 0x3fffffffffffff35ull,
    0x3ffffffffffffeefull, 0x3ffffffffffffee1ull, 0x3ffffffffffffec3ull,
    0x3ffffffffffffe45ull, 0x3ffffffffffffe1dull, 0x3ffffffffffffe11ull,
    0x3ffffffffffffdc1ull, 0x3ffffffffffffdbbull, 0x3ffffffffffffda5ull,
    0x3ffffffffffffd87ull, 0x3ffffffffffffd69ull, 0x3ffffffffffffd03ull,
    0x3ffffffffffffcfbull, 0x3ffffffffffffcf7ull, 0x3ffffffffffffce9ull,
    0x3ffffffffffffcd3ull, 0x3ffffffffffffcc1ull, 0x3ffffffffffffc65ull,
    0x3ffffffffffffc2bull, 0x3ffffffffffffc1full, 0x3ffffffffffffc17ull,
    0x3ffffffffffffc11ull, 0x3ffffffffffffc07ull, 0x3ffffffffffffb53ull,
    0x3ffffffffffffb27ull, 0x3ffffffffffffaf3ull, 0x3ffffffffffffab7ull,
    0x3ffffffffffffa67ull, 0x3ffffffffffffa15ull, 0x3ffffffffffff9efull,
    0x3ffffffffffff9d9ull, 0x3ffffffffffff9d3ull, 0x3ffffffffffff9c5ull,
    0x3ffffffffffff9afull, 0x3ffffffffffff977ull, 0x3ffffffffffff95full,
    0x3ffffffffffff95bull, 0x3ffffffffffff959ull, 0x3ffffffffffff8e1ull,
    0x3ffffffffffff8a7ull, 0x3ffffffffffff889ull, 0x3ffffffffffff87dull,
    0x3ffffffffffff805ull, 0x3ffffffffffff7e7ull, 0x3ffffffffffff7c9ull,
    0x3ffffffffffff7a3ull, 0x3ffffffffffff775ull, 0x3ffffffffffff757ull,
    0x3ffffffffffff739ull, 0x3ffffffffffff713ull, 0x3ffffffffffff6d1ull,
    0x3ffffffffffff6c1ull, 0x3ffffffffffff6b9ull, 0x3ffffffffffff6a3ull,
    0x3ffffffffffff68bull,
};

//! The high and low words of a * b.
inline void multiplyWide(std::uint64_t a, std::uint64_t b, std::uint64_t &high,
                         std::uint64_t &low) {
#if defined(__SIZEOF_INT128__)
  const unsigned __int128 product = static_cast<unsigned __int128>(a) * b;
  low = static_cast<std::uint64_t>(product);
  high = static_cast<std::uint64_t>(product >> 64);
#elif defined(_MSC_VER) && !defined(__clang__) && \
    (defined(_M_X64) || defined(_M_ARM64))
  low = a * b;
  high = __umulh(a, b);
#else
  // schoolbook multiplication in 32-bit halves
  const std::uint64_t aLow = a & 0xffffffffu, aHigh = a >> 32;
  const std::uint64_t bLow = b & 0xffffffffu, bHigh = b >> 32;
  const std::uint64_t lowLow = aLow * bLow, highLow = aHigh * bLow;
  const std::uint64_t lowHigh = aLow * bHigh, highHigh = aHigh * bHigh;
  const std::uint64_t middle =
      (lowLow >> 32) + (highLow & 0xffffffffu) + (lowHigh & 0xffffffffu);
  low = (middle << 32) | (lowLow & 0xffffffffu);
  high = highHigh + (highLow >> 32) + (lowHigh >> 32) + (middle >> 32);
#endif
}

//! Arithmetic modulo an odd prime p < 2^62 in Montgomery form (R = 2^64).
class MontgomeryField {
 public:
  explicit MontgomeryField(std::uint64_t prime) : d_p(prime) {
    // Newton's iteration for p^-1 mod 2^64; each step doubles the valid bits
    std::uint64_t inverse = d_p;
    for (int i = 0; i < 6; ++i) {
      inverse *= 2 - d_p * inverse;
    }
    d_negInverse = ~inverse + 1;
    // 2^64 mod p is the Montgomery form of 1; doubling it 64 more times gives
    // R^2 mod p, which converts into Montgomery form.
    d_one = (~d_p + 1) % d_p;
    d_rSquared = d_one;
    for (int i = 0; i < 64; ++i) {
      d_rSquared = add(d_rSquared, d_rSquared);
    }
  }

  std::uint64_t prime() const { return d_p; }
  std::uint64_t one() const { return d_one; }
  std::uint64_t to(std::uint64_t a) const { return mul(a % d_p, d_rSquared); }
  std::uint64_t from(std::uint64_t a) const { return reduce(0, a); }
  std::uint64_t add(std::uint64_t a, std::uint64_t b) const {
    const std::uint64_t s = a + b;  // < 2^63, no overflow
    return s >= d_p ? s - d_p : s;
  }
  std::uint64_t sub(std::uint64_t a, std::uint64_t b) const {
    return a >= b ? a - b : a + d_p - b;
  }
  std::uint64_t mul(std::uint64_t a, std::uint64_t b) const {
    std::uint64_t high, low;
    multiplyWide(a, b, high, low);
    return reduce(high, low);
  }
  //! a^-1 by Fermat's little theorem, a^(p-2)
  std::uint64_t inverse(std::uint64_t a) const {
    std::uint64_t result = d_one;
    for (std::uint64_t e = d_p - 2; e; e >>= 1, a = mul(a, a)) {
      if (e & 1) {
        result = mul(result, a);
      }
    }
    return result;
  }

 private:
  //! (high 2^64 + low) / R mod p, for an input below p R
  std::uint64_t reduce(std::uint64_t high, std::uint64_t low) const {
    // m p cancels the low word exactly, so the low words carry only when low
    // is nonzero; the result is below 2p < 2^63.
    const std::uint64_t m = low * d_negInverse;
    std::uint64_t mpHigh, mpLow;
    multiplyWide(m, d_p, mpHigh, mpLow);
    const std::uint64_t r = high + mpHigh + (low != 0 ? 1 : 0);
    return r >= d_p ? r - d_p : r;
  }

  std::uint64_t d_p, d_negInverse, d_one, d_rSquared;
};

//! A symmetric elimination order for xI - A and the entries each step touches.
/*!
  Entries 0..n-1 are the diagonal; the others are the off-diagonal pairs the
  elimination reaches, original edges and fill-in alike, shared by both
  triangles.
*/
struct EliminationPlan {
  struct Step {
    size_t pivot;                   // diagonal entry of the eliminated vertex
    std::vector<size_t> toLater;    // entries (v, a), a still to come
    std::vector<size_t> laterDiag;  // diagonal entries of those a
    std::vector<size_t> pairs;      // entries (a, b), a before b, row-major
  };
  size_t numVertices = 0;
  size_t numEntries = 0;
  std::vector<size_t> edges;  // entries holding -1 initially
  std::vector<Step> steps;
};

//! Minimum-degree order: on a near-tree it eliminates leaves first and fills
//! in little.
EliminationPlan planElimination(
    unsigned int numVertices,
    const std::vector<std::pair<unsigned int, unsigned int>> &edges) {
  EliminationPlan plan;
  plan.numVertices = numVertices;
  plan.numEntries = numVertices;
  // Molecular degrees are small, so adjacency and entry lookups are short
  // linear scans. Each list holds (vertex, entry) pairs.
  std::vector<std::vector<std::pair<size_t, size_t>>> adjacency(numVertices);
  const auto connect = [&](size_t a, size_t b) {
    for (const auto &[vertex, entry] : adjacency[a]) {
      if (vertex == b) {
        return entry;
      }
    }
    const size_t entry = plan.numEntries++;
    adjacency[a].emplace_back(b, entry);
    adjacency[b].emplace_back(a, entry);
    return entry;
  };
  for (const auto &[a, b] : edges) {
    plan.edges.push_back(connect(a, b));
  }

  std::vector<bool> eliminated(numVertices, false);
  for (size_t step = 0; step < numVertices; ++step) {
    size_t v = numVertices;
    for (size_t u = 0; u < numVertices; ++u) {
      if (!eliminated[u] &&
          (v == numVertices || adjacency[u].size() < adjacency[v].size())) {
        v = u;
      }
    }
    EliminationPlan::Step current{v, {}, {}, {}};
    const auto later = adjacency[v];
    for (const auto &[vertex, entry] : later) {
      current.toLater.push_back(entry);
      current.laterDiag.push_back(vertex);
    }
    for (size_t i = 0; i < later.size(); ++i) {
      for (size_t j = i + 1; j < later.size(); ++j) {
        current.pairs.push_back(connect(later[i].first, later[j].first));
      }
    }
    for (const auto &[vertex, entry] : later) {
      auto &list = adjacency[vertex];
      list.erase(std::find_if(list.begin(), list.end(), [v](const auto &item) {
        return item.first == v;
      }));
    }
    adjacency[v].clear();
    eliminated[v] = true;
    plan.steps.push_back(std::move(current));
  }
  return plan;
}

//! det(xI - A) mod p at \c count consecutive points from \c first, by LDL^T
//! elimination over all points at once, so that each pivot vector shares one
//! inversion. A pivot that vanishes mod p makes its point unusable; it is
//! flagged rather than used.
void evaluateDeterminants(const EliminationPlan &plan,
                          const std::vector<bool> &loops,
                          const MontgomeryField &field, std::uint64_t first,
                          size_t count,
                          std::vector<std::uint64_t> &determinants,
                          std::vector<bool> &usable) {
  std::vector<std::uint64_t> values(plan.numEntries * count, 0);
  std::vector<std::uint64_t> inverse(count), prefix(count), multipliers;
  for (size_t v = 0; v < plan.numVertices; ++v) {
    for (size_t j = 0; j < count; ++j) {
      const std::uint64_t x = field.to(first + j);
      values[v * count + j] = loops[v] ? field.sub(x, field.one()) : x;
    }
  }
  const std::uint64_t minusOne = field.sub(0, field.one());
  for (const size_t edge : plan.edges) {
    std::fill_n(&values[edge * count], count, minusOne);
  }
  determinants.assign(count, field.one());
  usable.assign(count, true);

  for (const auto &step : plan.steps) {
    std::uint64_t *const pivot = &values[step.pivot * count];
    for (size_t j = 0; j < count; ++j) {
      if (pivot[j] == 0) {
        usable[j] = false;
        pivot[j] = field.one();  // keeps the batch inversion defined
      }
    }
    // Montgomery's trick: one inversion for the whole pivot vector
    prefix[0] = pivot[0];
    for (size_t j = 1; j < count; ++j) {
      prefix[j] = field.mul(prefix[j - 1], pivot[j]);
    }
    std::uint64_t running = field.inverse(prefix[count - 1]);
    for (size_t j = count; j-- > 1;) {
      inverse[j] = field.mul(running, prefix[j - 1]);
      running = field.mul(running, pivot[j]);
    }
    inverse[0] = running;
    for (size_t j = 0; j < count; ++j) {
      determinants[j] = field.mul(determinants[j], pivot[j]);
    }

    // Schur complement: M[a][b] -= M[a][v] M[v][b] / M[v][v] over the later
    // neighbours a, b of v.
    const size_t later = step.toLater.size();
    multipliers.assign(later * count, 0);
    for (size_t a = 0; a < later; ++a) {
      const std::uint64_t *const column = &values[step.toLater[a] * count];
      std::uint64_t *const factor = &multipliers[a * count];
      std::uint64_t *const diagonal = &values[step.laterDiag[a] * count];
      for (size_t j = 0; j < count; ++j) {
        factor[j] = field.mul(column[j], inverse[j]);
        diagonal[j] = field.sub(diagonal[j], field.mul(factor[j], column[j]));
      }
    }
    size_t pair = 0;
    for (size_t a = 0; a < later; ++a) {
      for (size_t b = a + 1; b < later; ++b, ++pair) {
        const std::uint64_t *const factor = &multipliers[a * count];
        const std::uint64_t *const column = &values[step.toLater[b] * count];
        std::uint64_t *const target = &values[step.pairs[pair] * count];
        for (size_t j = 0; j < count; ++j) {
          target[j] = field.sub(target[j], field.mul(factor[j], column[j]));
        }
      }
    }
  }
}

//! Coefficients a_0 .. a_n of det(xI - A) mod p, a_i multiplying x^i,
//! interpolated from values at points above the largest row sum of A.
std::vector<std::uint64_t> coefficientsModPrime(const EliminationPlan &plan,
                                                const std::vector<bool> &loops,
                                                const MontgomeryField &field,
                                                std::uint64_t firstPoint) {
  const size_t n = plan.numVertices;

  // n + 1 values fix a polynomial of degree n; evaluating in blocks bounds
  // the memory, and a few spare points absorb unusable ones.
  constexpr size_t maxBlock = 256;
  std::vector<std::uint64_t> xs, ys, determinants;
  std::vector<bool> usable;
  for (std::uint64_t next = firstPoint; xs.size() < n + 1;) {
    const size_t count = std::min<size_t>(maxBlock, n + 1 - xs.size() + 4);
    evaluateDeterminants(plan, loops, field, next, count, determinants, usable);
    for (size_t j = 0; j < count && xs.size() < n + 1; ++j) {
      if (usable[j]) {
        xs.push_back(next + j);
        ys.push_back(determinants[j]);
      }
    }
    next += count;
  }

  // Newton's divided differences. The points are increasing integers, so
  // every gap is a small positive integer, and one table holds the inverses:
  // inv(i) = -(p / i) inv(p mod i).
  const std::uint64_t p = field.prime();
  const size_t span = static_cast<size_t>(xs.back() - xs.front());
  std::vector<std::uint64_t> smallInverse(span + 1, 0);
  if (span >= 1) {
    smallInverse[1] = field.one();
  }
  for (size_t i = 2; i <= span; ++i) {
    smallInverse[i] =
        field.sub(0, field.mul(field.to(p / i), smallInverse[p % i]));
  }
  for (size_t k = 1; k <= n; ++k) {
    for (size_t i = n; i >= k; --i) {
      ys[i] = field.mul(field.sub(ys[i], ys[i - 1]),
                        smallInverse[xs[i] - xs[i - k]]);
    }
  }

  // From Newton's form: P = ys[n], then P = P (x - xs[i]) + ys[i] for
  // i = n-1 .. 0.
  std::vector<std::uint64_t> polynomial(n + 1, 0);
  polynomial[0] = ys[n];
  for (size_t i = n, degree = 0; i-- > 0; ++degree) {
    const std::uint64_t xi = field.to(xs[i]);
    for (size_t t = degree + 1; t > 0; --t) {
      polynomial[t] =
          field.sub(polynomial[t - 1], field.mul(xi, polynomial[t]));
    }
    polynomial[0] =
        field.add(field.sub(0, field.mul(xi, polynomial[0])), ys[i]);
  }
  for (auto &c : polynomial) {
    c = field.from(c);
  }
  return polynomial;
}

//! log2 of a bound on every |c_k|: |c_k| = |e_k(lambda)| <= e_k(|lambda|)
//! <= C(n,k) (sum |lambda| / n)^k by Maclaurin's inequality, and
//! sum |lambda| <= sqrt(n tr A^2) by Cauchy-Schwarz.
double coefficientBitsBound(size_t n, size_t traceOfSquare) {
  const double halfLogDensity =
      0.5 * std::log2(static_cast<double>(traceOfSquare) / n);
  const double log2NFactorial = std::lgamma(n + 1.0) / std::log(2.0);
  double bound = 0.0;
  for (size_t k = 0; k <= n; ++k) {
    const double log2Binomial =
        log2NFactorial -
        (std::lgamma(k + 1.0) + std::lgamma(n - k + 1.0)) / std::log(2.0);
    bound = std::max(bound, log2Binomial + k * halfLogDensity);
  }
  return bound;
}

//! |value| rounded to the nearest double, ties to even.
double magnitude(BigInt value) {
  if (value < 0) {
    value = -value;
  }
  if (value == 0) {
    return 0.0;
  }
  const unsigned int top = boost::multiprecision::msb(value);
  if (top < 64) {
    return static_cast<double>(static_cast<std::uint64_t>(value));
  }
  // Keep the leading 64 bits and fold every dropped bit into the lowest one:
  // that bit lies below double's rounding position, so the conversion still
  // rounds exactly.
  const unsigned int shift = top - 63;
  auto head = static_cast<std::uint64_t>(value >> shift);
  if (boost::multiprecision::lsb(value) < shift) {
    head |= 1u;
  }
  return std::ldexp(static_cast<double>(head), static_cast<int>(shift));
}

}  // namespace

std::optional<std::vector<double>> characteristicPolynomialMagnitudes(
    unsigned int numVertices,
    const std::vector<std::pair<unsigned int, unsigned int>> &edges,
    const std::vector<bool> &loops) {
  const size_t n = numVertices;
  std::vector<double> result(n + 1, 0.0);
  result[0] = 1.0;
  std::vector<size_t> rowSums(n, 0);
  size_t traceOfSquare = 0;  // tr A^2 = 2 |edges| + |loops|
  for (const auto &[a, b] : edges) {
    ++rowSums[a];
    ++rowSums[b];
    traceOfSquare += 2;
  }
  for (size_t v = 0; v < n; ++v) {
    if (loops[v]) {
      ++rowSums[v];
      ++traceOfSquare;
    }
  }
  if (traceOfSquare == 0) {
    return result;  // A = 0, so det(xI - A) = x^n
  }

  // Symmetric residues need a modulus above twice the largest |c_k|; two
  // bits of slack cover the floating-point evaluation of the bound.
  const double neededBits = coefficientBitsBound(n, traceOfSquare) + 3.0;
  size_t numPrimes = 0;
  for (double covered = 0.0; covered < neededBits; ++numPrimes) {
    if (numPrimes == largePrimes.size()) {
      return std::nullopt;
    }
    covered += std::log2(static_cast<double>(largePrimes[numPrimes]));
  }

  // Every eigenvalue of every principal submatrix of A lies within the
  // largest row sum, so above it no pivot vanishes over the integers; modulo
  // p one does only by rare accident.
  const EliminationPlan plan = planElimination(numVertices, edges);
  const std::uint64_t first =
      *std::max_element(rowSums.begin(), rowSums.end()) + 1;
  const auto residues =
      coefficientsModPrime(plan, loops, MontgomeryField(largePrimes[0]), first);

  if (numPrimes == 1) {
    const std::uint64_t p = largePrimes[0];
    for (size_t k = 0; k <= n; ++k) {
      const std::uint64_t r = residues[n - k];  // c_k multiplies x^(n-k)
      result[k] = static_cast<double>(r > p / 2 ? p - r : r);
    }
    return result;
  }

  // Chinese remainder theorem, one prime at a time:
  // x += M ((r - x) M^-1 mod p).
  std::vector<BigInt> exact(residues.begin(), residues.end());
  BigInt modulus = largePrimes[0];
  for (size_t i = 1; i < numPrimes; ++i) {
    const std::uint64_t p = largePrimes[i];
    const MontgomeryField field(p);
    const auto r = coefficientsModPrime(plan, loops, field, first);
    const std::uint64_t inverseOfModulus =
        field.inverse(field.to(static_cast<std::uint64_t>(modulus % p)));
    for (size_t k = 0; k <= n; ++k) {
      const auto current = static_cast<std::uint64_t>(exact[k] % p);
      const std::uint64_t delta =
          r[k] >= current ? r[k] - current : r[k] + p - current;
      exact[k] +=
          modulus * field.from(field.mul(field.to(delta), inverseOfModulus));
    }
    modulus *= p;
  }
  const BigInt half = modulus / 2;
  for (size_t k = 0; k <= n; ++k) {
    const BigInt &x = exact[n - k];
    result[k] = magnitude(x > half ? BigInt(x - modulus) : x);
  }
  return result;
}

}  // namespace detail
}  // namespace Descriptors
}  // namespace RDKit
