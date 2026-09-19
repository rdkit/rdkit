//
//  Copyright (C) 2026 ETH Zurich
//  Created by: Katharina Buchthal
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RD_Z_MATRIX_UTILS_H
#define RD_Z_MATRIX_UTILS_H

#include <RDGeneral/Invariant.h>
#include <Numerics/SquareMatrix.h>
#include <variant>
#include <unordered_set>

namespace DistGeom {

// helper type for the visitor #4 (From
// https://en.cppreference.com/cpp/utility/variant/visit2)
template <class... Ts>
struct overloaded : Ts... {
  using Ts::operator()...;
};

namespace {
constexpr double TORSION_PRECISION =
    1e9;  // assuming torsions are between M_PI and -M_PI
inline int64_t quantize(double val) {
  val = std::fmod(val + M_PI, 2.0 * M_PI);
  if (val < 0.0) {
    val += 2.0 * M_PI;
  }
  val -= M_PI;

  return static_cast<int64_t>(val * TORSION_PRECISION);
}

inline double dequantize(int64_t val) {
  return static_cast<double>(val) / TORSION_PRECISION;
}

}  // namespace

struct TorsionRange {
  double lower, upper;

  TorsionRange(double lower, double upper) : lower(lower), upper(upper) {
    PRECONDITION(lower <= upper, "Invalid range");
    d_qLower = quantize(lower);
    d_qUpper = quantize(upper);
  }

  TorsionRange(int64_t qLower, int64_t qUpper)
      : d_qLower(qLower), d_qUpper(qUpper) {
    lower = dequantize(qLower);
    upper = dequantize(qUpper);
  }

  bool qContains(const int64_t value) const {
    if (d_qLower <= d_qUpper) {
      return value >= d_qLower && value <= d_qUpper;
    }
    // range across M_PI and -M_PI
    return value <= d_qLower && value >= d_qUpper;
  }

  double width() const {
    auto res = upper - lower;
    if (res < 0) {
      // range across M_PI - -M_PI
      return res + 2.0 * M_PI;
    }
    return res;
  }

  bool operator<(const TorsionRange &other) const {
    return this->width() < other.width();
  }

  double sample(RDKit::double_source_type &rng) const {
    return lower + (upper - lower) * (rng());
  }

 private:
  int64_t d_qLower, d_qUpper;
};

struct TorsionValues {
  using value_type = int64_t;
  using iterator = std::unordered_set<value_type>::iterator;
  using const_iterator = std::unordered_set<value_type>::const_iterator;

  TorsionValues() = default;

  TorsionValues(std::initializer_list<double> values) {
    content.reserve(values.size());
    // make sure values are within -M_PI and M_PI
    for (auto val : values) {
      content.emplace(quantize(val));
    }
  }

  bool empty() const { return content.empty(); }
  std::size_t size() const { return content.size(); }
  const_iterator begin() const { return content.begin(); }
  const_iterator end() const { return content.end(); }
  const_iterator cbegin() const { return content.cbegin(); }
  const_iterator cend() const { return content.cend(); }

  friend TorsionValues merge(TorsionValues lhs,  // copying on purpose!
                             TorsionValues rhs) {
    rhs.content.merge(lhs.content);  // lhs <- lhs n rhs; rhs <- union
    return lhs;
  }

  double sample(RDKit::double_source_type &rng) const {
    const auto idxT =
        static_cast<std::size_t>(rng() * static_cast<double>(this->size()));

    auto it = this->cbegin();
    std::advance(it, idxT);
    // in TorsionValues, the quantizied torsion (int_64)
    // is stored -> we need to convert it to double
    // [-M_PI, M_PI] again
    return dequantize(*it);
  }

  double firstValue() const { return dequantize(*content.begin()); }

  std::pair<iterator, bool> insert(double val) {
    return content.insert(quantize(val));
  }

  std::pair<iterator, bool> insertQ(int64_t val) { return content.insert(val); }

 private:
  std::unordered_set<int64_t> content;
};

using TorsionCandidates = std::variant<TorsionRange, TorsionValues>;

namespace {

inline TorsionCandidates merge(const TorsionRange &lhs,
                               const TorsionRange &rhs) {
  if (std::min(lhs.upper, rhs.upper) >= std::max(lhs.lower, rhs.lower)) {
    return TorsionRange{std::max(lhs.lower, rhs.lower),
                        std::min(lhs.upper, rhs.upper)};
  }
  return TorsionRange{
      std::min(lhs.lower, rhs.lower),
      std::max(lhs.upper, rhs.upper)};  // we can adapt the same logic
                                        // as for 1-4 ranges again but
                                        // it could get a bit tidious
}

inline TorsionCandidates merge(const TorsionRange &range,
                               const TorsionValues &values) {
  // we select only those fixed values that are within the range, if there is
  // non, we take the min/max from the range and the values
  TorsionValues mergedValues;

  for (const auto &val : values) {
    if (range.qContains(val)) {
      mergedValues.insertQ(val);
    }
  }

  if (!mergedValues.empty()) {
    return mergedValues;
  }

  const auto &[minIt, maxIt] = std::ranges::minmax_element(values);
  const TorsionRange minMaxValues{*minIt, *maxIt};

  return merge(range, minMaxValues);
}

inline TorsionCandidates merge(const TorsionValues &values,
                               const TorsionRange &range) {
  return merge(range, values);
}

inline bool less(const TorsionRange &lhs, const TorsionRange &rhs) {
  return lhs < rhs;
}

inline bool less(const TorsionRange &, const TorsionValues &) {
  return false;  // distinct values are more constraint than ranges
}

inline bool less(const TorsionValues &, const TorsionRange &) { return true; }

inline bool less(const TorsionValues &lhs, const TorsionValues &rhs) {
  return lhs.size() < rhs.size();
}

};  // namespace

inline TorsionCandidates merge(const TorsionCandidates &lhs,
                               const TorsionCandidates &rhs) {
  // for mixture of torsion range and fixed values

  return std::visit(
      [](const auto &l, const auto &r) -> TorsionCandidates {
        return merge(l, r);
      },
      lhs, rhs);
}

inline double sample(const TorsionCandidates &torsionCand,
                     RDKit::double_source_type &rng) {
  return std::visit(
      [&rng](const auto &torsionCand) -> double {
        return torsionCand.sample(rng);
      },
      torsionCand);
}

inline bool less(const TorsionCandidates &lhs, const TorsionCandidates &rhs) {
  // lhs < rhs => lhs is more constraint than rhs
  return std::visit(
      [](const auto &l, const auto &r) -> bool { return less(l, r); }, lhs,
      rhs);
}
}  // namespace DistGeom

#endif