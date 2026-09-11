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
#include "DistGeom/ZMatrixUtils.h"
#include "RDGeneral/Invariant.h"
#include <vector>
#include <optional>
#include <vector>
#include <ranges>
#include <algorithm>
#include "RDGeneral/Invariant.h"

#ifndef RD_BOUNDS_MATRIX_BUILDER_DETAILS_H
#define RD_BOUNDS_MATRIX_BUILDER_DETAILS_H

namespace RDKit {
namespace DGeomHelpers {
enum class TorsionType {
  CIS = 0,
  TRANS,
  FLEXIBLE,
  CUSTOM,
  CISTRANS
};

enum class Type14 {
  IN_CHAIN,
  IN_RING,
  TWO_IN_SAME_RING,
  TWO_IN_DIFF_RING,
  SHARE_RING_BOND,
  MACROCYCLE_TWO_IN_SAME_RING,
  MACROCYCLE_ALL_IN_SAME_RING
};

struct TorsionValue {
  TorsionType type = TorsionType::FLEXIBLE;
  std::optional<double> value = {};
  std::optional<double> extraDist = {};
  bool isForced = false;
};

inline DistGeom::TorsionRange ringTorsion(const std::size_t rSize) {
  double torsion = M_PI;
  switch (rSize) {
    case 4u:
      [[fallthrough]];
    case 5u:
      torsion = M_PI * 45.0 / 180.0;
      break;
    case 6u:
      torsion = M_PI * 60.0 / 180.0;
      break;
    case 7u:
      torsion = M_PI * 90.0 / 180.0;
      break;
    case 8u:
      torsion = M_PI * 100.0 / 180.0;
      break;
  }
  return {-torsion, torsion};
}

//! A structure used to store planar 14 paths - cis/trans
struct Path14Configuration {
  unsigned int bid1, bid2, bid3;
  unsigned int aid1, aid2, aid3, aid4;
  TorsionValue value;
  Type14 type14;
  std::size_t rSize = 0;

  DistGeom::TorsionCandidates toTorsionRange() const {
    switch (value.type) {
      case TorsionType::CIS:
        return DistGeom::TorsionValues{0.0};
      case TorsionType::TRANS:
        return DistGeom::TorsionValues{M_PI};
      case TorsionType::FLEXIBLE:
        return ringTorsion(rSize);
      case TorsionType::CISTRANS:
        return DistGeom::TorsionValues{0.0, M_PI};
      case TorsionType::CUSTOM:
        return DistGeom::TorsionValues{*value.value};
      default:
        break;
    }
    return ringTorsion(0);
  }
};

using PATH14_VECT = std::vector<Path14Configuration>;

inline std::size_t getUnifiedId(const unsigned int id1, const unsigned int id2,
                                const unsigned int n) {
  // returns an id for (id1, id2) independent of order within range (0, 2*n - 1)
  // assuming id1 < n and id2 < n
  return id1 < id2 ? (static_cast<std::size_t>(id1) * n + id2)
                   : (static_cast<std::size_t>(id2) * n + id1);
}

inline std::size_t getUnifiedId(const unsigned int id1, const unsigned int id2,
                                const unsigned int id3, const unsigned int n) {
  // returns an id for (id1, id2, id3) independent of order of id1, id3 within
  // range (0, 3*(n) - 1) assuming id1 < n, id2 < n and id3 < n
  return id1 < id3 ? (static_cast<std::size_t>(id1) * n * n + id2 * n + id3)
                   : (static_cast<std::size_t>(id3) * n * n + id2 * n + id1);
}

template <unsigned int numBondIds>
auto unifiedIdToBondIds(std::size_t id, const unsigned int n) {
  // implements TODO algorithm
  std::array<unsigned int, numBondIds> bondIds;

  for (auto i : std::views::iota(0u, numBondIds) | std::views::reverse) {
    bondIds[i] = static_cast<unsigned int>(id % n);
    id /= n;
  }

  return bondIds;
};

struct Bounds {
  double lower{1.0}, upper{-1.0};  // we start invalid
  unsigned int aid1{0}, aid4{0};

  inline bool valid() const { return lower <= upper; }

  bool operator==(const Bounds &) const = default;

  friend std::ostream &operator<<(std::ostream &os, const Bounds &b) {
    return os << "Bounds{"
              << "lower=" << b.lower << ", upper=" << b.upper
              << ", aid1=" << b.aid1 << ", aid4=" << b.aid4 << '}';
  }
};

inline Bounds merge(std::vector<Bounds> bounds) {
  PRECONDITION(bounds.size(), "Cannot merge empty list of bounds");

  std::ranges::sort(bounds, {}, &Bounds::lower);

  Bounds current = bounds.front();
  double componentUpper = current.upper;
  std::optional<double> resultLower;

  // What we are doing here:
  // U {i'=intersection(i_j,..,i_k) | {i_j, ..., i_k}\subset(I) ^ i` !=
  // \emptyset ^ !\exists(i_l): intersection(i`, i_l) != \emptyset}
  // or in other words:
  // we aim to find the union of all intersections that are maximal in a sense
  // that adding another arbitrary bounds to it, would lead into an empty set

  // we solve this by traversing the sorted bounds in a sweep manner while
  // keeping track on the current/active non-empty intersection
  // (currentIntersection), the largest upperBound that was reached so far
  // (this is needed since the currentIntersection.upper can be smaller than
  // that, losing track of potenial overlaps/intersections).
  // To avoid storing all maximal non-overlapping intersections (only the
  // first and last one is relevant), we store the lower bound of the first
  // maximal intersection in resultLower

  for (const auto &_bound : bounds | std::views::drop(1)) {
    if (_bound.lower <= current.upper) {
      // Case 1: _bounds intersects with currentIntersection => add to current
      // intersection
      //  we know that bounds are sorted by lower bounds =>
      // _bound.lower is always greater/equal currentIntersection.lower
      current.lower = _bound.lower;
      current.upper = std::min(current.upper, _bound.upper);
    } else {
      // Case 2: _bound is not overlapping with the current intersection => we
      // know that currentIntersection is maximal

      if (!resultLower) {
        resultLower = current.lower;
      }

      current.lower = _bound.lower;
      current.upper =
          _bound.lower <= componentUpper
              ? std::min(componentUpper,
                         _bound.upper)  // there is this at least former
                                        // bounds that is overlapping and
                                        // needs to be considered
              : _bound.upper;
    }

    componentUpper = std::max(componentUpper, _bound.upper);
  }
  return Bounds{.lower = resultLower.value_or(current.lower),
                .upper = current.upper,
                .aid1 = bounds.front().aid1,
                .aid4 = bounds.front().aid4};
}

}  // namespace DGeomHelpers
}  // namespace RDKit
#endif
