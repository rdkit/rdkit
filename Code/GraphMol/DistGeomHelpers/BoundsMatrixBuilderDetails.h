//
//  Copyright (C) 2026 Katharina Buchthal and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <vector>
#include "RDGeneral/Invariant.h"
#include <ranges>
#include <optional>
#include <algorithm>

#ifndef RD_BOUNDS_MATRIX_BUILDER_DETAILS_H
#define RD_BOUNDS_MATRIX_BUILDER_DETAILS_H

namespace RDKit {
namespace DGeomHelpers {

struct Bounds {
  double lower, upper;
  unsigned int aid1, aid4;

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
