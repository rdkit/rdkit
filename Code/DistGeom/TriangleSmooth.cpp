//
//  Copyright (C) 2004-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "BoundsMatrix.h"
#include "TriangleSmooth.h"

namespace DistGeom {
bool triangleSmoothBounds(BoundsMatPtr boundsMat, double tol) {
  return triangleSmoothBounds(boundsMat.get(), tol);
}
bool triangleSmoothBounds(BoundsMatrix *boundsMat, double tol) {
  auto npt = boundsMat->numRows();
  for (auto k = 0u; k < npt; k++) {
    for (auto i = 0u; i < npt - 1; i++) {
      if (i == k) {
        continue;
      }
      auto ii = i;
      auto ik = k;
      if (ii > ik) {
        std::swap(ii, ik);
      }

      const auto Uik = boundsMat->getValUnchecked(ii, ik);  // upper bound
      const auto Lik = boundsMat->getValUnchecked(ik, ii);  // lower bound
      for (auto j = i + 1; j < npt; j++) {
        if (j == k) {
          continue;
        }
        auto jj = j;
        auto jk = k;
        if (jj > jk) {
          std::swap(jj, jk);
        }
        const auto Ukj = boundsMat->getValUnchecked(jj, jk);  // upper bound
        const auto sumUikUkj = Uik + Ukj;
        if (boundsMat->getValUnchecked(i, j) > sumUikUkj) {
          // adjust the upper bound
          boundsMat->setValUnchecked(i, j, sumUikUkj);
        }

        const auto diffLikUjk = Lik - Ukj;
        const auto diffLjkUik = boundsMat->getValUnchecked(jk, jj) - Uik;
        if (boundsMat->getValUnchecked(j, i) < diffLikUjk) {
          // adjust the lower bound
          boundsMat->setValUnchecked(j, i, diffLikUjk);
        } else if (boundsMat->getValUnchecked(j, i) < diffLjkUik) {
          // adjust the lower bound
          boundsMat->setValUnchecked(j, i, diffLjkUik);
        }
        const auto lBound = boundsMat->getValUnchecked(j, i);
        const auto uBound = boundsMat->getValUnchecked(i, j);
        if (tol > 0. && (lBound - uBound) / lBound > 0. &&
            (lBound - uBound) / lBound < tol) {
          // adjust the upper bound
          boundsMat->setValUnchecked(i, j, lBound);
        } else if (lBound - uBound > 0.) {
          return false;
        }
      }
    }
  }
  return true;
}
}  // namespace DistGeom
