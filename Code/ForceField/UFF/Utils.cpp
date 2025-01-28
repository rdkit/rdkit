//
//  Copyright (C) 2013-2024 Paolo Tosco and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "Utils.h"
#include <cmath>

namespace ForceFields {
namespace UFF {
namespace Utils {

double calculateCosY(const RDGeom::Point3D &iPoint,
                     const RDGeom::Point3D &jPoint,
                     const RDGeom::Point3D &kPoint,
                     const RDGeom::Point3D &lPoint) {
  constexpr double zeroTol = 1.0e-16;
  RDGeom::Point3D rJI = iPoint - jPoint;
  RDGeom::Point3D rJK = kPoint - jPoint;
  RDGeom::Point3D rJL = lPoint - jPoint;
  auto l2JI = rJI.lengthSq();
  auto l2JK = rJK.lengthSq();
  auto l2JL = rJL.lengthSq();
  if (l2JI < zeroTol || l2JK < zeroTol || l2JL < zeroTol) {
    return 0.0;
  }
  RDGeom::Point3D n = rJI.crossProduct(rJK);
  n /= (sqrt(l2JI) * sqrt(l2JK));
  auto l2n = n.lengthSq();
  if (l2n < zeroTol) {
    return 0.0;
  }
  return n.dotProduct(rJL) / (sqrt(l2JL) * sqrt(l2n));
}

std::tuple<double, double, double, double>
calcInversionCoefficientsAndForceConstant(int at2AtomicNum, bool isCBoundToO) {
  double res = 0.0;
  double C0 = 0.0;
  double C1 = 0.0;
  double C2 = 0.0;
  // if the central atom is sp2 carbon, nitrogen or oxygen
  if ((at2AtomicNum == 6) || (at2AtomicNum == 7) || (at2AtomicNum == 8)) {
    C0 = 1.0;
    C1 = -1.0;
    C2 = 0.0;
    res = (isCBoundToO ? 50.0 : 6.0);
  } else {
    // group 5 elements are not clearly explained in the UFF paper
    // the following code was inspired by MCCCS Towhee's ffuff.F
    double w0 = M_PI / 180.0;
    switch (at2AtomicNum) {
      // if the central atom is phosphorous
      case 15:
        w0 *= 84.4339;
        break;

      // if the central atom is arsenic
      case 33:
        w0 *= 86.9735;
        break;

      // if the central atom is antimonium
      case 51:
        w0 *= 87.7047;
        break;

      // if the central atom is bismuth
      case 83:
        w0 *= 90.0;
        break;
    }
    C2 = 1.0;
    C1 = -4.0 * cos(w0);
    C0 = -(C1 * cos(w0) + C2 * cos(2.0 * w0));
    res = 22.0 / (C0 + C1 + C2);
  }
  res /= 3.0;

  return std::make_tuple(res, C0, C1, C2);
}

}  // end of namespace Utils
}  // namespace UFF
}  // namespace ForceFields