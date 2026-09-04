//
//  Copyright (C) 2004-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <Geometry/point.h>
#include <Geometry/Transform3D.h>
#include <Numerics/Vector.h>
#include <Numerics/Alignment/AlignPoints.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

void fillPointVec(nb::object seq, std::vector<RDGeom::Point3D> &pts) {
  if (!nb::isinstance<nb::sequence>(seq)) {
    throw nb::value_error("expected a non-empty sequence of points");
  }
  auto s = nb::cast<nb::sequence>(seq);
  if (nb::len(s) == 0) {
    throw nb::value_error("expected a non-empty sequence of points");
  }

  // Determine mode by trying to cast the first element to Point3D
  RDGeom::Point3D testPt;
  bool point3dMode = nb::try_cast<RDGeom::Point3D>(s[0], testPt);

  for (auto item : s) {
    if (point3dMode) {
      RDGeom::Point3D pt;
      if (!nb::try_cast<RDGeom::Point3D>(item, pt)) {
        throw nb::value_error("non-Point3D found in sequence of points");
      }
      pts.push_back(pt);
    } else {
      if (!nb::isinstance<nb::sequence>(item)) {
        throw nb::value_error("non-sequence found in sequence of points");
      }
      auto coords = nb::cast<nb::sequence>(item);
      if (nb::len(coords) != 3) {
        throw nb::value_error("each point must have 3 coordinates");
      }
      pts.emplace_back(nb::cast<double>(coords[0]), nb::cast<double>(coords[1]),
                       nb::cast<double>(coords[2]));
    }
  }
}

auto makeTransform4x4(const RDGeom::Transform3D &trans) {
  double *resData = new double[16];
  memcpy(resData, trans.getData(), 16 * sizeof(double));
  nb::capsule owner(resData, [](void *f) noexcept {
    delete[] reinterpret_cast<double *>(f);
  });
  return nb::ndarray<nb::numpy, double, nb::ndim<2>>(resData, {4, 4}, owner);
}

nb::tuple AlignPointPairs(nb::object refPoints, nb::object probePoints,
                          nb::object weights = nb::none(), bool reflect = false,
                          unsigned int maxIterations = 50) {
  std::vector<RDGeom::Point3D> refOwned, probeOwned;
  fillPointVec(refPoints, refOwned);
  fillPointVec(probePoints, probeOwned);

  if (refOwned.size() != probeOwned.size()) {
    throw nb::value_error("Mis-match in number of points");
  }

  RDGeom::Point3DConstPtrVect refPts, probePts;
  for (auto &p : refOwned) refPts.push_back(&p);
  for (auto &p : probeOwned) probePts.push_back(&p);

  std::unique_ptr<RDNumeric::DoubleVector> wtsVec;
  if (!weights.is_none()) {
    auto wseq = nb::cast<nb::sequence>(weights);
    auto nwts = nb::len(wseq);
    if (nwts != refOwned.size()) {
      throw nb::value_error(
          "Number of weights supplied do not match the number of points");
    }
    wtsVec = std::make_unique<RDNumeric::DoubleVector>(nwts);
    size_t i = 0;
    for (auto w : wseq) {
      wtsVec->setVal(i++, nb::cast<double>(w));
    }
  }

  RDGeom::Transform3D trans;
  double ssd = RDNumeric::Alignments::AlignPoints(
      refPts, probePts, trans, wtsVec.get(), reflect, maxIterations);

  return nb::make_tuple(ssd, makeTransform4x4(trans));
}

}  // namespace

NB_MODULE(rdAlignment, m) {
  m.doc() = "Module containing functions to align pairs of points in 3D";

  m.def(
      "GetAlignmentTransform", &AlignPointPairs, "refPoints"_a, "probePoints"_a,
      "weights"_a = nb::none(), "reflect"_a = false, "maxIterations"_a = 50,
      R"DOC(Compute the optimal alignment (minimum RMSD) between two set of points using the quaternion algorithm

ARGUMENTS:

  - refPoints : reference points specified as a sequence of 3-sequences or sequence of Point3Ds
  - probePoints : probe points to align to reference points - same format
    restrictions as reference points apply here
  - weights : optional sequence of weights to associate to each pair of points
  - reflect : reflect the probe points before attempting alignment
  - maxIterations : maximum number of iterations for the eigen solver

RETURNS:

  a 2-tuple:
    - SSD value for the alignment
    - the 4x4 transform matrix, as a list of lists
)DOC");
}
