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

  // Determine mode from the first element: Point3D objects (detected via
  // duck-typing) vs. coordinate sequences (lists, tuples, numpy rows, etc.)
  auto first = s[0];
  bool point3dMode =
      nb::hasattr(first, "x") && nb::hasattr(first, "y") && nb::hasattr(first, "z");

  for (auto item : s) {
    if (point3dMode) {
      if (!nb::hasattr(item, "x") || !nb::hasattr(item, "y") ||
          !nb::hasattr(item, "z")) {
        throw nb::value_error("non-Point3D found in sequence of points");
      }
      pts.emplace_back(nb::cast<double>(item.attr("x")),
                       nb::cast<double>(item.attr("y")),
                       nb::cast<double>(item.attr("z")));
    } else {
      if (!nb::isinstance<nb::sequence>(item)) {
        throw nb::value_error("non-sequence found in sequence of points");
      }
      auto s = nb::cast<nb::sequence>(item);
      if (nb::len(s) != 3) {
        throw nb::value_error("each point must have 3 coordinates");
      }
      pts.emplace_back(nb::cast<double>(s[0]), nb::cast<double>(s[1]),
                       nb::cast<double>(s[2]));
    }
  }
}

nb::tuple AlignPointPairs(nb::object refPoints, nb::object probePoints,
                          nb::object weights = nb::none(),
                          int reflect = 0,
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
  double ssd = RDNumeric::Alignments::AlignPoints(refPts, probePts, trans,
                                                   wtsVec.get(),
                                                   static_cast<bool>(reflect),
                                                   maxIterations);

  const double *tdata = trans.getData();
  nb::list mat;
  for (int i = 0; i < 4; i++) {
    nb::list row;
    for (int j = 0; j < 4; j++) {
      row.append(tdata[i * 4 + j]);
    }
    mat.append(row);
  }

  return nb::make_tuple(ssd, mat);
}

}  // namespace

NB_MODULE(rdAlignment, m) {
  m.doc() = "Module containing functions to align pairs of points in 3D";

  m.def(
      "GetAlignmentTransform", &AlignPointPairs, "refPoints"_a,
      "probePoints"_a, "weights"_a = nb::none(), "reflect"_a = 0,
      "maxIterations"_a = 50,
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
