//
//  Copyright (C) 2005-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/ndarray.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <Geometry/Transform3D.h>
#include <Geometry/point.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;
namespace {

RDGeom::Point3D computeCentroidHelper(const Conformer &conf, bool ignoreHs,
                                      nb::object weights) {
  std::vector<double> *weightsVecPtr = nullptr;
  std::vector<double> weightsVec;
  if (!weights.is_none()) {
    size_t numElements = nb::len(weights);
    weightsVec.resize(numElements);
    size_t i = 0;
    for (nb::handle h : weights) {
      weightsVec[i++] = nb::cast<double>(h);
    }
    weightsVecPtr = &weightsVec;
  }
  return MolTransforms::computeCentroid(conf, ignoreHs, weightsVecPtr);
}

auto computeCanonTrans(const Conformer &conf,
                       const RDGeom::Point3D *center = nullptr,
                       bool normalizeCovar = false, bool ignoreHs = true) {
  RDGeom::Transform3D *trans =
      MolTransforms::computeCanonicalTransform(conf, center, normalizeCovar,
                                               ignoreHs);
  double *resData = new double[4 * 4];
  const double *tdata = trans->getData();
  memcpy(static_cast<void *>(resData), static_cast<const void *>(tdata),
         4 * 4 * sizeof(double));
  delete trans;

  nb::capsule owner(resData, [](void *f) noexcept {
    double *data = reinterpret_cast<double *>(f);
    delete[] data;
  });
  return nb::ndarray<nb::numpy, double, nb::ndim<2>>(resData, {4, 4}, owner);
}

#ifdef RDK_HAS_EIGEN3
nb::object computePrincAxesMomentsHelper(
    bool func(const Conformer &, Eigen::Matrix3d &, Eigen::Vector3d &, bool,
              bool, const std::vector<double> *),
    const Conformer &conf, bool ignoreHs, nb::object weights) {
  Eigen::Matrix3d axes;
  Eigen::Vector3d moments;
  std::vector<double> *weightsVecPtr = nullptr;
  std::vector<double> weightsVec;
  if (!weights.is_none()) {
    size_t numElements = nb::len(weights);
    if (numElements != conf.getNumAtoms()) {
      throw ValueErrorException(
          "The Python container must have length equal to conf.GetNumAtoms()");
    }
    weightsVec.resize(numElements);
    size_t i = 0;
    for (nb::handle h : weights) {
      weightsVec[i++] = nb::cast<double>(h);
    }
    weightsVecPtr = &weightsVec;
  }
  bool success = func(conf, axes, moments, ignoreHs, true, weightsVecPtr);
  if (success) {
    double *axesData = new double[3 * 3];
    size_t i = 0;
    for (size_t y = 0; y < 3; ++y) {
      for (size_t x = 0; x < 3; ++x) {
        axesData[i++] = axes(y, x);
      }
    }
    nb::capsule axesOwner(axesData, [](void *f) noexcept {
      double *data = reinterpret_cast<double *>(f);
      delete[] data;
    });
    auto axesNb =
        nb::ndarray<nb::numpy, double, nb::ndim<2>>(axesData, {3, 3}, axesOwner);

    double *momentsData = new double[3];
    for (size_t y = 0; y < 3; ++y) {
      momentsData[y] = moments(y);
    }
    nb::capsule momentsOwner(momentsData, [](void *f) noexcept {
      double *data = reinterpret_cast<double *>(f);
      delete[] data;
    });
    auto momentsNb = nb::ndarray<nb::numpy, double, nb::ndim<1>>(
        momentsData, {3}, momentsOwner);

    return nb::make_tuple(axesNb, momentsNb);
  } else {
    return nb::make_tuple(nb::none(), nb::none());
  }
}

nb::object computePrincAxesMoments(const Conformer &conf, bool ignoreHs,
                                   nb::object weights) {
  return computePrincAxesMomentsHelper(
      MolTransforms::computePrincipalAxesAndMoments, conf, ignoreHs, weights);
}

nb::object computePrincAxesMomentsFromGyrationMatrix(const Conformer &conf,
                                                     bool ignoreHs,
                                                     nb::object weights) {
  return computePrincAxesMomentsHelper(
      MolTransforms::computePrincipalAxesAndMomentsFromGyrationMatrix, conf,
      ignoreHs, weights);
}
#endif

void transConformer(Conformer &conf,
                    const nb::ndarray<nb::numpy, const double,
                                      nb::shape<-1, -1>, nb::c_contig> &trans) {
  unsigned int nrows = trans.shape(0);
  unsigned int dSize = nrows * nrows;
  const double *inData = trans.data();
  RDGeom::Transform3D transform;
  double *tData = transform.getData();
  memcpy(static_cast<void *>(tData), static_cast<const void *>(inData),
         dSize * sizeof(double));
  MolTransforms::transformConformer(conf, transform);
}

}  // namespace

NB_MODULE(rdMolTransforms, m) {
  m.doc() = R"DOC(Module containing functions to perform 3D operations like rotate and
translate conformations)DOC";

  m.def(
      "ComputeCentroid", computeCentroidHelper, "conf"_a,
      "ignoreHs"_a = true, "weights"_a = nb::none(),
      R"DOC(Compute the centroid of the conformation - hydrogens are ignored and no attention
is paid to the difference in sizes of the heavy atoms; however,
an optional vector of weights can be passed.)DOC");

  m.def(
      "ComputeCanonicalTransform", computeCanonTrans, "conf"_a,
      "center"_a = static_cast<RDGeom::Point3D *>(nullptr),
      "normalizeCovar"_a = false, "ignoreHs"_a = true,
      R"DOC(Compute the transformation required to align a conformer so that
the principal axes align up with the x,y,z axes.
The conformer itself is left unchanged.
ARGUMENTS:
  - conf : the conformer of interest
  - center : optional center point to compute the principal axes around (defaults to the centroid)
  - normalizeCovar : optionally normalize the covariance matrix by the number of atoms)DOC");

#ifdef RDK_HAS_EIGEN3
  m.def(
      "ComputePrincipalAxesAndMoments", computePrincAxesMoments,
      "conf"_a, "ignoreHs"_a = true, "weights"_a = nb::none(),
      R"DOC(Compute principal axes and moments of inertia for a conformer.
These values are calculated from the inertia tensor:
  Iij = - sum_{s=1..N}(w_s * r_{si} * r_{sj}) i != j
  Iii = sum_{s=1..N} sum_{j!=i} (w_s * r_{sj} * r_{sj})
where the coordinates are relative to the center of mass.

ARGUMENTS:
  - conf : the conformer of interest
  - ignoreHs : if True, ignore hydrogen atoms
  - weights : if present, used to weight the atomic coordinates

Returns a (principal axes, principal moments) tuple)DOC");

  m.def(
      "ComputePrincipalAxesAndMomentsFromGyrationMatrix",
      computePrincAxesMomentsFromGyrationMatrix, "conf"_a,
      "ignoreHs"_a = true, "weights"_a = nb::none(),
      R"DOC(Compute principal axes and moments from the gyration matrix of a conformer.
These values are calculated from the gyration matrix/tensor:
  Iij = sum_{s=1..N}(w_s * r_{si} * r_{sj}) i != j
  Iii = sum_{s=1..N} sum_{t!=s}(w_s * r_{si} * r_{ti})
where the coordinates are relative to the center of mass.

ARGUMENTS:
  - conf : the conformer of interest
  - ignoreHs : if True, ignore hydrogen atoms
  - weights : if present, used to weight the atomic coordinates

Returns a (principal axes, principal moments) tuple)DOC");
#endif

  m.def("TransformConformer", transConformer, "conf"_a, "trans"_a,
        "Transform the coordinates of a conformer");

  m.def(
      "CanonicalizeConformer", MolTransforms::canonicalizeConformer, "conf"_a,
      "center"_a = static_cast<RDGeom::Point3D *>(nullptr),
      "normalizeCovar"_a = false, "ignoreHs"_a = true,
      R"DOC(Canonicalize the orientation of a conformer so that its principal axes
around the specified center point coincide with the x, y, z axes.

ARGUMENTS:
  - conf : conformer of interest
  - center : optionally center point about which the principal axes are computed;
if not specified the centroid of the conformer will be used
  - normalizeCovar : Optionally normalize the covariance matrix by the number of atoms)DOC");

  m.def("CanonicalizeMol", MolTransforms::canonicalizeMol, "mol"_a,
        "normalizeCovar"_a = false, "ignoreHs"_a = true,
        "Loop over the conformers in a molecule and canonicalize their "
        "orientation");

  m.def("GetBondLength", &MolTransforms::getBondLength, "conf"_a, "iAtomId"_a,
        "jAtomId"_a,
        "Returns the bond length in angstrom between atoms i, j\n");

  m.def("SetBondLength", &MolTransforms::setBondLength, "conf"_a, "iAtomId"_a,
        "jAtomId"_a, "value"_a,
        "Sets the bond length in angstrom between atoms i, j; "
        "all atoms bonded to atom j are moved\n");

  m.def("GetAngleRad", &MolTransforms::getAngleRad, "conf"_a, "iAtomId"_a,
        "jAtomId"_a, "kAtomId"_a,
        "Returns the angle in radians between atoms i, j, k\n");

  m.def("GetAngleDeg", &MolTransforms::getAngleDeg, "conf"_a, "iAtomId"_a,
        "jAtomId"_a, "kAtomId"_a,
        "Returns the angle in degrees between atoms i, j, k\n");

  m.def("SetAngleRad", &MolTransforms::setAngleRad, "conf"_a, "iAtomId"_a,
        "jAtomId"_a, "kAtomId"_a, "value"_a,
        "Sets the angle in radians between atoms i, j, k; "
        "all atoms bonded to atom k are moved\n");

  m.def("SetAngleDeg", &MolTransforms::setAngleDeg, "conf"_a, "iAtomId"_a,
        "jAtomId"_a, "kAtomId"_a, "value"_a,
        "Sets the angle in degrees between atoms i, j, k; "
        "all atoms bonded to atom k are moved\n");

  m.def("GetDihedralRad", &MolTransforms::getDihedralRad, "conf"_a,
        "iAtomId"_a, "jAtomId"_a, "kAtomId"_a, "lAtomId"_a,
        "Returns the dihedral angle in radians between atoms i, j, k, l\n");

  m.def("GetDihedralDeg", &MolTransforms::getDihedralDeg, "conf"_a,
        "iAtomId"_a, "jAtomId"_a, "kAtomId"_a, "lAtomId"_a,
        "Returns the dihedral angle in degrees between atoms i, j, k, l\n");

  m.def("SetDihedralRad", &MolTransforms::setDihedralRad, "conf"_a,
        "iAtomId"_a, "jAtomId"_a, "kAtomId"_a, "lAtomId"_a, "value"_a,
        "Sets the dihedral angle in radians between atoms i, j, k, l; "
        "all atoms bonded to atom l are moved\n");

  m.def("SetDihedralDeg", &MolTransforms::setDihedralDeg, "conf"_a,
        "iAtomId"_a, "jAtomId"_a, "kAtomId"_a, "lAtomId"_a, "value"_a,
        "Sets the dihedral angle in degrees between atoms i, j, k, l; "
        "all atoms bonded to atom l are moved\n");
}
