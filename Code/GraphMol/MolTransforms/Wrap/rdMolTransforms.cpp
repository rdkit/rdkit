// $Id$
//
//  Copyright (C) 2005-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define PY_ARRAY_UNIQUE_SYMBOL rdmoltransforms_array_API
#include <RDBoost/python.h>
#include <RDBoost/import_array.h>
#include "numpy/arrayobject.h"
#include <GraphMol/ROMol.h>
#include <RDBoost/Wrap.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <Geometry/Transform3D.h>
#include <Geometry/point.h>

namespace python = boost::python;

namespace RDKit {
PyObject *computeCanonTrans(const Conformer &conf,
                            const RDGeom::Point3D *center = nullptr,
                            bool normalizeCovar = false, bool ignoreHs = true) {
  RDGeom::Transform3D *trans;
  trans = MolTransforms::computeCanonicalTransform(conf, center, normalizeCovar,
                                                   ignoreHs);
  npy_intp dims[2];
  dims[0] = 4;
  dims[1] = 4;
  auto *res = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
  auto *resData = reinterpret_cast<double *>(PyArray_DATA(res));
  const double *tdata = trans->getData();
  memcpy(static_cast<void *>(resData), static_cast<const void *>(tdata),
         4 * 4 * sizeof(double));
  delete trans;
  return PyArray_Return(res);
}

#ifdef RDK_HAS_EIGEN3
PyObject *computePrincAxesMomentsHelper(
    bool func(const Conformer &, Eigen::Matrix3d &, Eigen::Vector3d &, bool,
              bool, const std::vector<double> *),
    const Conformer &conf, bool ignoreHs, const python::object &weights) {
  Eigen::Matrix3d axes;
  Eigen::Vector3d moments;
  std::vector<double> *weightsVecPtr = nullptr;
  std::vector<double> weightsVec;
  size_t i;
  if (weights != python::object()) {
    size_t numElements = python::extract<int>(weights.attr("__len__")());
    if (numElements != conf.getNumAtoms()) {
      throw ValueErrorException(
          "The Python container must have length equal to conf.GetNumAtoms()");
    }
    weightsVec.resize(numElements);
    for (i = 0; i < numElements; ++i) {
      weightsVec[i] = python::extract<double>(weights[i]);
    }
    weightsVecPtr = &weightsVec;
  }
  PyObject *res = PyTuple_New(2);
  bool success = func(conf, axes, moments, ignoreHs, true, weightsVecPtr);
  if (success) {
    npy_intp dims[2];
    dims[0] = 3;
    dims[1] = 3;
    auto *axesNpy = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    auto *axesNpyData = reinterpret_cast<double *>(PyArray_DATA(axesNpy));
    i = 0;
    for (size_t y = 0; y < 3; ++y) {
      for (size_t x = 0; x < 3; ++x) {
        axesNpyData[i++] = axes(y, x);
      }
    }
    auto *momentsNpy = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    auto *momentsNpyData = reinterpret_cast<double *>(PyArray_DATA(momentsNpy));
    for (size_t y = 0; y < 3; ++y) {
      momentsNpyData[y] = moments(y);
    }
    res = PyTuple_New(2);
    PyTuple_SetItem(res, 0, (PyObject *)axesNpy);
    PyTuple_SetItem(res, 1, (PyObject *)momentsNpy);
  } else {
    PyTuple_SetItem(res, 0, Py_None);
    PyTuple_SetItem(res, 1, Py_None);
  }
  return res;
}

PyObject *computePrincAxesMoments(const Conformer &conf, bool ignoreHs,
                                  const python::object &weights) {
  return computePrincAxesMomentsHelper(
      MolTransforms::computePrincipalAxesAndMoments, conf, ignoreHs, weights);
}

PyObject *computePrincAxesMomentsFromGyrationMatrix(
    const Conformer &conf, bool ignoreHs, const python::object &weights) {
  return computePrincAxesMomentsHelper(
      MolTransforms::computePrincipalAxesAndMomentsFromGyrationMatrix, conf,
      ignoreHs, weights);
}
#endif

void transConformer(Conformer &conf, python::object trans) {
  PyObject *transObj = trans.ptr();
  if (!PyArray_Check(transObj)) {
    throw_value_error("Expecting a numeric array for transformation");
  }
  auto *transMat = reinterpret_cast<PyArrayObject *>(transObj);
  unsigned int nrows = PyArray_DIM(transMat, 0);
  unsigned int dSize = nrows * nrows;
  auto *inData = reinterpret_cast<double *>(PyArray_DATA(transMat));
  RDGeom::Transform3D transform;
  double *tData = transform.getData();
  memcpy(static_cast<void *>(tData), static_cast<void *>(inData),
         dSize * sizeof(double));
  MolTransforms::transformConformer(conf, transform);
}
}  // namespace RDKit

BOOST_PYTHON_MODULE(rdMolTransforms) {
  python::scope().attr("__doc__") =
      "Module containing functions to perform 3D operations like rotate and "
      "translate conformations";

  rdkit_import_array();

  std::string docString =
      "Compute the centroid of the conformation - hydrogens are ignored and no attention\n\
                           is paid to the difference in sizes of the heavy atoms; however,\n\
                           an optional vector of weights can be passed.\n";
  python::def("ComputeCentroid", MolTransforms::computeCentroid,
              (python::arg("conf"), python::arg("ignoreHs") = true,
               python::arg("weights") =
                   static_cast<const std::vector<double> *>(nullptr)),
              docString.c_str());

  docString =
      "Compute the transformation required aligna conformer so that\n\
               the principal axes align up with the x,y, z axes\n\
               The conformer itself is left unchanged\n\
  ARGUMENTS:\n\
    - conf : the conformer of interest\n\
    - center : optional center point to compute the principal axes around (defaults to the centroid)\n\
    - normalizeCovar : optionally normalize the covariance matrix by the number of atoms\n";
  python::def(
      "ComputeCanonicalTransform", RDKit::computeCanonTrans,
      (python::arg("conf"),
       python::arg("center") = static_cast<RDGeom::Point3D *>(nullptr),
       python::arg("normalizeCovar") = false, python::arg("ignoreHs") = true),
      docString.c_str());

#ifdef RDK_HAS_EIGEN3
  docString =
      "Compute principal axes and moments of inertia for a conformer\n\
       These values are calculated from the inertia tensor:\n\
       Iij = - sum_{s=1..N}(w_s * r_{si} * r_{sj}) i != j\n\
       Iii = sum_{s=1..N} sum_{j!=i} (w_s * r_{sj} * r_{sj})\n\
       where the coordinates are relative to the center of mass.\n\
\n\
  ARGUMENTS:\n\
    - conf : the conformer of interest\n\
    - ignoreHs : if True, ignore hydrogen atoms\n\
    - weights : if present, used to weight the atomic coordinates\n\n\
  Returns a (principal axes, principal moments) tuple\n";
  python::def("ComputePrincipalAxesAndMoments", RDKit::computePrincAxesMoments,
              (python::arg("conf"), python::arg("ignoreHs") = true,
               python::arg("weights") = python::object()),
              docString.c_str());

  docString =
      "Compute principal axes and moments from the gyration matrix of a conformer\n\
       These values are calculated from the gyration matrix/tensor:\n\
       Iij = sum_{s=1..N}(w_s * r_{si} * r_{sj}) i != j\n\
       Iii = sum_{s=1..N} sum_{t!=s}(w_s * r_{si} * r_{ti})\n\
       where the coordinates are relative to the center of mass.\n\
\n\
  ARGUMENTS:\n\
    - conf : the conformer of interest\n\
    - ignoreHs : if True, ignore hydrogen atoms\n\
    - weights : if present, used to weight the atomic coordinates\n\n\
  Returns a (principal axes, principal moments) tuple\n";
  python::def("ComputePrincipalAxesAndMomentsFromGyrationMatrix",
              RDKit::computePrincAxesMomentsFromGyrationMatrix,
              (python::arg("conf"), python::arg("ignoreHs") = true,
               python::arg("weights") = python::object()),
              docString.c_str());
#endif

  python::def("TransformConformer", RDKit::transConformer,
              "Transform the coordinates of a conformer");

  docString =
      "Canonicalize the orientation of a conformer so that its principal axes\n\
               around the specified center point coincide with the x, y, z axes\n\
  \n\
  ARGUMENTS:\n\
    - conf : conformer of interest \n\
    - center : optionally center point about which the principal axes are computed \n\
                          if not specified the centroid of the conformer will be used\n\
    - normalizeCovar : Optionally normalize the covariance matrix by the number of atoms\n";
  python::def(
      "CanonicalizeConformer", MolTransforms::canonicalizeConformer,
      (python::arg("conf"), python::arg("center") = (RDGeom::Point3D *)nullptr,
       python::arg("normalizeCovar") = false, python::arg("ignoreHs") = true),
      docString.c_str());

  python::def("CanonicalizeMol", MolTransforms::canonicalizeMol,
              (python::arg("mol"), python::arg("normalizeCovar") = false,
               python::arg("ignoreHs") = true),
              "Loop over the conformers in a molecule and canonicalize their "
              "orientation");

  python::def(
      "GetBondLength", &MolTransforms::getBondLength,
      (python::arg("conf"), python::arg("iAtomId"), python::arg("jAtomId")),
      "Returns the bond length in angstrom between atoms i, j\n");
  python::def("SetBondLength", &MolTransforms::setBondLength,
              (python::arg("conf"), python::arg("iAtomId"),
               python::arg("jAtomId"), python::arg("value")),
              "Sets the bond length in angstrom between atoms i, j; "
              "all atoms bonded to atom j are moved\n");
  python::def("GetAngleRad", &MolTransforms::getAngleRad,
              (python::arg("conf"), python::arg("iAtomId"),
               python::arg("jAtomId"), python::arg("kAtomId")),
              "Returns the angle in radians between atoms i, j, k\n");
  python::def("GetAngleDeg", &MolTransforms::getAngleDeg,
              (python::arg("conf"), python::arg("iAtomId"),
               python::arg("jAtomId"), python::arg("kAtomId")),
              "Returns the angle in degrees between atoms i, j, k\n");
  python::def(
      "SetAngleRad", &MolTransforms::setAngleRad,
      (python::arg("conf"), python::arg("iAtomId"), python::arg("jAtomId"),
       python::arg("kAtomId"), python::arg("value")),
      "Sets the angle in radians between atoms i, j, k; "
      "all atoms bonded to atom k are moved\n");
  python::def(
      "SetAngleDeg", &MolTransforms::setAngleDeg,
      (python::arg("conf"), python::arg("iAtomId"), python::arg("jAtomId"),
       python::arg("kAtomId"), python::arg("value")),
      "Sets the angle in degrees between atoms i, j, k; "
      "all atoms bonded to atom k are moved\n");
  python::def(
      "GetDihedralRad", &MolTransforms::getDihedralRad,
      (python::arg("conf"), python::arg("iAtomId"), python::arg("jAtomId"),
       python::arg("kAtomId"), python::arg("lAtomId")),
      "Returns the dihedral angle in radians between atoms i, j, k, l\n");
  python::def(
      "GetDihedralDeg", &MolTransforms::getDihedralDeg,
      (python::arg("conf"), python::arg("iAtomId"), python::arg("jAtomId"),
       python::arg("kAtomId"), python::arg("lAtomId")),
      "Returns the dihedral angle in degrees between atoms i, j, k, l\n");
  python::def(
      "SetDihedralRad", &MolTransforms::setDihedralRad,
      (python::arg("conf"), python::arg("iAtomId"), python::arg("jAtomId"),
       python::arg("kAtomId"), python::arg("lAtomId"), python::arg("value")),
      "Sets the dihedral angle in radians between atoms i, j, k, l; "
      "all atoms bonded to atom l are moved\n");
  python::def("SetDihedralDeg", &MolTransforms::setDihedralDeg,
              "Sets the dihedral angle in degrees between atoms i, j, k, l; "
              "all atoms bonded to atom l are moved\n");
}
