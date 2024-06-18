// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>
#include <RDGeneral/Exceptions.h>
#include <GraphMol/GraphMol.h>
#include <ForceField/ForceField.h>
#include <ForceField/UFF/DistanceConstraint.h>
#include <ForceField/UFF/AngleConstraint.h>
#include <ForceField/UFF/TorsionConstraint.h>
#include <ForceField/UFF/PositionConstraint.h>
#include <ForceField/MMFF/DistanceConstraint.h>
#include <ForceField/MMFF/AngleConstraint.h>
#include <ForceField/MMFF/TorsionConstraint.h>
#include <ForceField/MMFF/PositionConstraint.h>
#include "PyForceField.h"

using namespace ForceFields;
namespace python = boost::python;

void ForceFieldAddDistanceConstraint(PyForceField *self, unsigned int idx1,
                                     unsigned int idx2, double minLen,
                                     double maxLen, double forceConstant) {
  UFF::DistanceConstraintContrib *constraint;
  constraint = new UFF::DistanceConstraintContrib(
      self->field.get(), idx1, idx2, minLen, maxLen, forceConstant);
  self->field->contribs().push_back(ForceFields::ContribPtr(constraint));
}

void ForceFieldAddFixedPoint(PyForceField *self, unsigned int idx) {
  self->field->fixedPoints().push_back(idx);
}

void UFFAddDistanceConstraint(PyForceField *self, unsigned int idx1,
                              unsigned int idx2, bool relative, double minLen,
                              double maxLen, double forceConstant) {
  auto *constraint = new UFF::DistanceConstraintContrib(
      self->field.get(), idx1, idx2, relative, minLen, maxLen, forceConstant);
  self->field->contribs().push_back(ForceFields::ContribPtr(constraint));
}

void UFFAddAngleConstraint(PyForceField *self, unsigned int idx1,
                           unsigned int idx2, unsigned int idx3, bool relative,
                           double minAngleDeg, double maxAngleDeg,
                           double forceConstant) {
  auto *constraint = new UFF::AngleConstraintContrib(
      self->field.get(), idx1, idx2, idx3, relative, minAngleDeg, maxAngleDeg,
      forceConstant);
  self->field->contribs().push_back(ForceFields::ContribPtr(constraint));
}

void UFFAddTorsionConstraint(PyForceField *self, unsigned int idx1,
                             unsigned int idx2, unsigned int idx3,
                             unsigned int idx4, bool relative,
                             double minDihedralDeg, double maxDihedralDeg,
                             double forceConstant) {
  auto *constraint = new UFF::TorsionConstraintContrib(
      self->field.get(), idx1, idx2, idx3, idx4, relative, minDihedralDeg,
      maxDihedralDeg, forceConstant);
  self->field->contribs().push_back(ForceFields::ContribPtr(constraint));
}

void UFFAddPositionConstraint(PyForceField *self, unsigned int idx,
                              double maxDispl, double forceConstant) {
  auto *constraint = new UFF::PositionConstraintContrib(
      self->field.get(), idx, maxDispl, forceConstant);
  self->field->contribs().push_back(ForceFields::ContribPtr(constraint));
}

void MMFFAddDistanceConstraint(PyForceField *self, unsigned int idx1,
                               unsigned int idx2, bool relative, double minLen,
                               double maxLen, double forceConstant) {
  auto *constraint = new MMFF::DistanceConstraintContrib(
      self->field.get(), idx1, idx2, relative, minLen, maxLen, forceConstant);
  self->field->contribs().push_back(ForceFields::ContribPtr(constraint));
}

void MMFFAddAngleConstraint(PyForceField *self, unsigned int idx1,
                            unsigned int idx2, unsigned int idx3, bool relative,
                            double minAngleDeg, double maxAngleDeg,
                            double forceConstant) {
  auto *constraint = new MMFF::AngleConstraintContrib(
      self->field.get(), idx1, idx2, idx3, relative, minAngleDeg, maxAngleDeg,
      forceConstant);
  self->field->contribs().push_back(ForceFields::ContribPtr(constraint));
}

void MMFFAddTorsionConstraint(PyForceField *self, unsigned int idx1,
                              unsigned int idx2, unsigned int idx3,
                              unsigned int idx4, bool relative,
                              double minDihedralDeg, double maxDihedralDeg,
                              double forceConstant) {
  auto *constraint = new MMFF::TorsionConstraintContrib(
      self->field.get(), idx1, idx2, idx3, idx4, relative, minDihedralDeg,
      maxDihedralDeg, forceConstant);
  self->field->contribs().push_back(ForceFields::ContribPtr(constraint));
}

void MMFFAddPositionConstraint(PyForceField *self, unsigned int idx,
                               double maxDispl, double forceConstant) {
  auto *constraint = new MMFF::PositionConstraintContrib(
      self->field.get(), idx, maxDispl, forceConstant);
  self->field->contribs().push_back(ForceFields::ContribPtr(constraint));
}

PyObject *ForceFieldGetExtraPointLoc(PyForceField *self, unsigned int idx) {
  if (idx >= self->extraPoints.size()) {
    throw IndexErrorException(idx);
  }
  PyObject *res = PyTuple_New(3);
  PyTuple_SetItem(res, 0, PyFloat_FromDouble(self->extraPoints[idx]->x));
  PyTuple_SetItem(res, 1, PyFloat_FromDouble(self->extraPoints[idx]->y));
  PyTuple_SetItem(res, 2, PyFloat_FromDouble(self->extraPoints[idx]->z));
  return res;
}

double PyForceField::calcEnergyWithPos(const python::object &pos) {
  PRECONDITION(this->field, "no force field");
  if (pos != python::object()) {
    size_t s = this->field->dimension() * this->field->numPoints();
    size_t numElements = python::extract<size_t>(pos.attr("__len__")());
    if (s != numElements) {
      throw ValueErrorException(
          "The Python container must have length equal to Dimension() * "
          "NumPoints()");
    }
    std::vector<double> c(s);
    for (size_t i = 0; i < s; ++i) {
      c[i] = python::extract<double>(pos[i]);
    }
    return this->field->calcEnergy(c.data());
  } else {
    return this->field->calcEnergy();
  }
}

PyObject *PyForceField::positions() {
  PRECONDITION(this->field, "no force field");
  size_t s = this->field->dimension() * this->field->numPoints();
  PyObject *coordTuple = PyTuple_New(s);
  const RDGeom::PointPtrVect &p = this->field->positions();
  size_t i = 0;
  PyObject *coordItem;
  for (const auto pptr : p) {
    for (size_t j = 0; j < 3; ++j) {
      coordItem = PyFloat_FromDouble((*pptr)[j]);
      PyTuple_SetItem(coordTuple, i++, coordItem);
    }
  }
  return coordTuple;
}

PyObject *PyForceField::calcGradWithPos(const python::object &pos) {
  PRECONDITION(this->field, "no force field");
  size_t s = this->field->dimension() * this->field->numPoints();
  std::vector<double> g(s, 0.0);
  PyObject *gradTuple = PyTuple_New(s);
  if (pos != python::object()) {
    size_t numElements = python::extract<size_t>(pos.attr("__len__")());
    if (s != numElements) {
      throw ValueErrorException(
          "The Python container must have length equal to Dimension() * "
          "NumPoints()");
    }
    std::vector<double> c(s);
    for (size_t i = 0; i < s; ++i) {
      c[i] = python::extract<double>(pos[i]);
    }
    this->field->calcGrad(c.data(), g.data());
  } else {
    this->field->calcGrad(g.data());
  }
  for (size_t i = 0; i < s; ++i) {
    PyObject *coordItem = PyFloat_FromDouble(g[i]);
    PyTuple_SetItem(gradTuple, i, coordItem);
  }
  return gradTuple;
}

python::tuple PyForceField::minimizeTrajectory(unsigned int snapshotFreq,
                                               int maxIts, double forceTol,
                                               double energyTol) {
  PRECONDITION(this->field, "no force field");
  RDKit::SnapshotVect snapshotVect;
  int resInt = this->field->minimize(snapshotFreq, &snapshotVect, maxIts,
                                     forceTol, energyTol);
  python::list l;
  for (const auto &it : snapshotVect) {
    l.append(new RDKit::Snapshot(it));
  }
  return python::make_tuple(resInt, l);
}

PyObject *PyMMFFMolProperties::getMMFFBondStretchParams(
    const RDKit::ROMol &mol, const unsigned int idx1, const unsigned int idx2) {
  PyObject *res = nullptr;
  unsigned int bondType;
  ForceFields::MMFF::MMFFBond mmffBondStretchParams;
  if (mmffMolProperties->getMMFFBondStretchParams(mol, idx1, idx2, bondType,
                                                  mmffBondStretchParams)) {
    res = PyTuple_New(3);
    PyTuple_SetItem(res, 0, PyInt_FromLong(bondType));
    PyTuple_SetItem(res, 1, PyFloat_FromDouble(mmffBondStretchParams.kb));
    PyTuple_SetItem(res, 2, PyFloat_FromDouble(mmffBondStretchParams.r0));
  }
  return res;
};

PyObject *PyMMFFMolProperties::getMMFFAngleBendParams(const RDKit::ROMol &mol,
                                                      const unsigned int idx1,
                                                      const unsigned int idx2,
                                                      const unsigned int idx3) {
  PyObject *res = nullptr;
  unsigned int angleType;
  ForceFields::MMFF::MMFFAngle mmffAngleBendParams;
  if (mmffMolProperties->getMMFFAngleBendParams(
          mol, idx1, idx2, idx3, angleType, mmffAngleBendParams)) {
    res = PyTuple_New(3);
    PyTuple_SetItem(res, 0, PyInt_FromLong(angleType));
    PyTuple_SetItem(res, 1, PyFloat_FromDouble(mmffAngleBendParams.ka));
    PyTuple_SetItem(res, 2, PyFloat_FromDouble(mmffAngleBendParams.theta0));
  }
  return res;
};

PyObject *PyMMFFMolProperties::getMMFFStretchBendParams(
    const RDKit::ROMol &mol, const unsigned int idx1, const unsigned int idx2,
    const unsigned int idx3) {
  PyObject *res = nullptr;
  unsigned int stretchBendType;
  ForceFields::MMFF::MMFFStbn mmffStretchBendParams;
  ForceFields::MMFF::MMFFBond mmffBondStretchParams[2];
  ForceFields::MMFF::MMFFAngle mmffAngleBendParams;
  if (mmffMolProperties->getMMFFStretchBendParams(
          mol, idx1, idx2, idx3, stretchBendType, mmffStretchBendParams,
          mmffBondStretchParams, mmffAngleBendParams)) {
    res = PyTuple_New(3);
    PyTuple_SetItem(res, 0, PyInt_FromLong(stretchBendType));
    PyTuple_SetItem(res, 1, PyFloat_FromDouble(mmffStretchBendParams.kbaIJK));
    PyTuple_SetItem(res, 2, PyFloat_FromDouble(mmffStretchBendParams.kbaKJI));
  }
  return res;
};

PyObject *PyMMFFMolProperties::getMMFFTorsionParams(const RDKit::ROMol &mol,
                                                    const unsigned int idx1,
                                                    const unsigned int idx2,
                                                    const unsigned int idx3,
                                                    const unsigned int idx4) {
  PyObject *res = nullptr;
  unsigned int torType;
  ForceFields::MMFF::MMFFTor mmffTorsionParams;
  if (mmffMolProperties->getMMFFTorsionParams(mol, idx1, idx2, idx3, idx4,
                                              torType, mmffTorsionParams)) {
    res = PyTuple_New(4);
    PyTuple_SetItem(res, 0, PyInt_FromLong(torType));
    PyTuple_SetItem(res, 1, PyFloat_FromDouble(mmffTorsionParams.V1));
    PyTuple_SetItem(res, 2, PyFloat_FromDouble(mmffTorsionParams.V2));
    PyTuple_SetItem(res, 3, PyFloat_FromDouble(mmffTorsionParams.V3));
  }
  return res;
};

PyObject *PyMMFFMolProperties::getMMFFOopBendParams(const RDKit::ROMol &mol,
                                                    const unsigned int idx1,
                                                    const unsigned int idx2,
                                                    const unsigned int idx3,
                                                    const unsigned int idx4) {
  PyObject *res = nullptr;
  ForceFields::MMFF::MMFFOop mmffOopBendParams;
  if (mmffMolProperties->getMMFFOopBendParams(mol, idx1, idx2, idx3, idx4,
                                              mmffOopBendParams)) {
    res = PyFloat_FromDouble(mmffOopBendParams.koop);
  }
  return res;
};

PyObject *PyMMFFMolProperties::getMMFFVdWParams(const unsigned int idx1,
                                                const unsigned int idx2) {
  PyObject *res = nullptr;
  ForceFields::MMFF::MMFFVdWRijstarEps mmffVdWParams;
  if (mmffMolProperties->getMMFFVdWParams(idx1, idx2, mmffVdWParams)) {
    res = PyTuple_New(4);
    PyTuple_SetItem(res, 0,
                    PyFloat_FromDouble(mmffVdWParams.R_ij_starUnscaled));
    PyTuple_SetItem(res, 1, PyFloat_FromDouble(mmffVdWParams.epsilonUnscaled));
    PyTuple_SetItem(res, 2, PyFloat_FromDouble(mmffVdWParams.R_ij_star));
    PyTuple_SetItem(res, 3, PyFloat_FromDouble(mmffVdWParams.epsilon));
  }
  return res;
};

BOOST_PYTHON_MODULE(rdForceField) {
  python::scope().attr("__doc__") = "Exposes the ForceField class";

  std::string docString;

  python::class_<PyForceField>("ForceField", "A force field", python::no_init)
      .def("CalcEnergy",
           (double(PyForceField::*)(const python::object &) const) &
               PyForceField::calcEnergyWithPos,
           ((python::arg("self"), python::arg("pos") = python::object())),
           "Returns the energy (in kcal/mol) of the current arrangement\n"
           "or of the supplied coordinate list (if non-empty)")
      .def("CalcGrad", &PyForceField::calcGradWithPos,
           ((python::arg("self"), python::arg("pos") = python::object())),
           "Returns a tuple filled with the per-coordinate gradients\n"
           "of the current arrangement or of the supplied coordinate list "
           "(if non-empty)")
      .def("Positions", &PyForceField::positions, python::args("self"),
           "Returns a tuple filled with the coordinates of the\n"
           "points the ForceField is handling")
      .def("Dimension",
           (unsigned int (PyForceField::*)() const) & PyForceField::dimension,
           python::args("self"), "Returns the dimension of the ForceField")
      .def("NumPoints",
           (unsigned int (PyForceField::*)() const) & PyForceField::numPoints,
           python::args("self"),
           "Returns the number of points the ForceField is handling")
      .def("Minimize", &PyForceField::minimize,
           ((python::arg("self"), python::arg("maxIts") = 200),
            python::arg("forceTol") = 1e-4, python::arg("energyTol") = 1e-6),
           "Runs some minimization iterations.\n\n  Returns 0 if the "
           "minimization succeeded.")
      .def("MinimizeTrajectory", &PyForceField::minimizeTrajectory,
           ((python::arg("self"), python::arg("snapshotFreq")),
            python::arg("maxIts") = 200, python::arg("forceTol") = 1e-4,
            python::arg("energyTol") = 1e-6),
           "Runs some minimization iterations, recording the minimization "
           "trajectory every snapshotFreq steps.\n\n"
           "Returns a (int, []) tuple; the int is 0 if the minimization "
           "succeeded, "
           "while the list contains Snapshot objects.")
      .def("AddDistanceConstraint", ForceFieldAddDistanceConstraint,
           (python::arg("self"), python::arg("idx1"), python::arg("idx2"),
            python::arg("minLen"), python::arg("maxLen"),
            python::arg("forceConstant")),
           "Adds a distance constraint to the UFF force field "
           "(deprecated, use UFFAddDistanceConstraint instead).")
      .def("AddFixedPoint", ForceFieldAddFixedPoint,
           (python::arg("self"), python::arg("idx")),
           "Adds a fixed point to the force field.")
      .def("UFFAddDistanceConstraint", UFFAddDistanceConstraint,
           (python::arg("self"), python::arg("idx1"), python::arg("idx2"),
            python::arg("relative"), python::arg("minLen"),
            python::arg("maxLen"), python::arg("forceConstant")),
           "Adds a distance constraint to the UFF force field; if relative == "
           "True, "
           "then minLen and maxLen are intended as relative to the current "
           "distance.")
      .def("UFFAddAngleConstraint", UFFAddAngleConstraint,
           (python::arg("self"), python::arg("idx1"), python::arg("idx2"),
            python::arg("idx3"), python::arg("relative"),
            python::arg("minAngleDeg"), python::arg("maxAngleDeg"),
            python::arg("forceConstant")),
           "Adds an angle constraint to the UFF force field; if relative == "
           "True, "
           "then minAngleDeg and maxAngleDeg are intended as relative to the "
           "current angle.")
      .def("UFFAddTorsionConstraint", UFFAddTorsionConstraint,
           (python::arg("self"), python::arg("idx1"), python::arg("idx2"),
            python::arg("idx3"), python::arg("idx4"), python::arg("relative"),
            python::arg("minDihedralDeg"), python::arg("maxDihedralDeg"),
            python::arg("forceConstant")),
           "Adds a dihedral angle constraint to the UFF force field; if "
           "relative == True, "
           "then minDihedralDeg and maxDihedralDeg are intended as relative to "
           "the current "
           "dihedral angle.")
      .def("UFFAddPositionConstraint", UFFAddPositionConstraint,
           (python::arg("self"), python::arg("idx"), python::arg("maxDispl"),
            python::arg("forceConstant")),
           "Adds a position constraint to the UFF force field.")
      .def("MMFFAddDistanceConstraint", MMFFAddDistanceConstraint,
           (python::arg("self"), python::arg("idx1"), python::arg("idx2"),
            python::arg("relative"), python::arg("minLen"),
            python::arg("maxLen"), python::arg("forceConstant")),
           "Adds a distance constraint to the MMFF force field; if relative == "
           "True, "
           "then minLen and maxLen are intended as relative to the current "
           "distance.")
      .def("MMFFAddAngleConstraint", MMFFAddAngleConstraint,
           (python::arg("self"), python::arg("idx1"), python::arg("idx2"),
            python::arg("idx3"), python::arg("relative"),
            python::arg("minAngleDeg"), python::arg("maxAngleDeg"),
            python::arg("forceConstant")),
           "Adds an angle constraint to the MMFF force field; if relative == "
           "True, "
           "then minAngleDeg and maxAngleDeg are intended as relative to the "
           "current angle.")
      .def("MMFFAddTorsionConstraint", MMFFAddTorsionConstraint,
           (python::arg("self"), python::arg("idx1"), python::arg("idx2"),
            python::arg("idx3"), python::arg("idx4"), python::arg("relative"),
            python::arg("minDihedralDeg"), python::arg("maxDihedralDeg"),
            python::arg("forceConstant")),
           "Adds a dihedral angle constraint to the MMFF force field; if "
           "relative == True, "
           "then minDihedralDeg and maxDihedralDeg are intended as relative to "
           "the current "
           "dihedral angle.")
      .def("MMFFAddPositionConstraint", MMFFAddPositionConstraint,
           (python::arg("self"), python::arg("idx"), python::arg("maxDispl"),
            python::arg("forceConstant")),
           "Adds a position constraint to the MMFF force field.")
      .def("Initialize", &PyForceField::initialize, python::args("self"),
           "initializes the force field (call this before minimizing)")
      .def("AddExtraPoint", &PyForceField::addExtraPoint,
           (python::arg("self"), python::arg("x"), python::arg("y"),
            python::arg("z"), python::arg("fixed") = true),
           "Adds an extra point, this can be useful for adding constraints.")
      .def("GetExtraPointPos", ForceFieldGetExtraPointLoc,
           (python::arg("self"), python::arg("idx")),
           "returns the location of an extra point as a tuple");
  python::class_<PyMMFFMolProperties>(
      "MMFFMolProperties", "MMFF molecular properties", python::no_init)
      .def("GetMMFFAtomType", &PyMMFFMolProperties::getMMFFAtomType,
           (python::arg("self"), python::arg("idx")),
           "Retrieves MMFF atom type for atom with index idx")
      .def("GetMMFFFormalCharge", &PyMMFFMolProperties::getMMFFFormalCharge,
           (python::arg("self"), python::arg("idx")),
           "Retrieves MMFF formal charge for atom with index idx")
      .def("GetMMFFPartialCharge", &PyMMFFMolProperties::getMMFFPartialCharge,
           (python::arg("self"), python::arg("idx")),
           "Retrieves MMFF partial charge for atom with index idx")
      .def("GetMMFFBondStretchParams",
           &PyMMFFMolProperties::getMMFFBondStretchParams,
           (python::arg("self"), python::arg("mol"), python::arg("idx1"),
            python::arg("idx2")),
           "Retrieves MMFF bond stretch parameters for atoms with indexes "
           "idx1, idx2 "
           "as a (bondType, kb, r0) tuple, or None if no parameters could be "
           "found")
      .def("GetMMFFAngleBendParams",
           &PyMMFFMolProperties::getMMFFAngleBendParams,
           (python::arg("self"), python::arg("mol"), python::arg("idx1"),
            python::arg("idx2"), python::arg("idx3")),
           "Retrieves MMFF angle bend parameters for atoms with indexes idx1, "
           "idx2, idx3 "
           "as a (angleType, ka, theta0) tuple, or None if no parameters could "
           "be found")
      .def("GetMMFFStretchBendParams",
           &PyMMFFMolProperties::getMMFFStretchBendParams,
           (python::arg("self"), python::arg("mol"), python::arg("idx1"),
            python::arg("idx2"), python::arg("idx3")),
           "Retrieves MMFF stretch-bend parameters for atoms with indexes "
           "idx1, idx2, idx3 "
           "as a (stretchBendType, kbaIJK, kbaKJI) tuple, or None if no "
           "parameters could be found")
      .def("GetMMFFTorsionParams", &PyMMFFMolProperties::getMMFFTorsionParams,
           (python::arg("self"), python::arg("mol"), python::arg("idx1"),
            python::arg("idx2"), python::arg("idx3"), python::arg("idx4")),
           "Retrieves MMFF torsion parameters for atoms with indexes idx1, "
           "idx2, idx3, idx4 "
           "as a (torsionType, V1, V2, V3) tuple, or None if no parameters "
           "could be found")
      .def("GetMMFFOopBendParams", &PyMMFFMolProperties::getMMFFOopBendParams,
           (python::arg("self"), python::arg("mol"), python::arg("idx1"),
            python::arg("idx2"), python::arg("idx3"), python::arg("idx4")),
           "Retrieves MMFF out-of-plane bending force constant for atoms with "
           "indexes "
           "idx1, idx2, idx3, idx4 as a koop float value")
      .def("GetMMFFVdWParams", &PyMMFFMolProperties::getMMFFVdWParams,
           (python::arg("self"), python::arg("idx1"), python::arg("idx2")),
           "Retrieves MMFF van der Waals parameters for atoms with indexes "
           "idx1, idx2 as a (R_ij_starUnscaled, epsilonUnscaled, R_ij_star, "
           "epsilon) tuple, "
           "or None if no parameters could be found")
      .def("SetMMFFDielectricModel",
           &PyMMFFMolProperties::setMMFFDielectricModel,
           (python::arg("self"), python::arg("dielModel") = 1),
           "Sets the DielModel MMFF property (1: constant; 2: "
           "distance-dependent; "
           "defaults to constant)")
      .def("SetMMFFDielectricConstant",
           &PyMMFFMolProperties::setMMFFDielectricConstant,
           (python::arg("self"), python::arg("dielConst") = 1.0),
           "Sets the DielConst MMFF property (defaults to 1.0)")
      .def("SetMMFFBondTerm", &PyMMFFMolProperties::setMMFFBondTerm,
           (python::arg("self"), python::arg("state") = true),
           "Sets the bond term to be included in the MMFF equation (defaults "
           "to True)")
      .def("SetMMFFAngleTerm", &PyMMFFMolProperties::setMMFFAngleTerm,
           (python::arg("self"), python::arg("state") = true),
           "Sets the angle term to be included in the MMFF equation (defaults "
           "to True)")
      .def("SetMMFFStretchBendTerm",
           &PyMMFFMolProperties::setMMFFStretchBendTerm,
           (python::arg("self"), python::arg("state") = true),
           "Sets the stretch-bend term to be included in the MMFF equation "
           "(defaults to True)")
      .def("SetMMFFOopTerm", &PyMMFFMolProperties::setMMFFOopTerm,
           (python::arg("self"), python::arg("state") = true),
           "Sets the out-of-plane bend term to be included in the MMFF "
           "equation (defaults to True)")
      .def("SetMMFFTorsionTerm", &PyMMFFMolProperties::setMMFFTorsionTerm,
           (python::arg("self"), python::arg("state") = true),
           "Sets the torsional term to be included in the MMFF equation "
           "(defaults to True)")
      .def("SetMMFFVdWTerm", &PyMMFFMolProperties::setMMFFVdWTerm,
           (python::arg("self"), python::arg("state") = true),
           "Sets the Van der Waals term to be included in the MMFF equation "
           "(defaults to True)")
      .def("SetMMFFEleTerm", &PyMMFFMolProperties::setMMFFEleTerm,
           (python::arg("self"), python::arg("state") = true),
           "Sets the electrostatic term to be included in the MMFF equation "
           "(defaults to True)")
      .def("SetMMFFVariant", &PyMMFFMolProperties::setMMFFVariant,
           (python::arg("self"), python::arg("mmffVariant") = "MMFF94"),
           "Sets the MMFF variant to be used (\"MMFF94\" or \"MMFF94s\"; "
           "defaults to \"MMFF94\")")
      .def("SetMMFFVerbosity", &PyMMFFMolProperties::setMMFFVerbosity,
           (python::arg("self"), python::arg("verbosity") = 0),
           "Sets the MMFF verbosity (0: none; 1: low; 2: high; defaults to 0)");
}
/*
    (python::arg("self"), python::arg("mol"), python::arg("idx1"),
   python::arg("idx2")),
    "Retrieves MMFF bond stretch parameters for atoms with indexes idx1, idx2; "
    "as a tuple (bondType, kb, r0)")
*/
