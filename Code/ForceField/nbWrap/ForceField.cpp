//
//  Copyright (C) 2004-2026 Rational Discovery LLC and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <RDBoost/Wrap_nb.h>
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
#include <GraphMol/Trajectory/Snapshot.h>
#include "PyForceField.h"

using namespace ForceFields;
namespace nb = nanobind;
using namespace nb::literals;

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

nb::tuple ForceFieldGetExtraPointLoc(PyForceField *self, unsigned int idx) {
  if (idx >= self->extraPoints.size()) {
    throw IndexErrorException(idx);
  }
  return nb::make_tuple(self->extraPoints[idx]->x, self->extraPoints[idx]->y,
                        self->extraPoints[idx]->z);
}

double PyForceField::calcEnergyWithPos(nb::object pos) {
  PRECONDITION(this->field, "no force field");
  if (!pos.is_none()) {
    size_t s = this->field->dimension() * this->field->numPoints();
    size_t numElements = nb::len(pos);
    if (s != numElements) {
      throw ValueErrorException(
          "The Python container must have length equal to Dimension() * "
          "NumPoints()");
    }
    std::vector<double> c(s);
    for (size_t i = 0; i < s; ++i) {
      c[i] = nb::cast<double>(pos[nb::cast(i)]);
    }
    return this->field->calcEnergy(c.data());
  } else {
    return this->field->calcEnergy();
  }
}

nb::tuple PyForceField::positions() {
  PRECONDITION(this->field, "no force field");
  const RDGeom::PointPtrVect &p = this->field->positions();
  nb::list coordList;
  for (const auto pptr : p) {
    for (size_t j = 0; j < 3; ++j) {
      coordList.append((*pptr)[j]);
    }
  }
  return nb::tuple(coordList);
}

nb::tuple PyForceField::calcGradWithPos(nb::object pos) {
  PRECONDITION(this->field, "no force field");
  size_t s = this->field->dimension() * this->field->numPoints();
  std::vector<double> g(s, 0.0);
  if (!pos.is_none()) {
    size_t numElements = nb::len(pos);
    if (s != numElements) {
      throw ValueErrorException(
          "The Python container must have length equal to Dimension() * "
          "NumPoints()");
    }
    std::vector<double> c(s);
    for (size_t i = 0; i < s; ++i) {
      c[i] = nb::cast<double>(pos[nb::cast(i)]);
    }
    this->field->calcGrad(c.data(), g.data());
  } else {
    this->field->calcGrad(g.data());
  }
  nb::list gradList;
  for (size_t i = 0; i < s; ++i) {
    gradList.append(g[i]);
  }
  return nb::tuple(gradList);
}

nb::tuple PyForceField::minimizeTrajectory(unsigned int snapshotFreq,
                                           int maxIts, double forceTol,
                                           double energyTol) {
  PRECONDITION(this->field, "no force field");
  RDKit::SnapshotVect snapshotVect;
  int resInt = this->field->minimize(snapshotFreq, &snapshotVect, maxIts,
                                     forceTol, energyTol);
  nb::list l;
  for (const auto &it : snapshotVect) {
    l.append(nb::cast(new RDKit::Snapshot(it), nb::rv_policy::take_ownership));
  }
  return nb::make_tuple(resInt, l);
}

nb::object PyMMFFMolProperties::getMMFFBondStretchParams(
    const RDKit::ROMol &mol, const unsigned int idx1,
    const unsigned int idx2) const {
  unsigned int bondType;
  ForceFields::MMFF::MMFFBond mmffBondStretchParams;
  if (mmffMolProperties->getMMFFBondStretchParams(mol, idx1, idx2, bondType,
                                                  mmffBondStretchParams)) {
    return nb::cast(nb::make_tuple((int)bondType, mmffBondStretchParams.kb,
                                   mmffBondStretchParams.r0));
  }
  return nb::none();
}

nb::object PyMMFFMolProperties::getMMFFAngleBendParams(
    const RDKit::ROMol &mol, const unsigned int idx1, const unsigned int idx2,
    const unsigned int idx3) const {
  unsigned int angleType;
  ForceFields::MMFF::MMFFAngle mmffAngleBendParams;
  if (mmffMolProperties->getMMFFAngleBendParams(
          mol, idx1, idx2, idx3, angleType, mmffAngleBendParams)) {
    return nb::cast(nb::make_tuple((int)angleType, mmffAngleBendParams.ka,
                                   mmffAngleBendParams.theta0));
  }
  return nb::none();
}

nb::object PyMMFFMolProperties::getMMFFStretchBendParams(
    const RDKit::ROMol &mol, const unsigned int idx1, const unsigned int idx2,
    const unsigned int idx3) const {
  unsigned int stretchBendType;
  ForceFields::MMFF::MMFFStbn mmffStretchBendParams;
  ForceFields::MMFF::MMFFBond mmffBondStretchParams[2];
  ForceFields::MMFF::MMFFAngle mmffAngleBendParams;
  if (mmffMolProperties->getMMFFStretchBendParams(
          mol, idx1, idx2, idx3, stretchBendType, mmffStretchBendParams,
          mmffBondStretchParams, mmffAngleBendParams)) {
    return nb::cast(nb::make_tuple((int)stretchBendType,
                                   mmffStretchBendParams.kbaIJK,
                                   mmffStretchBendParams.kbaKJI));
  }
  return nb::none();
}

nb::object PyMMFFMolProperties::getMMFFTorsionParams(
    const RDKit::ROMol &mol, const unsigned int idx1, const unsigned int idx2,
    const unsigned int idx3, const unsigned int idx4) const {
  unsigned int torType;
  ForceFields::MMFF::MMFFTor mmffTorsionParams;
  if (mmffMolProperties->getMMFFTorsionParams(mol, idx1, idx2, idx3, idx4,
                                              torType, mmffTorsionParams)) {
    return nb::cast(nb::make_tuple((int)torType, mmffTorsionParams.V1,
                                   mmffTorsionParams.V2, mmffTorsionParams.V3));
  }
  return nb::none();
}

nb::object PyMMFFMolProperties::getMMFFOopBendParams(
    const RDKit::ROMol &mol, const unsigned int idx1, const unsigned int idx2,
    const unsigned int idx3, const unsigned int idx4) const {
  ForceFields::MMFF::MMFFOop mmffOopBendParams;
  if (mmffMolProperties->getMMFFOopBendParams(mol, idx1, idx2, idx3, idx4,
                                              mmffOopBendParams)) {
    return nb::cast(mmffOopBendParams.koop);
  }
  return nb::none();
}

nb::object PyMMFFMolProperties::getMMFFVdWParams(const unsigned int idx1,
                                                 const unsigned int idx2) const {
  ForceFields::MMFF::MMFFVdWRijstarEps mmffVdWParams;
  if (mmffMolProperties->getMMFFVdWParams(idx1, idx2, mmffVdWParams)) {
    return nb::cast(
        nb::make_tuple(mmffVdWParams.R_ij_starUnscaled,
                       mmffVdWParams.epsilonUnscaled, mmffVdWParams.R_ij_star,
                       mmffVdWParams.epsilon));
  }
  return nb::none();
}

NB_MODULE(rdForceField, m) {
  m.doc() = "Exposes the ForceField class";

  // Minimal Snapshot binding needed for MinimizeTrajectory return value.
  // Full Snapshot bindings live in the Trajectory nbWrap (not yet migrated).
  nb::class_<RDKit::Snapshot>(m, "Snapshot",
                               "A snapshot of atomic coordinates from a "
                               "minimization trajectory")
      .def("GetPoint2D", &RDKit::Snapshot::getPoint2D, "pointNum"_a,
           "Returns the coordinates at pointNum as a Point2D object; "
           "requires the Trajectory dimension to be == 2")
      .def("GetPoint3D", &RDKit::Snapshot::getPoint3D, "pointNum"_a,
           "Returns the coordinates at pointNum as a Point3D object; "
           "requires the Trajectory dimension to be >= 2")
      .def("GetEnergy", &RDKit::Snapshot::getEnergy,
           "Returns the energy for this Snapshot")
      .def("SetEnergy", &RDKit::Snapshot::setEnergy, "energy"_a,
           "Sets the energy for this Snapshot");

  nb::class_<PyForceField>(m, "ForceField", "A force field")
      .def("CalcEnergy", &PyForceField::calcEnergyWithPos,
           "pos"_a = nb::none(),
           R"DOC(Returns the energy (in kcal/mol) of the current arrangement
or of the supplied coordinate list (if non-empty))DOC")
      .def("CalcGrad", &PyForceField::calcGradWithPos, "pos"_a = nb::none(),
           R"DOC(Returns a tuple filled with the per-coordinate gradients
of the current arrangement or of the supplied coordinate list (if non-empty))DOC")
      .def("Positions", &PyForceField::positions,
           R"DOC(Returns a tuple filled with the coordinates of the
points the ForceField is handling)DOC")
      .def("Dimension", &PyForceField::dimension,
           "Returns the dimension of the ForceField")
      .def("NumPoints", &PyForceField::numPoints,
           "Returns the number of points the ForceField is handling")
      .def("Minimize", &PyForceField::minimize, "maxIts"_a = 200,
           "forceTol"_a = 1e-4, "energyTol"_a = 1e-6,
           "Runs some minimization iterations.\n\n  Returns 0 if the "
           "minimization succeeded.")
      .def("MinimizeTrajectory", &PyForceField::minimizeTrajectory,
           "snapshotFreq"_a, "maxIts"_a = 200, "forceTol"_a = 1e-4,
           "energyTol"_a = 1e-6,
           R"DOC(Runs some minimization iterations, recording the minimization
trajectory every snapshotFreq steps.

Returns a (int, []) tuple; the int is 0 if the minimization succeeded,
while the list contains Snapshot objects.)DOC")
      .def("AddDistanceConstraint", ForceFieldAddDistanceConstraint,
           "idx1"_a, "idx2"_a, "minLen"_a, "maxLen"_a, "forceConstant"_a,
           "Adds a distance constraint to the UFF force field "
           "(deprecated, use UFFAddDistanceConstraint instead).")
      .def("AddFixedPoint", ForceFieldAddFixedPoint, "idx"_a,
           "Adds a fixed point to the force field.")
      .def("UFFAddDistanceConstraint", UFFAddDistanceConstraint,
           "idx1"_a, "idx2"_a, "relative"_a, "minLen"_a, "maxLen"_a,
           "forceConstant"_a,
           "Adds a distance constraint to the UFF force field; if relative == "
           "True, then minLen and maxLen are intended as relative to the "
           "current distance.")
      .def("UFFAddAngleConstraint", UFFAddAngleConstraint, "idx1"_a,
           "idx2"_a, "idx3"_a, "relative"_a, "minAngleDeg"_a, "maxAngleDeg"_a,
           "forceConstant"_a,
           "Adds an angle constraint to the UFF force field; if relative == "
           "True, then minAngleDeg and maxAngleDeg are intended as relative to "
           "the current angle.")
      .def("UFFAddTorsionConstraint", UFFAddTorsionConstraint,
           "idx1"_a, "idx2"_a, "idx3"_a, "idx4"_a, "relative"_a,
           "minDihedralDeg"_a, "maxDihedralDeg"_a, "forceConstant"_a,
           "Adds a dihedral angle constraint to the UFF force field; if "
           "relative == True, then minDihedralDeg and maxDihedralDeg are "
           "intended as relative to the current dihedral angle.")
      .def("UFFAddPositionConstraint", UFFAddPositionConstraint,
           "idx"_a, "maxDispl"_a, "forceConstant"_a,
           "Adds a position constraint to the UFF force field.")
      .def("MMFFAddDistanceConstraint", MMFFAddDistanceConstraint,
           "idx1"_a, "idx2"_a, "relative"_a, "minLen"_a, "maxLen"_a,
           "forceConstant"_a,
           "Adds a distance constraint to the MMFF force field; if relative == "
           "True, then minLen and maxLen are intended as relative to the "
           "current distance.")
      .def("MMFFAddAngleConstraint", MMFFAddAngleConstraint, "idx1"_a,
           "idx2"_a, "idx3"_a, "relative"_a, "minAngleDeg"_a, "maxAngleDeg"_a,
           "forceConstant"_a,
           "Adds an angle constraint to the MMFF force field; if relative == "
           "True, then minAngleDeg and maxAngleDeg are intended as relative to "
           "the current angle.")
      .def("MMFFAddTorsionConstraint", MMFFAddTorsionConstraint,
           "idx1"_a, "idx2"_a, "idx3"_a, "idx4"_a, "relative"_a,
           "minDihedralDeg"_a, "maxDihedralDeg"_a, "forceConstant"_a,
           "Adds a dihedral angle constraint to the MMFF force field; if "
           "relative == True, then minDihedralDeg and maxDihedralDeg are "
           "intended as relative to the current dihedral angle.")
      .def("MMFFAddPositionConstraint", MMFFAddPositionConstraint,
           "idx"_a, "maxDispl"_a, "forceConstant"_a,
           "Adds a position constraint to the MMFF force field.")
      .def("Initialize", &PyForceField::initialize,
           "initializes the force field (call this before minimizing)")
      .def("AddExtraPoint", &PyForceField::addExtraPoint, "x"_a, "y"_a, "z"_a,
           "fixed"_a = true,
           "Adds an extra point, this can be useful for adding constraints.")
      .def("GetExtraPointPos", ForceFieldGetExtraPointLoc, "idx"_a,
           "returns the location of an extra point as a tuple");

  nb::class_<PyMMFFMolProperties>(m, "MMFFMolProperties",
                                  "MMFF molecular properties")
      .def("GetMMFFAtomType", &PyMMFFMolProperties::getMMFFAtomType, "idx"_a,
           "Retrieves MMFF atom type for atom with index idx")
      .def("GetMMFFFormalCharge", &PyMMFFMolProperties::getMMFFFormalCharge,
           "idx"_a, "Retrieves MMFF formal charge for atom with index idx")
      .def("GetMMFFPartialCharge", &PyMMFFMolProperties::getMMFFPartialCharge,
           "idx"_a, "Retrieves MMFF partial charge for atom with index idx")
      .def("GetMMFFBondStretchParams",
           &PyMMFFMolProperties::getMMFFBondStretchParams, "mol"_a, "idx1"_a,
           "idx2"_a,
           "Retrieves MMFF bond stretch parameters for atoms with indexes "
           "idx1, idx2 as a (bondType, kb, r0) tuple, or None if no "
           "parameters could be found")
      .def("GetMMFFAngleBendParams",
           &PyMMFFMolProperties::getMMFFAngleBendParams, "mol"_a, "idx1"_a,
           "idx2"_a, "idx3"_a,
           "Retrieves MMFF angle bend parameters for atoms with indexes idx1, "
           "idx2, idx3 as a (angleType, ka, theta0) tuple, or None if no "
           "parameters could be found")
      .def("GetMMFFStretchBendParams",
           &PyMMFFMolProperties::getMMFFStretchBendParams, "mol"_a, "idx1"_a,
           "idx2"_a, "idx3"_a,
           "Retrieves MMFF stretch-bend parameters for atoms with indexes "
           "idx1, idx2, idx3 as a (stretchBendType, kbaIJK, kbaKJI) tuple, "
           "or None if no parameters could be found")
      .def("GetMMFFTorsionParams", &PyMMFFMolProperties::getMMFFTorsionParams,
           "mol"_a, "idx1"_a, "idx2"_a, "idx3"_a, "idx4"_a,
           "Retrieves MMFF torsion parameters for atoms with indexes idx1, "
           "idx2, idx3, idx4 as a (torsionType, V1, V2, V3) tuple, or None "
           "if no parameters could be found")
      .def("GetMMFFOopBendParams", &PyMMFFMolProperties::getMMFFOopBendParams,
           "mol"_a, "idx1"_a, "idx2"_a, "idx3"_a, "idx4"_a,
           "Retrieves MMFF out-of-plane bending force constant for atoms with "
           "indexes idx1, idx2, idx3, idx4 as a koop float value")
      .def("GetMMFFVdWParams", &PyMMFFMolProperties::getMMFFVdWParams, "idx1"_a,
           "idx2"_a,
           "Retrieves MMFF van der Waals parameters for atoms with indexes "
           "idx1, idx2 as a (R_ij_starUnscaled, epsilonUnscaled, R_ij_star, "
           "epsilon) tuple, or None if no parameters could be found")
      .def("SetMMFFDielectricModel",
           &PyMMFFMolProperties::setMMFFDielectricModel, "dielModel"_a = 1,
           "Sets the DielModel MMFF property (1: constant; 2: "
           "distance-dependent; defaults to constant)")
      .def("GetMMFFDielectricModel",
           &PyMMFFMolProperties::getMMFFDielectricModel,
           "Returns the currently configured MMFF dielectric model "
           "(1: constant; 2: distance-dependent).")
      .def("SetMMFFDielectricConstant",
           &PyMMFFMolProperties::setMMFFDielectricConstant, "dielConst"_a = 1.0,
           "Sets the DielConst MMFF property (defaults to 1.0)")
      .def("GetMMFFDielectricConstant",
           &PyMMFFMolProperties::getMMFFDielectricConstant,
           "Returns the currently configured MMFF dielectric constant.")
      .def("SetMMFFBondTerm", &PyMMFFMolProperties::setMMFFBondTerm,
           "state"_a = true,
           "Sets the bond term to be included in the MMFF equation "
           "(defaults to True)")
      .def("GetMMFFBondTerm", &PyMMFFMolProperties::getMMFFBondTerm,
           "Returns whether the bond term is included in the MMFF equation.")
      .def("SetMMFFAngleTerm", &PyMMFFMolProperties::setMMFFAngleTerm,
           "state"_a = true,
           "Sets the angle term to be included in the MMFF equation "
           "(defaults to True)")
      .def("GetMMFFAngleTerm", &PyMMFFMolProperties::getMMFFAngleTerm,
           "Returns whether the angle term is included in the MMFF equation.")
      .def("SetMMFFStretchBendTerm",
           &PyMMFFMolProperties::setMMFFStretchBendTerm, "state"_a = true,
           "Sets the stretch-bend term to be included in the MMFF equation "
           "(defaults to True)")
      .def("GetMMFFStretchBendTerm",
           &PyMMFFMolProperties::getMMFFStretchBendTerm,
           "Returns whether the stretch-bend term is included in the MMFF "
           "equation.")
      .def("SetMMFFOopTerm", &PyMMFFMolProperties::setMMFFOopTerm,
           "state"_a = true,
           "Sets the out-of-plane bend term to be included in the MMFF "
           "equation (defaults to True)")
      .def("GetMMFFOopTerm", &PyMMFFMolProperties::getMMFFOopTerm,
           "Returns whether the out-of-plane bend term is included in the "
           "MMFF equation.")
      .def("SetMMFFTorsionTerm", &PyMMFFMolProperties::setMMFFTorsionTerm,
           "state"_a = true,
           "Sets the torsional term to be included in the MMFF equation "
           "(defaults to True)")
      .def("GetMMFFTorsionTerm", &PyMMFFMolProperties::getMMFFTorsionTerm,
           "Returns whether the torsional term is included in the MMFF "
           "equation.")
      .def("SetMMFFVdWTerm", &PyMMFFMolProperties::setMMFFVdWTerm,
           "state"_a = true,
           "Sets the Van der Waals term to be included in the MMFF equation "
           "(defaults to True)")
      .def("GetMMFFVdWTerm", &PyMMFFMolProperties::getMMFFVdWTerm,
           "Returns whether the Van der Waals term is included in the MMFF "
           "equation.")
      .def("SetMMFFEleTerm", &PyMMFFMolProperties::setMMFFEleTerm,
           "state"_a = true,
           "Sets the electrostatic term to be included in the MMFF equation "
           "(defaults to True)")
      .def("GetMMFFEleTerm", &PyMMFFMolProperties::getMMFFEleTerm,
           "Returns whether the electrostatic term is included in the MMFF "
           "equation.")
      .def("SetMMFFVariant", &PyMMFFMolProperties::setMMFFVariant,
           "mmffVariant"_a = "MMFF94",
           "Sets the MMFF variant to be used (\"MMFF94\" or \"MMFF94s\"; "
           "defaults to \"MMFF94\")")
      .def("GetMMFFVariant", &PyMMFFMolProperties::getMMFFVariant,
           "Returns the currently configured MMFF variant "
           "(\"MMFF94\" or \"MMFF94s\").")
      .def("SetMMFFVerbosity", &PyMMFFMolProperties::setMMFFVerbosity,
           "verbosity"_a = 0,
           "Sets the MMFF verbosity (0: none; 1: low; 2: high; defaults to 0)");
}
