// $Id$
//
//  Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#include <ForceField/ForceField.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <ForceField/MMFF/Params.h>
#include <GraphMol/Trajectory/Snapshot.h>
#include <boost/python/tuple.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>
#include <algorithm>
#include <Geometry/point.h>

namespace python = boost::python;
namespace ForceFields {
class PyForceField {
 public:
  PyForceField(ForceField *f) : field(f) {}

  ~PyForceField() {
    // std::cerr << " *** destroy PyForce field " << std::endl;
    field.reset();
    // std::cerr << " ***       reset DONE" << std::endl;
    extraPoints.clear();
    // std::cerr << " *** destroy PyForce field DONE" << std::endl;
  }

  int addExtraPoint(double x, double y, double z, bool fixed = true) {
    PRECONDITION(this->field, "no force field");
    RDGeom::Point3D *pt = new RDGeom::Point3D(x, y, z);
    this->extraPoints.push_back(boost::shared_ptr<RDGeom::Point3D>(pt));
    unsigned int ptIdx = this->extraPoints.size() - 1;
    RDGeom::Point3D *ptr = this->extraPoints[ptIdx].get();
    this->field->positions().push_back(ptr);
    int idx = this->field->positions().size();
    if (fixed) {
      this->field->fixedPoints().push_back(idx - 1);
    }
    return idx;
  }

  double calcEnergyWithPos(const python::object &pos = python::object());

  double calcEnergy() { return calcEnergyWithPos(); }

  PyObject *calcGradWithPos(const python::object &pos = python::object());

  PyObject *positions();

  int minimize(int maxIts, double forceTol, double energyTol) {
    PRECONDITION(this->field, "no force field");
    return this->field->minimize(maxIts, forceTol, energyTol);
  }

  boost::python::tuple minimizeTrajectory(unsigned int snapshotFreq, int maxIts,
                                          double forceTol, double energyTol);

  void initialize() {
    PRECONDITION(this->field, "no force field");
    this->field->initialize();
  }

  unsigned int dimension() {
    PRECONDITION(this->field, "no force field");
    return this->field->dimension();
  }

  unsigned int numPoints() {
    PRECONDITION(this->field, "no force field");
    return this->field->numPoints();
  }

  // private:
  std::vector<boost::shared_ptr<RDGeom::Point3D>> extraPoints;
  boost::shared_ptr<ForceField> field;
};

class PyMMFFMolProperties {
 public:
  PyMMFFMolProperties(RDKit::MMFF::MMFFMolProperties *mp)
      : mmffMolProperties(mp) {}
  ~PyMMFFMolProperties() = default;

  unsigned int getMMFFAtomType(unsigned int idx) const {
    return (unsigned int)(mmffMolProperties->getMMFFAtomType(idx));
  }
  double getMMFFFormalCharge(unsigned int idx) const {
    return mmffMolProperties->getMMFFFormalCharge(idx);
  }
  double getMMFFPartialCharge(unsigned int idx) const {
    return mmffMolProperties->getMMFFPartialCharge(idx);
  }
  PyObject *getMMFFBondStretchParams(const RDKit::ROMol &mol,
                                     const unsigned int idx1,
                                     const unsigned int idx2) const;
  PyObject *getMMFFAngleBendParams(const RDKit::ROMol &mol,
                                   const unsigned int idx1,
                                   const unsigned int idx2,
                                   const unsigned int idx3) const;
  PyObject *getMMFFStretchBendParams(const RDKit::ROMol &mol,
                                     const unsigned int idx1,
                                     const unsigned int idx2,
                                     const unsigned int idx3) const;
  PyObject *getMMFFTorsionParams(const RDKit::ROMol &mol,
                                 const unsigned int idx1,
                                 const unsigned int idx2,
                                 const unsigned int idx3,
                                 const unsigned int idx4) const;
  PyObject *getMMFFOopBendParams(const RDKit::ROMol &mol,
                                 const unsigned int idx1,
                                 const unsigned int idx2,
                                 const unsigned int idx3,
                                 const unsigned int idx4) const;
  PyObject *getMMFFVdWParams(const unsigned int idx1,
                             const unsigned int idx2) const;
  void setMMFFDielectricModel(std::uint8_t dielModel) {
    mmffMolProperties->setMMFFDielectricModel(dielModel);
  }
  std::uint8_t getMMFFDielectricModel() const {
    return mmffMolProperties->getMMFFDielectricModel();
  }
  void setMMFFDielectricConstant(double dielConst) {
    mmffMolProperties->setMMFFDielectricConstant(dielConst);
  }
  double getMMFFDielectricConstant() const {
    return mmffMolProperties->getMMFFDielectricConstant();
  }
  void setMMFFBondTerm(bool state) {
    mmffMolProperties->setMMFFBondTerm(state);
  }
  bool getMMFFBondTerm() const { return mmffMolProperties->getMMFFBondTerm(); }
  void setMMFFAngleTerm(const bool state) {
    mmffMolProperties->setMMFFAngleTerm(state);
  }
  bool getMMFFAngleTerm() const {
    return mmffMolProperties->getMMFFAngleTerm();
  }
  void setMMFFStretchBendTerm(const bool state) {
    mmffMolProperties->setMMFFStretchBendTerm(state);
  }
  bool getMMFFStretchBendTerm() const {
    return mmffMolProperties->getMMFFStretchBendTerm();
  }
  void setMMFFOopTerm(const bool state) {
    mmffMolProperties->setMMFFOopTerm(state);
  }
  bool getMMFFOopTerm() const { return mmffMolProperties->getMMFFOopTerm(); }
  void setMMFFTorsionTerm(const bool state) {
    mmffMolProperties->setMMFFTorsionTerm(state);
  }
  bool getMMFFTorsionTerm() const {
    return mmffMolProperties->getMMFFTorsionTerm();
  }
  void setMMFFVdWTerm(const bool state) {
    mmffMolProperties->setMMFFVdWTerm(state);
  }
  bool getMMFFVdWTerm() const { return mmffMolProperties->getMMFFVdWTerm(); }
  void setMMFFEleTerm(const bool state) {
    mmffMolProperties->setMMFFEleTerm(state);
  }
  bool getMMFFEleTerm() const { return mmffMolProperties->getMMFFEleTerm(); }
  void setMMFFVariant(const std::string &mmffVariant) {
    mmffMolProperties->setMMFFVariant(mmffVariant);
  }
  std::string getMMFFVariant() const {
    return mmffMolProperties->getMMFFVariant();
  }
  void setMMFFVerbosity(unsigned int verbosity) {
    mmffMolProperties->setMMFFVerbosity(verbosity);
  }
  boost::shared_ptr<RDKit::MMFF::MMFFMolProperties> mmffMolProperties;
};
PyObject *getUFFBondStretchParams(const RDKit::ROMol &mol,
                                  const unsigned int idx1,
                                  const unsigned int idx2);
PyObject *getUFFAngleBendParams(const RDKit::ROMol &mol,
                                const unsigned int idx1,
                                const unsigned int idx2,
                                const unsigned int idx3);
PyObject *getUFFTorsionParams(const RDKit::ROMol &mol, const unsigned int idx1,
                              const unsigned int idx2, const unsigned int idx3,
                              const unsigned int idx4);
PyObject *getUFFInversionParams(const RDKit::ROMol &mol,
                                const unsigned int idx1,
                                const unsigned int idx2,
                                const unsigned int idx3,
                                const unsigned int idx4);
PyObject *getUFFVdWParams(const RDKit::ROMol &mol, const unsigned int idx1,
                          const unsigned int idx2);
}  // namespace ForceFields
