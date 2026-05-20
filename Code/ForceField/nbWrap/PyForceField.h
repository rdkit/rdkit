//
//  Copyright (C) 2005-2026 Rational Discovery LLC and other RDKit contributors
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
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <memory>
#include <vector>
#include <algorithm>
#include <Geometry/point.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace ForceFields {
class PyForceField {
 public:
  PyForceField(std::unique_ptr<ForceField> f) : field(std::move(f)) {}

  ~PyForceField() = default;

  int addExtraPoint(double x, double y, double z, bool fixed = true) {
    PRECONDITION(this->field, "no force field");
    this->extraPoints.push_back(std::make_unique<RDGeom::Point3D>(x, y, z));
    unsigned int ptIdx = this->extraPoints.size() - 1;
    RDGeom::Point3D *ptr = this->extraPoints[ptIdx].get();
    this->field->positions().push_back(ptr);
    int idx = this->field->positions().size();
    if (fixed) {
      this->field->fixedPoints().push_back(idx - 1);
    }
    return idx;
  }

  double calcEnergyWithPos(nb::object pos = nb::none());

  double calcEnergy() { return calcEnergyWithPos(); }

  nb::tuple calcGradWithPos(nb::object pos = nb::none());

  nb::tuple positions();

  int minimize(int maxIts, double forceTol, double energyTol) {
    PRECONDITION(this->field, "no force field");
    return this->field->minimize(maxIts, forceTol, energyTol);
  }

  nb::tuple minimizeTrajectory(unsigned int snapshotFreq, int maxIts,
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
  std::vector<std::unique_ptr<RDGeom::Point3D>> extraPoints;
  std::unique_ptr<ForceField> field;
};

class PyMMFFMolProperties {
 public:
  PyMMFFMolProperties(std::unique_ptr<RDKit::MMFF::MMFFMolProperties> mp)
      : mmffMolProperties(std::move(mp)) {}
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
  nb::object getMMFFBondStretchParams(const RDKit::ROMol &mol,
                                      const unsigned int idx1,
                                      const unsigned int idx2) const;
  nb::object getMMFFAngleBendParams(const RDKit::ROMol &mol,
                                    const unsigned int idx1,
                                    const unsigned int idx2,
                                    const unsigned int idx3) const;
  nb::object getMMFFStretchBendParams(const RDKit::ROMol &mol,
                                      const unsigned int idx1,
                                      const unsigned int idx2,
                                      const unsigned int idx3) const;
  nb::object getMMFFTorsionParams(const RDKit::ROMol &mol,
                                  const unsigned int idx1,
                                  const unsigned int idx2,
                                  const unsigned int idx3,
                                  const unsigned int idx4) const;
  nb::object getMMFFOopBendParams(const RDKit::ROMol &mol,
                                  const unsigned int idx1,
                                  const unsigned int idx2,
                                  const unsigned int idx3,
                                  const unsigned int idx4) const;
  nb::object getMMFFVdWParams(const unsigned int idx1,
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
  std::unique_ptr<RDKit::MMFF::MMFFMolProperties> mmffMolProperties;
};
nb::object getUFFBondStretchParams(const RDKit::ROMol &mol,
                                   const unsigned int idx1,
                                   const unsigned int idx2);
nb::object getUFFAngleBendParams(const RDKit::ROMol &mol,
                                 const unsigned int idx1,
                                 const unsigned int idx2,
                                 const unsigned int idx3);
nb::object getUFFTorsionParams(const RDKit::ROMol &mol,
                               const unsigned int idx1,
                               const unsigned int idx2,
                               const unsigned int idx3,
                               const unsigned int idx4);
nb::object getUFFInversionParams(const RDKit::ROMol &mol,
                                 const unsigned int idx1,
                                 const unsigned int idx2,
                                 const unsigned int idx3,
                                 const unsigned int idx4);
nb::object getUFFVdWParams(const RDKit::ROMol &mol, const unsigned int idx1,
                           const unsigned int idx2);
}  // namespace ForceFields
