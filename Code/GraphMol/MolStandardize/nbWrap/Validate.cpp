//
//  Copyright (C) 2018-2026 Susan H. Leung and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/trampoline.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolStandardize/Validate.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

namespace {

struct ValidationMethodTrampoline : MolStandardize::ValidationMethod {
  NB_TRAMPOLINE(MolStandardize::ValidationMethod, 2);

  std::vector<MolStandardize::ValidationErrorInfo> validate(
      const ROMol &mol, bool reportAllFailures) const override {
    NB_OVERRIDE_PURE(validate, mol, reportAllFailures);
  }

  std::shared_ptr<MolStandardize::ValidationMethod> copy() const override {
    NB_OVERRIDE_PURE(copy);
  }
};

// Wrap ValidationMethod::validate and convert the returned
// vector into a Python list of strings
nb::list pythonValidateMethod(const MolStandardize::ValidationMethod &self,
                              const ROMol &mol, bool reportAllFailures) {
  nb::list res;
  std::vector<MolStandardize::ValidationErrorInfo> errout =
      self.validate(mol, reportAllFailures);
  for (const auto &msg : errout) {
    res.append(msg);
  }
  return res;
}

MolStandardize::MolVSValidation *getMolVSValidation(nb::object validations) {
  std::vector<std::shared_ptr<MolStandardize::ValidationMethod>> vs;
  for (nb::handle h : validations) {
    auto v = nb::cast<std::shared_ptr<MolStandardize::ValidationMethod>>(h);
    vs.push_back(v->copy());
  }
  if (vs.empty()) {
    throw nb::value_error("validations argument must be non-empty");
  }
  return new MolStandardize::MolVSValidation(vs);
}

MolStandardize::AllowedAtomsValidation *getAllowedAtomsValidation(
    nb::object atoms) {
  std::vector<std::shared_ptr<Atom>> satoms;
  for (nb::handle h : atoms) {
    auto ap = nb::cast<Atom *>(h);
    satoms.push_back(std::shared_ptr<Atom>(ap->copy()));
  }
  if (satoms.empty()) {
    throw nb::value_error("allowedAtoms argument must be non-empty");
  }
  return new MolStandardize::AllowedAtomsValidation(satoms);
}

MolStandardize::DisallowedAtomsValidation *getDisallowedAtomsValidation(
    nb::object atoms) {
  std::vector<std::shared_ptr<Atom>> satoms;
  for (nb::handle h : atoms) {
    auto ap = nb::cast<Atom *>(h);
    satoms.push_back(std::shared_ptr<Atom>(ap->copy()));
  }
  if (satoms.empty()) {
    throw nb::value_error("disallowedAtoms must be non-empty");
  }
  return new MolStandardize::DisallowedAtomsValidation(satoms);
}

nb::list standardizeSmilesHelper(const std::string &smiles) {
  nb::list res;
  std::vector<MolStandardize::ValidationErrorInfo> errout =
      MolStandardize::validateSmiles(smiles);
  for (const auto &msg : errout) {
    res.append(msg);
  }
  return res;
}

}  // namespace

void wrap_validate(nb::module_ &m) {
  nb::class_<MolStandardize::ValidationMethod, ValidationMethodTrampoline>(
      m, "ValidationMethod")
      .def(nb::init<>())
      .def("validate", pythonValidateMethod, "mol"_a,
           "reportAllFailures"_a = false, "");

  nb::class_<MolStandardize::RDKitValidation,
             MolStandardize::ValidationMethod>(m, "RDKitValidation")
      .def(nb::init<bool>(), "allowEmptyMolecules"_a = false)
      .def_rw("allowEmptyMolecules",
               &MolStandardize::RDKitValidation::allowEmptyMolecules);

  nb::class_<MolStandardize::NoAtomValidation,
             MolStandardize::ValidationMethod>(m, "NoAtomValidation")
      .def(nb::init<>());

  nb::class_<MolStandardize::FragmentValidation,
             MolStandardize::ValidationMethod>(m, "FragmentValidation")
      .def(nb::init<>());

  nb::class_<MolStandardize::NeutralValidation,
             MolStandardize::ValidationMethod>(m, "NeutralValidation")
      .def(nb::init<>());

  nb::class_<MolStandardize::IsotopeValidation,
             MolStandardize::ValidationMethod>(m, "IsotopeValidation")
      .def(nb::init<bool>(), "strict"_a = false)
      .def_rw("strict", &MolStandardize::IsotopeValidation::strict);

  nb::class_<MolStandardize::MolVSValidation,
             MolStandardize::ValidationMethod>(m, "MolVSValidation")
      .def(nb::init<>())
      .def("__init__",
           [](MolStandardize::MolVSValidation *self, nb::object validations) {
             std::unique_ptr<MolStandardize::MolVSValidation> v(
                 getMolVSValidation(validations));
             new (self) MolStandardize::MolVSValidation(*v);
           },
           "validations"_a);

  nb::class_<MolStandardize::AllowedAtomsValidation,
             MolStandardize::ValidationMethod>(m, "AllowedAtomsValidation")
      .def("__init__",
           [](MolStandardize::AllowedAtomsValidation *self,
              nb::object atoms) {
             std::unique_ptr<MolStandardize::AllowedAtomsValidation> v(
                 getAllowedAtomsValidation(atoms));
             new (self) MolStandardize::AllowedAtomsValidation(*v);
           },
           "atoms"_a);

  nb::class_<MolStandardize::DisallowedAtomsValidation,
             MolStandardize::ValidationMethod>(m, "DisallowedAtomsValidation")
      .def("__init__",
           [](MolStandardize::DisallowedAtomsValidation *self,
              nb::object atoms) {
             std::unique_ptr<MolStandardize::DisallowedAtomsValidation> v(
                 getDisallowedAtomsValidation(atoms));
             new (self) MolStandardize::DisallowedAtomsValidation(*v);
           },
           "atoms"_a);

  nb::class_<MolStandardize::FeaturesValidation,
             MolStandardize::ValidationMethod>(m, "FeaturesValidation")
      .def(nb::init<bool, bool, bool, bool, bool, bool>(),
           "allowEnhancedStereo"_a = false, "allowAromaticBondType"_a = false,
           "allowDativeBondType"_a = false, "allowQueries"_a = false,
           "allowDummies"_a = false, "allowAtomAliases"_a = false)
      .def_rw("allowEnhancedStereo",
               &MolStandardize::FeaturesValidation::allowEnhancedStereo)
      .def_rw("allowAromaticBondType",
               &MolStandardize::FeaturesValidation::allowAromaticBondType)
      .def_rw("allowDativeBondType",
               &MolStandardize::FeaturesValidation::allowDativeBondType)
      .def_rw("allowQueries",
               &MolStandardize::FeaturesValidation::allowQueries)
      .def_rw("allowDummies",
               &MolStandardize::FeaturesValidation::allowDummies)
      .def_rw("allowAtomAliases",
               &MolStandardize::FeaturesValidation::allowAtomAliases);

  nb::class_<MolStandardize::DisallowedRadicalValidation,
             MolStandardize::ValidationMethod>(m,
                                               "DisallowedRadicalValidation")
      .def(nb::init<>());

  nb::class_<MolStandardize::Is2DValidation,
             MolStandardize::ValidationMethod>(m, "Is2DValidation")
      .def(nb::init<double>(), "threshold"_a = 1e-3)
      .def_rw("threshold", &MolStandardize::Is2DValidation::threshold);

  nb::class_<MolStandardize::Layout2DValidation,
             MolStandardize::ValidationMethod>(m, "Layout2DValidation")
      .def(nb::init<double, double, bool, bool, double>(),
           "clashLimit"_a = 0.15, "bondLengthLimit"_a = 25.,
           "allowLongBondsInRings"_a = true,
           "allowAtomBondClashExemption"_a = true,
           "minMedianBondLength"_a = false)
      .def_rw("clashLimit", &MolStandardize::Layout2DValidation::clashLimit)
      .def_rw("bondLengthLimit",
               &MolStandardize::Layout2DValidation::bondLengthLimit)
      .def_rw("allowLongBondsInRings",
               &MolStandardize::Layout2DValidation::allowLongBondsInRings)
      .def_rw(
          "allowAtomBondClashExemption",
          &MolStandardize::Layout2DValidation::allowAtomBondClashExemption)
      .def_rw("minMedianBondLength",
               &MolStandardize::Layout2DValidation::minMedianBondLength);

  nb::class_<MolStandardize::StereoValidation,
             MolStandardize::ValidationMethod>(m, "StereoValidation")
      .def(nb::init<>());

  m.def("ValidateSmiles", standardizeSmilesHelper, "mol"_a, "");
}
