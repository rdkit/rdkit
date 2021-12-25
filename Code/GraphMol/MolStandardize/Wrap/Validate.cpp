//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/Wrap.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolStandardize/Validate.h>

namespace python = boost::python;
using namespace RDKit;

namespace {

python::list rdkitValidate(MolStandardize::RDKitValidation &self,
                           const ROMol &mol, const bool reportAllFailures) {
  python::list res;
  std::vector<MolStandardize::ValidationErrorInfo> errout =
      self.validate(mol, reportAllFailures);
  for (auto &query : errout) {
    std::string msg = query.what();
    res.append(msg);
  }
  return res;
}

MolStandardize::MolVSValidation *getMolVSValidation(
    python::object validations) {
  std::vector<boost::shared_ptr<MolStandardize::MolVSValidations>> vs;

  auto pvect =
      pythonObjectToVect<boost::shared_ptr<MolStandardize::MolVSValidations>>(
          validations);
  if (!pvect) {
    throw_value_error("validations argument must be non-empty");
  }
  for (auto v : *pvect) {
    vs.push_back(v->copy());
  }
  return new MolStandardize::MolVSValidation(vs);
}

python::list molVSvalidateHelper(MolStandardize::MolVSValidation &self,
                                 const ROMol &mol, bool reportAllFailures) {
  python::list s;
  std::vector<MolStandardize::ValidationErrorInfo> errout =
      self.validate(mol, reportAllFailures);
  for (auto &query : errout) {
    s.append(query.what());
  }
  return s;
}

MolStandardize::AllowedAtomsValidation *getAllowedAtomsValidation(
    python::object atoms) {
  auto p_atomList = pythonObjectToVect<Atom *>(atoms);
  if (!p_atomList) {
    throw_value_error("allowedAtoms argument must be non-empty");
  }
  std::vector<std::shared_ptr<Atom>> satoms;
  for (auto ap : *p_atomList) {
    satoms.push_back(std::shared_ptr<Atom>(ap->copy()));
  }
  return new MolStandardize::AllowedAtomsValidation(satoms);
}

python::list allowedAtomsValidate(MolStandardize::AllowedAtomsValidation &self,
                                  const ROMol &mol,
                                  const bool reportAllFailures) {
  python::list res;
  std::vector<MolStandardize::ValidationErrorInfo> errout =
      self.validate(mol, reportAllFailures);
  for (auto &query : errout) {
    std::string msg = query.what();
    res.append(msg);
  }
  return res;
}

MolStandardize::DisallowedAtomsValidation *getDisallowedAtomsValidation(
    python::object atoms) {
  auto p_atomList = pythonObjectToVect<Atom *>(atoms);
  if (!p_atomList) {
    throw_value_error("disallowedAtoms must be non-empty");
  }
  std::vector<std::shared_ptr<Atom>> satoms;
  for (auto ap : *p_atomList) {
    satoms.push_back(std::shared_ptr<Atom>(ap->copy()));
  }
  return new MolStandardize::DisallowedAtomsValidation(satoms);
}

python::list disallowedAtomsValidate(
    MolStandardize::DisallowedAtomsValidation &self, const ROMol &mol,
    const bool reportAllFailures) {
  python::list res;
  std::vector<MolStandardize::ValidationErrorInfo> errout =
      self.validate(mol, reportAllFailures);
  for (auto &query : errout) {
    std::string msg = query.what();
    res.append(msg);
  }
  return res;
}

python::list standardizeSmilesHelper(const std::string &smiles) {
  python::list res;
  std::vector<MolStandardize::ValidationErrorInfo> errout =
      MolStandardize::validateSmiles(smiles);
  for (auto &query : errout) {
    std::string msg = query.what();
    res.append(msg);
  }
  return res;
}

}  // namespace

struct validate_wrapper {
  static void wrap() {
    std::string docString = "";

    python::class_<MolStandardize::RDKitValidation, boost::noncopyable>(
        "RDKitValidation", python::init<>())
        .def("validate", rdkitValidate,
             (python::arg("self"), python::arg("mol"),
              python::arg("reportAllFailures") = false),
             "");

    python::class_<MolStandardize::MolVSValidations, boost::noncopyable>(
        "MolVSValidations", python::no_init)
        .def("run", &MolStandardize::MolVSValidations::run,
             (python::arg("self"), python::arg("mol"),
              python::arg("reportAllFailures"), python::arg("errors")),
             "");

    python::class_<MolStandardize::NoAtomValidation,
                   python::bases<MolStandardize::MolVSValidations>>(
        "NoAtomValidation", python::init<>())
        .def("run", &MolStandardize::NoAtomValidation::run,
             (python::arg("self"), python::arg("mol"),
              python::arg("reportAllFailures"), python::arg("errors")),
             "");

    python::class_<MolStandardize::FragmentValidation,
                   python::bases<MolStandardize::MolVSValidations>>(
        "FragmentValidation", python::init<>())
        .def("run", &MolStandardize::FragmentValidation::run,
             (python::arg("self"), python::arg("mol"),
              python::arg("reportAllFailures"), python::arg("errors")),
             "");
    python::class_<MolStandardize::NeutralValidation,
                   python::bases<MolStandardize::MolVSValidations>>(
        "NeutralValidation", python::init<>())
        .def("run", &MolStandardize::NeutralValidation::run,
             (python::arg("self"), python::arg("mol"),
              python::arg("reportAllFailures"), python::arg("errors")),
             "");
    python::class_<MolStandardize::IsotopeValidation,
                   python::bases<MolStandardize::MolVSValidations>>(
        "IsotopeValidation", python::init<>())
        .def("run", &MolStandardize::IsotopeValidation::run,
             (python::arg("self"), python::arg("mol"),
              python::arg("reportAllFailures"), python::arg("errors")),
             "");

    python::class_<MolStandardize::MolVSValidation, boost::noncopyable>(
        "MolVSValidation")
        .def("__init__", python::make_constructor(&getMolVSValidation))
        .def("validate", molVSvalidateHelper,
             (python::arg("self"), python::arg("mol"),
              python::arg("reportAllFailures") = false),
             "");

    python::class_<MolStandardize::AllowedAtomsValidation, boost::noncopyable>(
        "AllowedAtomsValidation", python::no_init)
        .def("__init__", python::make_constructor(&getAllowedAtomsValidation))
        .def("validate", allowedAtomsValidate,
             (python::arg("self"), python::arg("mol"),
              python::arg("reportAllFailures") = false),
             "");
    ;

    python::class_<MolStandardize::DisallowedAtomsValidation,
                   boost::noncopyable>("DisallowedAtomsValidation",
                                       python::no_init)
        .def("__init__",
             python::make_constructor(&getDisallowedAtomsValidation))
        .def("validate", disallowedAtomsValidate,
             (python::arg("self"), python::arg("mol"),
              python::arg("reportAllFailures") = false),
             "");

    python::def("ValidateSmiles", standardizeSmilesHelper, (python::arg("mol")),
                docString.c_str());
  }
};

void wrap_validate() { validate_wrapper::wrap(); }
