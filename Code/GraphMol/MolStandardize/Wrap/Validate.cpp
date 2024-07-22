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

struct ValidationMethodWrap
    : MolStandardize::ValidationMethod,
      python::wrapper<MolStandardize::ValidationMethod> {
  std::vector<MolStandardize::ValidationErrorInfo> validate(
      const ROMol &mol, bool reportAllFailures) const override {
    return this->get_override("validate")(mol, reportAllFailures);
  }

  std::shared_ptr<MolStandardize::ValidationMethod> copy() const override {
    return this->get_override("copy")();
  }
};

// Wrap ValidationMethod::validate and convert the returned
// vector into a python list of strings
python::list pythonValidateMethod(const MolStandardize::ValidationMethod &self,
                                  const ROMol &mol, bool reportAllFailures) {
  python::list res;
  std::vector<MolStandardize::ValidationErrorInfo> errout =
      self.validate(mol, reportAllFailures);
  for (const auto &msg : errout) {
    res.append(msg);
  }
  return res;
}

MolStandardize::MolVSValidation *getMolVSValidation(
    python::object validations) {
  std::vector<std::shared_ptr<MolStandardize::ValidationMethod>> vs;

  auto pvect =
      pythonObjectToVect<std::shared_ptr<MolStandardize::ValidationMethod>>(
          validations);
  if (!pvect) {
    throw_value_error("validations argument must be non-empty");
  }
  for (auto v : *pvect) {
    vs.push_back(v->copy());
  }
  return new MolStandardize::MolVSValidation(vs);
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

python::list standardizeSmilesHelper(const std::string &smiles) {
  python::list res;
  std::vector<MolStandardize::ValidationErrorInfo> errout =
      MolStandardize::validateSmiles(smiles);
  for (const auto &msg : errout) {
    res.append(msg);
  }
  return res;
}

}  // namespace

struct validate_wrapper {
  static void wrap() {
    std::string docString = "";

    python::class_<ValidationMethodWrap, boost::noncopyable>("ValidationMethod")
        .def("validate", pythonValidateMethod,
             (python::arg("self"), python::arg("mol"),
              python::arg("reportAllFailures") = false),
             "");

    python::class_<MolStandardize::RDKitValidation,
                   python::bases<MolStandardize::ValidationMethod>,
                   boost::noncopyable>("RDKitValidation");

    python::class_<MolStandardize::NoAtomValidation,
                   python::bases<MolStandardize::ValidationMethod>,
                   boost::noncopyable>("NoAtomValidation");

    python::class_<MolStandardize::FragmentValidation,
                   python::bases<MolStandardize::ValidationMethod>,
                   boost::noncopyable>("FragmentValidation");

    python::class_<MolStandardize::NeutralValidation,
                   python::bases<MolStandardize::ValidationMethod>,
                   boost::noncopyable>("NeutralValidation");

    python::class_<MolStandardize::IsotopeValidation,
                   python::bases<MolStandardize::ValidationMethod>,
                   boost::noncopyable>("IsotopeValidation")
        .def(python::init<bool>(python::arg("strict") = false))
        .def_readwrite("strict", &MolStandardize::IsotopeValidation::strict);

    python::class_<MolStandardize::MolVSValidation,
                   python::bases<MolStandardize::ValidationMethod>,
                   boost::noncopyable>("MolVSValidation")
        .def("__init__", python::make_constructor(&getMolVSValidation));

    python::class_<MolStandardize::AllowedAtomsValidation,
                   python::bases<MolStandardize::ValidationMethod>,
                   boost::noncopyable>("AllowedAtomsValidation",
                                       python::no_init)
        .def("__init__", python::make_constructor(&getAllowedAtomsValidation));

    python::class_<MolStandardize::DisallowedAtomsValidation,
                   python::bases<MolStandardize::ValidationMethod>,
                   boost::noncopyable>("DisallowedAtomsValidation",
                                       python::no_init)
        .def("__init__",
             python::make_constructor(&getDisallowedAtomsValidation));

    python::class_<MolStandardize::FeaturesValidation,
                   python::bases<MolStandardize::ValidationMethod>>(
        "FeaturesValidation")
        .def(python::init<bool, bool, bool, bool, bool, bool>(
            (python::arg("allowEnhancedStereo") = false,
             python::arg("allowAromaticBondType") = false,
             python::arg("allowDativeBondType") = false,
             python::arg("allowQueries") = false,
             python::arg("allowDummmies") = false,
             python::arg("allowAtomAliases") = false)))
        .def_readwrite("allowEnhancedStereo",
                       &MolStandardize::FeaturesValidation::allowEnhancedStereo)
        .def_readwrite(
            "allowAromaticBondType",
            &MolStandardize::FeaturesValidation::allowAromaticBondType)
        .def_readwrite("allowDativeBondType",
                       &MolStandardize::FeaturesValidation::allowDativeBondType)
        .def_readwrite("allowQueries",
                       &MolStandardize::FeaturesValidation::allowQueries)
        .def_readwrite("allowDummies",
                       &MolStandardize::FeaturesValidation::allowDummies)
        .def_readwrite("allowAtomAliases",
                       &MolStandardize::FeaturesValidation::allowAtomAliases);

    python::class_<MolStandardize::DisallowedRadicalValidation,
                   python::bases<MolStandardize::ValidationMethod>,
                   boost::noncopyable>("DisallowedRadicalValidation");

    python::class_<MolStandardize::Is2DValidation,
                   python::bases<MolStandardize::ValidationMethod>,
                   boost::noncopyable>("Is2DValidation")
        .def(python::init<double>(python::arg("threshold") = 1e-3))
        .def_readwrite("threshold", &MolStandardize::Is2DValidation::threshold);

    python::class_<MolStandardize::Layout2DValidation,
                   python::bases<MolStandardize::ValidationMethod>,
                   boost::noncopyable>("Layout2DValidation")
        .def(python::init<double, double, bool, bool, double>(
            (python::arg("clashLimit") = 0.15,
             python::arg("bondLengthLimit") = 25.,
             python::arg("allowLongBondsInRings") = true,
             python::arg("allowAtomBondClashExemption") = true,
             python::arg("minMedianBondLength") = false)))
        .def_readwrite("clashLimit",
                       &MolStandardize::Layout2DValidation::clashLimit)
        .def_readwrite("bondLengthLimit",
                       &MolStandardize::Layout2DValidation::bondLengthLimit)
        .def_readwrite(
            "allowLongBondsInRings",
            &MolStandardize::Layout2DValidation::allowLongBondsInRings)
        .def_readwrite(
            "allowAtomBondClashExemption",
            &MolStandardize::Layout2DValidation::allowAtomBondClashExemption)
        .def_readwrite(
            "minMedianBondLength",
            &MolStandardize::Layout2DValidation::minMedianBondLength);

    python::class_<MolStandardize::StereoValidation,
                   python::bases<MolStandardize::ValidationMethod>,
                   boost::noncopyable>("StereoValidation");

    python::def("ValidateSmiles", standardizeSmilesHelper, (python::arg("mol")),
                docString.c_str());
  }
};

void wrap_validate() { validate_wrapper::wrap(); }
