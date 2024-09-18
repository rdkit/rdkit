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
#include <GraphMol/MolStandardize/Charge.h>

namespace python = boost::python;
using namespace RDKit;

namespace {

std::vector<MolStandardize::ChargeCorrection> defaultChargeCorrections() {
  return MolStandardize::CHARGE_CORRECTIONS;
}

ROMol *reionizeHelper(MolStandardize::Reionizer &self, const ROMol &mol) {
  return self.reionize(mol);
}

void reionizeInPlaceHelper(MolStandardize::Reionizer &self, ROMol &mol) {
  self.reionizeInPlace(static_cast<RWMol &>(mol));
}
MolStandardize::Reionizer *reionizerFromData(const std::string &data,
                                             python::object chargeCorrections) {
  std::istringstream sstr(data);
  auto corrections =
      pythonObjectToVect<MolStandardize::ChargeCorrection>(chargeCorrections);
  MolStandardize::Reionizer *res;
  if (corrections) {
    res = new MolStandardize::Reionizer(sstr, *corrections);
  } else {
    res = new MolStandardize::Reionizer(
        sstr, std::vector<MolStandardize::ChargeCorrection>());
  }
  return res;
}

void unchargeInPlaceHelper(MolStandardize::Uncharger &self, ROMol &mol) {
  self.unchargeInPlace(static_cast<RWMol &>(mol));
}

}  // namespace

struct charge_wrapper {
  static void wrap() {
    python::scope().attr("__doc__") =
        "Module containing functions for charge corrections";

    std::string docString = "";

    python::class_<MolStandardize::ChargeCorrection, boost::noncopyable>(
        "ChargeCorrection",
        python::init<std::string, std::string, int>(
            python::args("self", "name", "smarts", "charge")))
        .def_readwrite("Name", &MolStandardize::ChargeCorrection::Name)
        .def_readwrite("Smarts", &MolStandardize::ChargeCorrection::Smarts)
        .def_readwrite("Charge", &MolStandardize::ChargeCorrection::Charge);

    python::def("CHARGE_CORRECTIONS", defaultChargeCorrections);

    python::class_<MolStandardize::Reionizer, boost::noncopyable>(
        "Reionizer", python::init<>(python::args("self")))
        .def(python::init<std::string>(python::args("self", "acidbaseFile")))
        .def(python::init<std::string,
                          std::vector<MolStandardize::ChargeCorrection>>(
            python::args("self", "acidbaseFile", "ccs")))
        .def("reionize", &reionizeHelper,
             (python::arg("self"), python::arg("mol")), "",
             python::return_value_policy<python::manage_new_object>())
        .def("reionizeInPlace", reionizeInPlaceHelper,
             (python::arg("self"), python::arg("mol")),
             "modifies the input molecule");

    python::def("ReionizerFromData", &reionizerFromData,
                (python::arg("paramData"),
                 python::arg("chargeCorrections") = python::list()),
                "creates a reionizer from a string containing parameter data "
                "and a list of charge corrections",
                python::return_value_policy<python::manage_new_object>());
    python::class_<MolStandardize::Uncharger, boost::noncopyable>(
        "Uncharger", python::init<bool, bool, bool>((python::arg("self"),
                                         python::arg("canonicalOrder") = true,
                                         python::arg("force") = false,
                                         python::arg("protonationOnly") = false)))
        .def("uncharge", &MolStandardize::Uncharger::uncharge,
             (python::arg("self"), python::arg("mol")), "",
             python::return_value_policy<python::manage_new_object>())
        .def("unchargeInPlace", unchargeInPlaceHelper,
             (python::arg("self"), python::arg("mol")),
             "modifies the input molecule");
  }
};

void wrap_charge() { charge_wrapper::wrap(); }
