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

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolStandardize/Charge.h>

#include <sstream>

namespace nb = nanobind;
using namespace nb::literals;
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
                                             nb::object chargeCorrections) {
  std::istringstream sstr(data);
  std::vector<MolStandardize::ChargeCorrection> corrections;
  if (!chargeCorrections.is_none()) {
    for (nb::handle h : chargeCorrections) {
      corrections.push_back(nb::cast<MolStandardize::ChargeCorrection>(h));
    }
  }
  return new MolStandardize::Reionizer(sstr, corrections);
}

void unchargeInPlaceHelper(MolStandardize::Uncharger &self, ROMol &mol) {
  self.unchargeInPlace(static_cast<RWMol &>(mol));
}

}  // namespace

void wrap_charge(nb::module_ &m) {
  nb::class_<MolStandardize::ChargeCorrection>(m, "ChargeCorrection")
      .def(nb::init<std::string, std::string, int>(), "name"_a, "smarts"_a,
           "charge"_a)
      .def_rw("Name", &MolStandardize::ChargeCorrection::Name)
      .def_rw("Smarts", &MolStandardize::ChargeCorrection::Smarts)
      .def_rw("Charge", &MolStandardize::ChargeCorrection::Charge);

  m.def("CHARGE_CORRECTIONS", defaultChargeCorrections);

  nb::class_<MolStandardize::Reionizer>(m, "Reionizer")
      .def(nb::init<>())
      .def(nb::init<std::string>(), "acidbaseFile"_a)
      .def(nb::init<std::string,
                    std::vector<MolStandardize::ChargeCorrection>>(),
           "acidbaseFile"_a, "ccs"_a)
      .def("reionize", &reionizeHelper, "mol"_a, "",
           nb::rv_policy::take_ownership)
      .def("reionizeInPlace", reionizeInPlaceHelper, "mol"_a,
           "modifies the input molecule");

  m.def("ReionizerFromData", &reionizerFromData, "paramData"_a,
        "chargeCorrections"_a = nb::none(),
        "creates a reionizer from a string containing parameter data "
        "and a list of charge corrections",
        nb::rv_policy::take_ownership);

  nb::class_<MolStandardize::Uncharger>(m, "Uncharger")
      .def(nb::init<bool, bool, bool>(), "canonicalOrder"_a = true,
           "force"_a = false, "protonationOnly"_a = false)
      .def("uncharge", &MolStandardize::Uncharger::uncharge, "mol"_a, "",
           nb::rv_policy::take_ownership)
      .def("unchargeInPlace", unchargeInPlaceHelper, "mol"_a,
           "modifies the input molecule");
}
