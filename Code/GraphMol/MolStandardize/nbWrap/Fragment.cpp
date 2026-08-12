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

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolStandardize/Fragment.h>

#include <sstream>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

namespace {

ROMol *removeHelper(MolStandardize::FragmentRemover &self, const ROMol &mol) {
  return self.remove(mol);
}

void removeInPlaceHelper(MolStandardize::FragmentRemover &self, ROMol &mol) {
  self.removeInPlace(static_cast<RWMol &>(mol));
}

ROMol *chooseHelper(MolStandardize::LargestFragmentChooser &self,
                    const ROMol &mol) {
  return self.choose(mol);
}

void chooseInPlaceHelper(MolStandardize::LargestFragmentChooser &self,
                         ROMol &mol) {
  self.chooseInPlace(static_cast<RWMol &>(mol));
}

MolStandardize::FragmentRemover *removerFromParams(const std::string &data,
                                                   bool leave_last,
                                                   bool skip_if_all_match) {
  std::istringstream sstr(data);
  return new MolStandardize::FragmentRemover(sstr, leave_last,
                                             skip_if_all_match);
}

}  // namespace

void wrap_fragment(nb::module_ &m) {
  nb::class_<MolStandardize::FragmentRemover>(m, "FragmentRemover")
      .def(nb::init<>())
      .def(nb::init<std::string, bool, bool>(), "fragmentFilename"_a = "",
           "leave_last"_a = true, "skip_if_all_match"_a = false)
      .def("remove", &removeHelper, "mol"_a, "",
           nb::rv_policy::take_ownership)
      .def("removeInPlace", &removeInPlaceHelper, "mol"_a,
           "modifies the molecule in place");

  m.def("FragmentRemoverFromData", &removerFromParams, "fragmentData"_a,
        "leave_last"_a = true, "skip_if_all_match"_a = false,
        "creates a FragmentRemover from a string containing parameter data",
        nb::rv_policy::take_ownership);

  nb::class_<MolStandardize::LargestFragmentChooser>(m,
                                                      "LargestFragmentChooser")
      .def(nb::init<bool>(), "preferOrganic"_a = false)
      .def(nb::init<const MolStandardize::CleanupParameters &>(), "params"_a)
      .def("choose", &chooseHelper, "mol"_a, "",
           nb::rv_policy::take_ownership)
      .def("chooseInPlace", &chooseInPlaceHelper, "mol"_a,
           "modifies the molecule in place");
}
