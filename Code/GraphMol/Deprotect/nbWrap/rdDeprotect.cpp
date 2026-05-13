//
//  Copyright (C) 2020-2026 Brian P Kelley and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unique_ptr.h>
#include <nanobind/stl/vector.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Deprotect/Deprotect.h>
#include <RDBoost/Wrap_nb.h>

using namespace nb::literals;
using namespace RDKit;

NB_MODULE(rdDeprotect, m) {
  m.doc() = R"DOC(Module containing the Deprotect functionality for removing
protecting groups from molecules.)DOC";

  nb::class_<Deprotect::DeprotectData>(m, "DeprotectData",
                                       R"DOC(DeprotectData class, contains a single deprotection reaction and information

 deprotectdata.deprotection_class - functional group being protected
 deprotectdata.reaction_smarts - reaction smarts used for deprotection
 deprotectdata.abbreviation - common abbreviation for the protecting group
 deprotectdata.full_name - full name for the protecting group

)DOC")
      .def(nb::init<std::string, std::string, std::string, std::string>(),
           "deprotection_class"_a, "reaction_smarts"_a, "abbreviation"_a,
           "full_name"_a,
           R"DOC(Construct a new DeprotectData instance.
  >>> reaction_class = "amine"
  >>> reaction_smarts = "[C;R0][C;R0]([C;R0])([O;R0][C;R0](=[O;R0])[NX3;H0,H1:1])C>>[N:1]"
  >>> abbreviation = "Boc"
  >>> full_name = "tert-butyloxycarbonyl"
  >>> data = DeprotectData(reaction_class, reaction_smarts, abbreviation, full_name)
  >>> assert data.isValid()
)DOC")
      .def_ro("deprotection_class",
              &Deprotect::DeprotectData::deprotection_class)
      .def_ro("full_name", &Deprotect::DeprotectData::full_name)
      .def_ro("abbreviation", &Deprotect::DeprotectData::abbreviation)
      .def_ro("reaction_smarts", &Deprotect::DeprotectData::reaction_smarts)
      .def_ro("example", &Deprotect::DeprotectData::example)
      .def("isValid", &Deprotect::DeprotectData::isValid,
           "Returns True if the DeprotectData has a valid reaction");

  m.def(
      "GetDeprotections",
      []() { return Deprotect::getDeprotections(); },
      "Return the default list of deprotections");

  m.def(
      "Deprotect",
      [](const ROMol &mol, const nb::object &deprotections) {
        if (!deprotections.is_none()) {
          std::vector<Deprotect::DeprotectData> deprots;
          pythonObjectToVect<Deprotect::DeprotectData>(deprotections, deprots);
          return Deprotect::deprotect(mol, deprots);
        } else {
          return Deprotect::deprotect(mol, Deprotect::getDeprotections());
        }
      },
      "mol"_a, "deprotections"_a = nb::none(),
      "Return the deprotected version of the molecule.");

  m.def(
      "DeprotectInPlace",
      [](ROMol &mol, const nb::object &deprotections) {
        RWMol &rwmol = static_cast<RWMol &>(mol);
        if (!deprotections.is_none()) {
          std::vector<Deprotect::DeprotectData> deprots;
          pythonObjectToVect<Deprotect::DeprotectData>(deprotections, deprots);
          return Deprotect::deprotectInPlace(rwmol, deprots);
        } else {
          return Deprotect::deprotectInPlace(rwmol,
                                             Deprotect::getDeprotections());
        }
      },
      "mol"_a, "deprotections"_a = nb::none(),
      "Deprotects the molecule in place.");
}
