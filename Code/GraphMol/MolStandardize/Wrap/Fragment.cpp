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
#include <GraphMol/MolStandardize/Fragment.h>

namespace python = boost::python;
using namespace RDKit;

namespace {

ROMol *removeHelper(MolStandardize::FragmentRemover &self, const ROMol &mol) {
  return self.remove(mol);
}

ROMol *chooseHelper(MolStandardize::LargestFragmentChooser &self,
                    const ROMol &mol) {
  return self.choose(mol);
}
MolStandardize::FragmentRemover *removerFromParams(const std::string &data,
                                                   bool leave_last,
                                                   bool skip_if_all_match) {
  std::istringstream sstr(data);
  return new MolStandardize::FragmentRemover(sstr, leave_last,
                                             skip_if_all_match);
}

}  // namespace

struct fragment_wrapper {
  static void wrap() {
    python::scope().attr("__doc__") =
        "Module containing tools for dealing with molecules with more than \
					one covalently bonded unit";

    std::string docString = "";

    python::class_<MolStandardize::FragmentRemover, boost::noncopyable>(
        "FragmentRemover", python::init<>())
        .def(python::init<std::string, bool, bool>(
            (python::arg("fragmentFilename") = "",
             python::arg("leave_last") = true,
             python::arg("skip_if_all_match") = false)))
        .def("remove", &removeHelper, (python::arg("self"), python::arg("mol")),
             "", python::return_value_policy<python::manage_new_object>());

    python::def(
        "FragmentRemoverFromData", &removerFromParams,
        (python::arg("fragmentData"), python::arg("leave_last") = true,
         python::arg("skip_if_all_match") = false),
        "creates a FragmentRemover from a string containing parameter data",
        python::return_value_policy<python::manage_new_object>());

    python::class_<MolStandardize::LargestFragmentChooser, boost::noncopyable>(
        "LargestFragmentChooser",
        python::init<bool>(python::arg("preferOrganic") = false))
        .def(python::init<const MolStandardize::CleanupParameters &>(
            python::arg("params")))
        .def("choose", &chooseHelper, (python::arg("self"), python::arg("mol")),
             "", python::return_value_policy<python::manage_new_object>());
  }
};

void wrap_fragment() { fragment_wrapper::wrap(); }
