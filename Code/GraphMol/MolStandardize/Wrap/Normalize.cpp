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
#include <GraphMol/MolStandardize/Normalize.h>

namespace python = boost::python;
using namespace RDKit;

namespace {

ROMol *normalizeHelper(MolStandardize::Normalizer &self, const ROMol &mol) {
  return self.normalize(mol);
}

}  // namespace

struct normalize_wrapper {
  static void wrap() {
    python::scope().attr("__doc__") =
        "Module containing tools for normalizing molecules defined by SMARTS "
        "patterns";

    std::string docString = "";

    python::class_<MolStandardize::Normalizer, boost::noncopyable>(
        "Normalizer", python::init<>())
        .def(python::init<std::string, unsigned int>())
        .def("normalize", &normalizeHelper,
             (python::arg("self"), python::arg("mol")), "",
             python::return_value_policy<python::manage_new_object>());
  }
};

void wrap_normalize() { normalize_wrapper::wrap(); }
