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
#include <sstream>

namespace python = boost::python;
using namespace RDKit;

namespace {

ROMol *normalizeHelper(MolStandardize::Normalizer &self, const ROMol &mol) {
  return self.normalize(mol);
}

void normalizeInPlaceHelper(MolStandardize::Normalizer &self, ROMol &mol) {
  self.normalizeInPlace(static_cast<RWMol &>(mol));
}

MolStandardize::Normalizer *normalizerFromDataAndParams(
    const std::string &data, const MolStandardize::CleanupParameters &params) {
  std::istringstream sstr(data);
  return new MolStandardize::Normalizer(sstr, params.maxRestarts);
}

}  // namespace

struct normalize_wrapper {
  static void wrap() {
    python::scope().attr("__doc__") =
        "Module containing tools for normalizing molecules defined by SMARTS "
        "patterns";

    std::string docString = "";

    python::class_<MolStandardize::Normalizer, boost::noncopyable>(
        "Normalizer", python::init<>(python::args("self")))
        .def(python::init<std::string, unsigned int>(
            python::args("self", "normalizeFilename", "maxRestarts")))
        .def("normalize", &normalizeHelper,
             (python::arg("self"), python::arg("mol")), "",
             python::return_value_policy<python::manage_new_object>())
        .def("normalizeInPlace", &normalizeInPlaceHelper,
             (python::arg("self"), python::arg("mol")),
             "modifies the input molecule");
    python::def(
        "NormalizerFromData", &normalizerFromDataAndParams,
        (python::arg("paramData"), python::arg("params")),
        "creates a Normalizer from a string containing normalization SMARTS",
        python::return_value_policy<python::manage_new_object>());
    python::def("NormalizerFromParams", &MolStandardize::normalizerFromParams,
                (python::arg("params")),
                "creates a Normalizer from CleanupParameters",
                python::return_value_policy<python::manage_new_object>());
  }
};

void wrap_normalize() { normalize_wrapper::wrap(); }
