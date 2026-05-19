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
#include <GraphMol/MolStandardize/Normalize.h>
#include <sstream>

namespace nb = nanobind;
using namespace nb::literals;
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

void wrap_normalize(nb::module_ &m) {
  nb::class_<MolStandardize::Normalizer>(m, "Normalizer")
      .def(nb::init<>())
      .def(nb::init<std::string, unsigned int>(), "normalizeFilename"_a,
           "maxRestarts"_a)
      .def("normalize", &normalizeHelper, "mol"_a, "",
           nb::rv_policy::take_ownership)
      .def("normalizeInPlace", &normalizeInPlaceHelper, "mol"_a,
           "modifies the input molecule");

  m.def("NormalizerFromData", &normalizerFromDataAndParams, "paramData"_a,
        "params"_a,
        "creates a Normalizer from a string containing normalization SMARTS",
        nb::rv_policy::take_ownership);

  m.def("NormalizerFromParams", &MolStandardize::normalizerFromParams,
        "params"_a, "creates a Normalizer from CleanupParameters",
        nb::rv_policy::take_ownership);
}
