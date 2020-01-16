//
//  Copyright (C) 2020 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/Wrap.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolStandardize/MolStandardize.h>
#include <GraphMol/MolStandardize/Tautomer.h>
#include <sstream>

namespace python = boost::python;
using namespace RDKit;

namespace {

// ROMol *normalizeHelper(MolStandardize::Normalizer &self, const ROMol &mol) {
//   return self.normalize(mol);
// }

MolStandardize::TautomerEnumerator *EnumeratorFromParams(
    const MolStandardize::CleanupParameters &params) {
  auto tautparams =
      new MolStandardize::TautomerCatalogParams(params.tautomerTransforms);
  auto cat = new MolStandardize::TautomerCatalog(tautparams);
  return new MolStandardize::TautomerEnumerator(cat);
}

MolStandardize::TautomerEnumerator *createDefaultEnumerator() {
  MolStandardize::CleanupParameters ps;
  return EnumeratorFromParams(ps);
}

ROMol *canonicalizeHelper(MolStandardize::TautomerEnumerator &self,
                          const ROMol &mol) {
  return self.canonicalize(mol);
}
}  // namespace

struct tautomer_wrapper {
  static void wrap() {
    std::string docString = "";

    python::class_<MolStandardize::TautomerEnumerator, boost::noncopyable>(
        "TautomerEnumerator", python::no_init)
        .def("__init__", python::make_constructor(createDefaultEnumerator))
        .def("__init__", python::make_constructor(EnumeratorFromParams))
        .def("Enumerate", &MolStandardize::TautomerEnumerator::enumerate,
             (python::arg("self"), python::arg("mol")), "")
        .def("Canonicalize", &canonicalizeHelper,
             (python::arg("self"), python::arg("mol")), "",
             python::return_value_policy<python::manage_new_object>());
    ;
  }
};

void wrap_tautomer() { tautomer_wrapper::wrap(); }
