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

class pyobjFunctor {
 public:
  pyobjFunctor(python::object obj) : dp_obj(std::move(obj)) {}
  ~pyobjFunctor() {}
  int operator()(const ROMol &m) {
    return python::extract<int>(dp_obj(boost::ref(m)));
  }

 private:
  python::object dp_obj;
};

ROMol *canonicalizeHelper2(MolStandardize::TautomerEnumerator &self,
                           const ROMol &mol, python::object scoreFunc) {
  pyobjFunctor ftor(scoreFunc);
  return self.canonicalize(mol, ftor);
}
}  // namespace

struct tautomer_wrapper {
  static void wrap() {
    python::class_<MolStandardize::TautomerEnumerator, boost::noncopyable>(
        "TautomerEnumerator", python::no_init)
        .def("__init__", python::make_constructor(createDefaultEnumerator))
        .def("__init__", python::make_constructor(EnumeratorFromParams))
        .def("Enumerate", &MolStandardize::TautomerEnumerator::enumerate,
             (python::arg("self"), python::arg("mol")),
             "generates the tautomers for a molecule")
        .def("Canonicalize", &canonicalizeHelper,
             (python::arg("self"), python::arg("mol")),
             "returns the canonical tautomer for a molecule",
             python::return_value_policy<python::manage_new_object>())
        .def(
            "Canonicalize", &canonicalizeHelper2,
            (python::arg("self"), python::arg("mol"), python::arg("scoreFunc")),
            "returns the canonical tautomer for a molecul using a custom "
            "scoring function",
            python::return_value_policy<python::manage_new_object>());
    ;
  }
};

void wrap_tautomer() { tautomer_wrapper::wrap(); }
