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

ROMol *canonicalizeHelper(const MolStandardize::TautomerEnumerator &self,
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

ROMol *canonicalizeHelper2(const MolStandardize::TautomerEnumerator &self,
                           const ROMol &mol, python::object scoreFunc) {
  pyobjFunctor ftor(scoreFunc);
  return self.canonicalize(mol, ftor);
}

double scoreHelper(const MolStandardize::TautomerEnumerator &self,
                   const ROMol &mol) {
  RDUNUSED_PARAM(self);
  return MolStandardize::TautomerScoringFunctions::scoreTautomer(mol);
}

std::vector<ROMOL_SPTR> enumerateHelper(
    const MolStandardize::TautomerEnumerator &self, const ROMol &mol,
    python::object pyModAtoms, python::object pyModBonds) {
  if (pyModAtoms != python::object() || pyModBonds != python::object()) {
    boost::dynamic_bitset<> modifiedAtoms(mol.getNumAtoms());
    boost::dynamic_bitset<> modifiedBonds(mol.getNumBonds());
    auto res = self.enumerate(mol, &modifiedAtoms, &modifiedBonds);
    if (pyModAtoms != python::object()) {
      python::list atList = python::extract<python::list>(pyModAtoms);
      for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        if (modifiedAtoms[i]) {
          atList.append(i);
        }
      }
    }
    if (pyModBonds != python::object()) {
      python::list bndList = python::extract<python::list>(pyModBonds);
      for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
        if (modifiedBonds[i]) {
          bndList.append(i);
        }
      }
    }
    return res;
  } else {
    return self.enumerate(mol);
  }
}

}  // namespace

struct tautomer_wrapper {
  static void wrap() {
    python::class_<MolStandardize::TautomerEnumerator, boost::noncopyable>(
        "TautomerEnumerator", python::no_init)
        .def("__init__", python::make_constructor(createDefaultEnumerator))
        .def("__init__", python::make_constructor(EnumeratorFromParams))
        .def("Enumerate", &enumerateHelper,
             (python::arg("self"), python::arg("mol"),
              python::arg("modifiedAtoms") = python::object(),
              python::arg("modifiedBonds") = python::object()),
             R"DOC(Generates the tautomers for a molecule.
             
  The enumeration rules are inspired by the publication:
  M. Sitzmann et al., “Tautomerism in Large Databases.”, JCAMD 24:521 (2010)
  https://doi.org/10.1007/s10822-010-9346-4
  
  Note: the definitions used here are that the atoms modified during
  tautomerization are the atoms at the beginning and end of each tautomer
  transform (the H "donor" and H "acceptor" in the transform) and the bonds
  modified during transformation are any bonds whose order is changed during
  the tautomer transform (these are the bonds between the "donor" and the
  "acceptor"))DOC")
        .def("Canonicalize", &canonicalizeHelper,
             (python::arg("self"), python::arg("mol")),
             R"DOC(Returns the canonical tautomer for a molecule.

  The default scoring scheme is inspired by the publication:
  M. Sitzmann et al., “Tautomerism in Large Databases.”, JCAMD 24:521 (2010)
  https://doi.org/10.1007/s10822-010-9346-4

  Note that the canonical tautomer is very likely not the most stable tautomer
  for any given conditions. The default scoring rules are designed to produce
  "reasonable" tautomers, but the primary concern is that the results are
  canonical: you always get the same canonical tautomer for a molecule
  regardless of what the input tautomer or atom ordering were.)DOC",
             python::return_value_policy<python::manage_new_object>())
        .def(
            "Canonicalize", &canonicalizeHelper2,
            (python::arg("self"), python::arg("mol"), python::arg("scoreFunc")),
            "returns the canonical tautomer for a molecule using a custom "
            "scoring function",
            python::return_value_policy<python::manage_new_object>())
        .def("ScoreTautomer", &scoreHelper,
             (python::arg("self"), python::arg("mol")),
             "returns the score for a tautomer using the default scoring "
             "scheme.")
        .def_readonly(
            "tautomerScoreVersion",
            MolStandardize::TautomerScoringFunctions::tautomerScoringVersion);
  }
};

void wrap_tautomer() { tautomer_wrapper::wrap(); }
