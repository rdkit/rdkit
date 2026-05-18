//
//  Copyright (C) 2006-2026 Rational Discovery LLC and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>

#include <boost/dynamic_bitset.hpp>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureFactory.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

namespace {
nb::object GetAtomMatch(nb::object featMatch, int maxAts = 1024) {
  nb::list res;
  auto nEntries = nb::len(featMatch);

  boost::dynamic_bitset<> indices(maxAts);
  for (size_t i = 0; i < static_cast<size_t>(nEntries); ++i) {
    auto feat = nb::cast<MolChemicalFeature *>(featMatch[i]);
    nb::list local;
    for (const auto *atom : feat->getAtoms()) {
      unsigned int idx = atom->getIdx();
      if (indices[idx]) {
        return nb::list();
      } else {
        indices[idx] = 1;
      }
      local.append(idx);
    }
    res.append(local);
  }
  return res;
}
}  // namespace

void wrap_ChemicalFeatureUtils(nb::module_ &m) {
  m.def("GetAtomMatch", GetAtomMatch, "featMatch"_a, "maxAts"_a = 1024,
        R"DOC(Returns an empty list if any of the features passed in share an atom.
Otherwise a list of lists of atom indices is returned.
)DOC");
}
