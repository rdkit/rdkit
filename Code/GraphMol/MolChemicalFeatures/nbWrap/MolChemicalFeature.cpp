//
//  Copyright (C) 2004-2026 Rational Discovery LLC and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureFactory.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

namespace {
nb::tuple getFeatAtomIds(const MolChemicalFeature &feat) {
  nb::list res;
  for (const auto *atom : feat.getAtoms()) {
    res.append(atom->getIdx());
  }
  return nb::tuple(res);
}
}  // namespace

void wrap_MolChemicalFeat(nb::module_ &m) {
  nb::class_<MolChemicalFeature>(m, "MolChemicalFeature",
      R"DOC(Class to represent a chemical feature.
These chemical features may or may not have been derived from molecule object;
i.e. it is possible to have a chemical feature that was created just from its type
and location.
)DOC")
      .def("GetId", &MolChemicalFeature::getId,
           "Returns the identifier of the feature")
      .def("GetFamily", &MolChemicalFeature::getFamily,
           "Get the family to which the feature belongs; donor, acceptor, etc.")
      .def("GetType", &MolChemicalFeature::getType,
           "Get the specific type for the feature")
      .def("GetPos",
           nb::overload_cast<int>(&MolChemicalFeature::getPos, nb::const_),
           "confId"_a, "Get the location of the chemical feature")
      .def("GetPos",
           nb::overload_cast<>(&MolChemicalFeature::getPos, nb::const_),
           "Get the location of the default chemical feature (first position)")
      .def("GetAtomIds", getFeatAtomIds,
           "Get the IDs of the atoms that participate in the feature")
      .def("GetMol", &MolChemicalFeature::getMol,
           "Get the molecule used to derive the features",
           nb::rv_policy::reference)
      .def("GetFactory", &MolChemicalFeature::getFactory,
           "Get the factory used to generate this feature",
           nb::rv_policy::reference)
      .def("ClearCache", &MolChemicalFeature::clearCache,
           "Clears the cache used to store position information.")
      .def("SetActiveConformer", &MolChemicalFeature::setActiveConformer,
           "confId"_a,
           "Sets the conformer to use (must be associated with a molecule).")
      .def("GetActiveConformer", &MolChemicalFeature::getActiveConformer,
           "Gets the conformer to use.");
}
