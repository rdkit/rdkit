//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <string>
#include <tuple>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

// ours
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolBundle.h>
#include <RDBoost/Wrap_nb.h>

#include "substructmethods.h"
namespace nb = nanobind;
using namespace nb::literals;

#ifdef RDK_USE_BOOST_SERIALIZATION
nb::bytes BundleToBinary(const RDKit::MolBundle &self) {
  auto res = self.serialize();
  return nb::bytes(res.c_str(), res.length());
}

#else
nb::bytes BundleToBinary(const RDKit::MolBundle &) {
  throw ValueErrorException("Pickling of MolBundle instances is not enabled");
}
#endif

namespace RDKit {

namespace {
const ROMol *getMolHelper(const MolBundle &self, unsigned int idx) {
  return self.getMol(idx).get();
}

size_t addMolHelper(MolBundle &self, const ROMol &nmol) {
  return self.addMol(ROMOL_SPTR(new ROMol(nmol)));
}
}  // namespace

std::string molBundleClassDoc =
    R"DOC(A class for storing groups of related molecules.

)DOC";
struct molbundle_wrap {
  static void wrap(nb::module_ &m) {
    nb::class_<MolBundle>(m, "MolBundle", molBundleClassDoc.c_str())
        .def(nb::init<>())
        .def(nb::init<const std::string &>(), "pkl"_a)
        .def("ToBinary", BundleToBinary,
             R"DOC(Returns a binary string representation of the MolBundle.
)DOC")

        .def("__getitem__", getMolHelper, nb::rv_policy::reference_internal,
             "idx"_a)
        .def("__len__", &MolBundle::size)
        .def("AddMol", addMolHelper, "nmol"_a)
        .def("GetMol", getMolHelper, nb::rv_policy::reference_internal, "idx"_a)
        .def("Size", &MolBundle::size)

        // substructures
        .def(
            "HasSubstructMatch",
            (bool (*)(const MolBundle &m, const ROMol &query, bool, bool,
                      bool))HasSubstructMatch,
            "query"_a, "recursionPossible"_a = true, "useChirality"_a = false,
            "useQueryQueryMatches"_a = false,
            R"DOC(Queries whether or not any molecule in the bundle contains a particular substructure.

     ARGUMENTS:
          - query: a Molecule

          - recursionPossible: (optional)

          - useChirality: enables the use of stereochemistry in the matching

          - useQueryQueryMatches: use query-query matching logic

     RETURNS: True or False
)DOC")
        .def(
            "GetSubstructMatch",
            (std::vector<int> (*)(const MolBundle &m, const ROMol &query, bool,
                                  bool))GetSubstructMatch,
            "query"_a, "useChirality"_a = false,
            "useQueryQueryMatches"_a = false,
            R"DOC(Returns the indices of the atoms from the first molecule in a bundle that matches a substructure query.

     ARGUMENTS:
          - query: a Molecule

          - useChirality: enables the use of stereochemistry in the matching

          - useQueryQueryMatches: use query-query matching logic

     RETURNS: a tuple of integers

     NOTES:
           - only a single match is returned
           - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
)DOC")
        .def(
            "GetSubstructMatches",
            (std::vector<std::vector<int>> (*)(
                const MolBundle &m, const ROMol &query, bool, bool, bool,
                unsigned int))GetSubstructMatches,
            "query"_a, "uniquify"_a = true, "useChirality"_a = false,
            "useQueryQueryMatches"_a = false, "maxMatches"_a = 1000,
            R"DOC(Returns tuple of all indices of the atoms from the first molecule in a bundle that matches a substructure query.

     ARGUMENTS:
          - query: a molecule.
          - uniquify: (optional) determines whether or not the matches are uniquified.
                                        Defaults to 1.

          - useChirality: enables the use of stereochemistry in the matching

          - useQueryQueryMatches: use query-query matching logic

          - maxMatches: The maximum number of matches that will be returned.
                                             In high-symmetry cases with medium-sized molecules, it is
                                             very easy to end up with a combinatorial explosion in the
                                             number of possible matches. This argument prevents that from
                                             having unintended consequences

     RETURNS: a tuple of tuples of integers

     NOTE:
           - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
)DOC")
        .def(
            "HasSubstructMatch",
            (bool (*)(const MolBundle &m, const MolBundle &query, bool, bool,
                      bool))HasSubstructMatch,
            "query"_a, "recursionPossible"_a = true, "useChirality"_a = false,
            "useQueryQueryMatches"_a = false,
            R"DOC(Queries whether or not any molecule in the first bundle matches any molecule in the second bundle.

     ARGUMENTS:
          - query: a MolBundle

          - recursionPossible: (optional)

          - useChirality: enables the use of stereochemistry in the matching

          - useQueryQueryMatches: use query-query matching logic

     RETURNS: True or False
)DOC")
        .def(
            "GetSubstructMatch",
            (std::vector<int> (*)(const MolBundle &m, const MolBundle &query,
                                  bool, bool))GetSubstructMatch,
            "query"_a, "useChirality"_a = false,
            "useQueryQueryMatches"_a = false,
            R"DOC(Returns the indices of the atoms from the first molecule in a bundle that matches a substructure query from a bundle.

     ARGUMENTS:
          - query: a MolBundle

          - useChirality: enables the use of stereochemistry in the matching

          - useQueryQueryMatches: use query-query matching logic

     RETURNS: a tuple of integers

     NOTES:
           - only a single match is returned
           - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
)DOC")

        .def(
            "GetSubstructMatches",
            (std::vector<std::vector<int>> (*)(
                const MolBundle &m, const MolBundle &query, bool, bool, bool,
                unsigned int))GetSubstructMatches,
            "query"_a, "uniquify"_a = true, "useChirality"_a = false,
            "useQueryQueryMatches"_a = false, "maxMatches"_a = 1000,
            R"DOC(Returns tuple of all indices of the atoms from the first molecule in a bundle that matches a substructure query from the second bundle.

     ARGUMENTS:
          - query: a MolBundle.
          - uniquify: (optional) determines whether or not the matches are uniquified.
                                        Defaults to 1.

          - useChirality: enables the use of stereochemistry in the matching

          - useQueryQueryMatches: use query-query matching logic

          - maxMatches: The maximum number of matches that will be returned.
                                             In high-symmetry cases with medium-sized molecules, it is
                                             very easy to end up with a combinatorial explosion in the
                                             number of possible matches. This argument prevents that from
                                             having unintended consequences

     RETURNS: a tuple of tuples of integers

     NOTE:
           - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
)DOC")

        // ------------------------------------------------
        .def(
            "HasSubstructMatch",
            (bool (*)(const MolBundle &m, const ROMol &query,
                      const std::optional<SubstructMatchParameters>))
                helpHasSubstructMatch,
            "query"_a, "params"_a,
            R"DOC(Queries whether or not any molecule in the bundle contains a particular substructure.

     ARGUMENTS:
          - query: a Molecule

          - params: parameters controlling the substructure match

     RETURNS: True or False
)DOC")
        .def(
            "GetSubstructMatch",
            (std::vector<int> (*)(
                const MolBundle &m, const ROMol &query,
                const std::optional<SubstructMatchParameters>))
                helpGetSubstructMatch,
            "query"_a, "params"_a,
            R"DOC(Returns the indices of the atoms from the first molecule in a bundle that matches a substructure query.

     ARGUMENTS:
          - query: a Molecule

          - params: parameters controlling the substructure match

     RETURNS: a tuple of integers

     NOTES:
           - only a single match is returned
           - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
)DOC")
        .def(
            "GetSubstructMatches",
            (std::vector<std::vector<int>> (*)(
                const MolBundle &m, const ROMol &query,
                const std::optional<SubstructMatchParameters>))
                helpGetSubstructMatches,
            "query"_a, "params"_a,
            R"DOC(Returns tuple of all indices of the atoms from the first molecule in a bundle that matches a substructure query.

     ARGUMENTS:
          - query: a molecule.
          - params: parameters controlling the substructure match

     RETURNS: a tuple of tuples of integers

     NOTE:
           - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
)DOC")
        .def(
            "HasSubstructMatch",
            (bool (*)(const MolBundle &m, const MolBundle &query,
                      const std::optional<SubstructMatchParameters>))
                helpHasSubstructMatch,
            "query"_a, "params"_a,
            R"DOC(Queries whether or not any molecule in the first bundle matches any molecule in the second bundle.

     ARGUMENTS:
          - query: a MolBundle

          - params: parameters controlling the substructure match

     RETURNS: True or False
)DOC")
        .def(
            "GetSubstructMatch",
            (std::vector<int> (*)(
                const MolBundle &m, const MolBundle &query,
                const std::optional<SubstructMatchParameters>))
                helpGetSubstructMatch,
            "query"_a, "params"_a,
            R"DOC(Returns the indices of the atoms from the first molecule in a bundle that matches a substructure query from a bundle.

     ARGUMENTS:
          - query: a MolBundle

          - params: parameters controlling the substructure match

     RETURNS: a tuple of integers

     NOTES:
           - only a single match is returned
           - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
)DOC")

        .def(
            "GetSubstructMatches",
            (std::vector<std::vector<int>> (*)(
                const MolBundle &m, const MolBundle &query,
                const std::optional<SubstructMatchParameters>))
                helpGetSubstructMatches,
            "query"_a, "params"_a,
            R"DOC(Returns tuple of all indices of the atoms from the first molecule in a bundle that matches a substructure query from the second bundle.

     ARGUMENTS:
          - query: a MolBundle.
          - params: parameters controlling the substructure match

     RETURNS: a tuple of tuples of integers

     NOTE:
           - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
)DOC")

        .def("__getstate__",
             [](const MolBundle &bundle) {
               const auto pkl = BundleToBinary(bundle);
               return std::make_tuple(pkl);
             })
        .def("__setstate__",
             [](MolBundle &bundle, const std::tuple<nb::bytes> &state) {
               std::string pkl = std::string(
                   static_cast<const char *>(std::get<0>(state).data()),
                   static_cast<size_t>(std::get<0>(state).size()));
               new (&bundle) MolBundle(pkl);
             })
        .doc() = molBundleClassDoc.c_str();

    molBundleClassDoc =
        R"DOC(A class for storing groups of related molecules.
          Here related means that the molecules have to have the same number of atoms.

)DOC";
    nb::class_<FixedMolSizeMolBundle, MolBundle>(m, "FixedMolSizeMolBundle",
                                                 molBundleClassDoc.c_str())
        .def(nb::init<>());

    m.def(
        "MolBundleCanSerialize", MolBundleCanSerialize,
        R"DOC(Returns True if the MolBundle is serializable (requires boost serialization))DOC");
  };
};
}  // namespace RDKit

void wrap_molbundle(nb::module_ &m) { RDKit::molbundle_wrap::wrap(m); }
