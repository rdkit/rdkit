//
//  Copyright (C) 2017 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define NO_IMPORT_ARRAY
#include <boost/python.hpp>
#include <string>

// ours
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolBundle.h>
#include <RDBoost/PySequenceHolder.h>

#include "substructmethods.h"
namespace python = boost::python;

#ifdef RDK_USE_BOOST_SERIALIZATION
struct molbundle_pickle_suite : rdkit_pickle_suite {
  static python::tuple getinitargs(const RDKit::MolBundle &self) {
    auto res = self.serialize();
    return python::make_tuple(python::object(python::handle<>(
        PyBytes_FromStringAndSize(res.c_str(), res.length()))));
  };
};

python::object BundleToBinary(const RDKit::MolBundle &self) {
  auto res = self.serialize();
  python::object retval = python::object(
      python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
  return retval;
}

#else
struct molbundle_pickle_suite : rdkit_pickle_suite {
  static python::tuple getinitargs(const RDKit::MolBundle &) {
    throw_runtime_error("Pickling of MolBundle instances is not enabled");
    return python::tuple();  // warning suppression, we never get here
  };
};

python::object BundleToBinary(const RDKit::MolBundle &) {
  throw_runtime_error("Pickling of MolBundle instances is not enabled");
  return python::object();  // warning suppression, we never get here
}
#endif

namespace RDKit {

std::string molBundleClassDoc =
    "A class for storing groups of related molecules.\n\
\n";
struct molbundle_wrap {
  static void wrap() {
    python::class_<MolBundle, boost::noncopyable>(
        "MolBundle", molBundleClassDoc.c_str(),
        python::init<>(python::args("self")))
        .def(python::init<const std::string &>(python::args("self", "pkl")))
        .def_pickle(molbundle_pickle_suite())
        .def("ToBinary", BundleToBinary, python::args("self"),
             "Returns a binary string representation of the MolBundle.\n")

        .def("__getitem__", &MolBundle::getMol, python::args("self", "idx"))
        .def("__len__", &MolBundle::size, python::args("self"))
        .def("AddMol", &MolBundle::addMol, python::args("self", "nmol"))
        .def("GetMol", &MolBundle::getMol, python::args("self", "idx"))
        .def("Size", &MolBundle::size, python::args("self"))

        // substructures
        .def("HasSubstructMatch",
             (bool (*)(const MolBundle &m, const ROMol &query, bool, bool,
                       bool))HasSubstructMatch,
             (python::arg("self"), python::arg("query"),
              python::arg("recursionPossible") = true,
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false),
             "Queries whether or not any molecule in the bundle contains a "
             "particular substructure.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a Molecule\n\n"
             "    - recursionPossible: (optional)\n\n"
             "    - useChirality: enables the use of stereochemistry in the "
             "matching\n\n"
             "    - useQueryQueryMatches: use query-query matching logic\n\n"
             "  RETURNS: True or False\n")
        .def("GetSubstructMatch",
             (PyObject * (*)(const MolBundle &m, const ROMol &query, bool,
                             bool)) GetSubstructMatch,
             (python::arg("self"), python::arg("query"),
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false),
             "Returns the indices of the atoms from the first molecule in a "
             "bundle that matches a substructure query.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a Molecule\n\n"
             "    - useChirality: enables the use of stereochemistry in the "
             "matching\n\n"
             "    - useQueryQueryMatches: use query-query matching logic\n\n"
             "  RETURNS: a tuple of integers\n\n"
             "  NOTES:\n"
             "     - only a single match is returned\n"
             "     - the ordering of the indices corresponds to the atom "
             "ordering\n"
             "         in the query. For example, the first index is for the "
             "atom in\n"
             "         this molecule that matches the first atom in the "
             "query.\n")
        .def("GetSubstructMatches",
             (PyObject * (*)(const MolBundle &m, const ROMol &query, bool, bool,
                             bool, unsigned int)) GetSubstructMatches,
             (python::arg("self"), python::arg("query"),
              python::arg("uniquify") = true,
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false,
              python::arg("maxMatches") = 1000),
             "Returns tuple of all indices of the atoms from the first "
             "molecule in a bundle that matches a substructure query.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a molecule.\n"
             "    - uniquify: (optional) determines whether or not the "
             "matches "
             "are uniquified.\n"
             "                Defaults to 1.\n\n"
             "    - useChirality: enables the use of stereochemistry in the "
             "matching\n\n"
             "    - useQueryQueryMatches: use query-query matching logic\n\n"
             "    - maxMatches: The maximum number of matches that will be "
             "returned.\n"
             "                  In high-symmetry cases with medium-sized "
             "molecules, it is\n"
             "                  very easy to end up with a combinatorial "
             "explosion in the\n"
             "                  number of possible matches. This argument "
             "prevents that from\n"
             "                  having unintended consequences\n\n"
             "  RETURNS: a tuple of tuples of integers\n\n"
             "  NOTE:\n"
             "     - the ordering of the indices corresponds to the atom "
             "ordering\n"
             "         in the query. For example, the first index is for the "
             "atom in\n"
             "         this molecule that matches the first atom in the "
             "query.\n")
        .def("HasSubstructMatch",
             (bool (*)(const MolBundle &m, const MolBundle &query, bool, bool,
                       bool))HasSubstructMatch,
             (python::arg("self"), python::arg("query"),
              python::arg("recursionPossible") = true,
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false),
             "Queries whether or not any molecule in the first bundle matches "
             "any molecule in the second bundle.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a MolBundle\n\n"
             "    - recursionPossible: (optional)\n\n"
             "    - useChirality: enables the use of stereochemistry in the "
             "matching\n\n"
             "    - useQueryQueryMatches: use query-query matching logic\n\n"
             "  RETURNS: True or False\n")
        .def("GetSubstructMatch",
             (PyObject * (*)(const MolBundle &m, const MolBundle &query, bool,
                             bool)) GetSubstructMatch,
             (python::arg("self"), python::arg("query"),
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false),
             "Returns the indices of the atoms from the first molecule in a "
             "bundle that matches a substructure query from a bundle.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a MolBundle\n\n"
             "    - useChirality: enables the use of stereochemistry in the "
             "matching\n\n"
             "    - useQueryQueryMatches: use query-query matching logic\n\n"
             "  RETURNS: a tuple of integers\n\n"
             "  NOTES:\n"
             "     - only a single match is returned\n"
             "     - the ordering of the indices corresponds to the atom "
             "ordering\n"
             "         in the query. For example, the first index is for the "
             "atom in\n"
             "         this molecule that matches the first atom in the "
             "query.\n")

        .def("GetSubstructMatches",
             (PyObject * (*)(const MolBundle &m, const MolBundle &query, bool,
                             bool, bool, unsigned int)) GetSubstructMatches,
             (python::arg("self"), python::arg("query"),
              python::arg("uniquify") = true,
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false,
              python::arg("maxMatches") = 1000),
             "Returns tuple of all indices of the atoms from the first "
             "molecule in a bundle that matches a substructure query from the "
             "second bundle.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a MolBundle.\n"
             "    - uniquify: (optional) determines whether or not the "
             "matches "
             "are uniquified.\n"
             "                Defaults to 1.\n\n"
             "    - useChirality: enables the use of stereochemistry in the "
             "matching\n\n"
             "    - useQueryQueryMatches: use query-query matching logic\n\n"
             "    - maxMatches: The maximum number of matches that will be "
             "returned.\n"
             "                  In high-symmetry cases with medium-sized "
             "molecules, it is\n"
             "                  very easy to end up with a combinatorial "
             "explosion in the\n"
             "                  number of possible matches. This argument "
             "prevents that from\n"
             "                  having unintended consequences\n\n"
             "  RETURNS: a tuple of tuples of integers\n\n"
             "  NOTE:\n"
             "     - the ordering of the indices corresponds to the atom "
             "ordering\n"
             "         in the query. For example, the first index is for the "
             "atom in\n"
             "         this molecule that matches the first atom in the "
             "query.\n")

        // ------------------------------------------------
        .def("HasSubstructMatch",
             (bool (*)(const MolBundle &m, const ROMol &query,
                       const SubstructMatchParameters &))helpHasSubstructMatch,
             (python::arg("self"), python::arg("query"), python::arg("params")),
             "Queries whether or not any molecule in the bundle contains a "
             "particular substructure.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a Molecule\n\n"
             "    - params: parameters controlling the substructure match\n\n"
             "matching\n\n"
             "    - useQueryQueryMatches: use query-query matching logic\n\n"
             "  RETURNS: True or False\n")
        .def("GetSubstructMatch",
             (PyObject * (*)(const MolBundle &m, const ROMol &query,
                             const SubstructMatchParameters &))
                 helpGetSubstructMatch,
             (python::arg("self"), python::arg("query"), python::arg("params")),
             "Returns the indices of the atoms from the first molecule in a "
             "bundle that matches a substructure query.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a Molecule\n\n"
             "    - params: parameters controlling the substructure match\n\n"
             "  RETURNS: a tuple of integers\n\n"
             "  NOTES:\n"
             "     - only a single match is returned\n"
             "     - the ordering of the indices corresponds to the atom "
             "ordering\n"
             "         in the query. For example, the first index is for the "
             "atom in\n"
             "         this molecule that matches the first atom in the "
             "query.\n")
        .def("GetSubstructMatches",
             (PyObject * (*)(const MolBundle &m, const ROMol &query,
                             const SubstructMatchParameters &))
                 helpGetSubstructMatches,
             (python::arg("self"), python::arg("query"), python::arg("params")),
             "Returns tuple of all indices of the atoms from the first "
             "molecule in a bundle that matches a substructure query.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a molecule.\n"
             "    - params: parameters controlling the substructure match\n\n"
             "  RETURNS: a tuple of tuples of integers\n\n"
             "  NOTE:\n"
             "     - the ordering of the indices corresponds to the atom "
             "ordering\n"
             "         in the query. For example, the first index is for the "
             "atom in\n"
             "         this molecule that matches the first atom in the "
             "query.\n")
        .def("HasSubstructMatch",
             (bool (*)(const MolBundle &m, const MolBundle &query,
                       const SubstructMatchParameters &))helpHasSubstructMatch,
             (python::arg("self"), python::arg("query"), python::arg("params")),
             "Queries whether or not any molecule in the first bundle matches "
             "any molecule in the second bundle.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a MolBundle\n\n"
             "    - params: parameters controlling the substructure match\n\n"
             "  RETURNS: True or False\n")
        .def("GetSubstructMatch",
             (PyObject * (*)(const MolBundle &m, const MolBundle &query,
                             const SubstructMatchParameters &))
                 helpGetSubstructMatch,
             (python::arg("self"), python::arg("query"), python::arg("params")),
             "Returns the indices of the atoms from the first molecule in a "
             "bundle that matches a substructure query from a bundle.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a MolBundle\n\n"
             "    - params: parameters controlling the substructure match\n\n"
             "  RETURNS: a tuple of integers\n\n"
             "  NOTES:\n"
             "     - only a single match is returned\n"
             "     - the ordering of the indices corresponds to the atom "
             "ordering\n"
             "         in the query. For example, the first index is for the "
             "atom in\n"
             "         this molecule that matches the first atom in the "
             "query.\n")

        .def("GetSubstructMatches",
             (PyObject * (*)(const MolBundle &m, const MolBundle &query,
                             const SubstructMatchParameters &))
                 helpGetSubstructMatches,
             (python::arg("self"), python::arg("query"), python::arg("params")),
             "Returns tuple of all indices of the atoms from the first "
             "molecule in a bundle that matches a substructure query from the "
             "second bundle.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a MolBundle.\n"
             "    - params: parameters controlling the substructure match\n\n"
             "  RETURNS: a tuple of tuples of integers\n\n"
             "  NOTE:\n"
             "     - the ordering of the indices corresponds to the atom "
             "ordering\n"
             "         in the query. For example, the first index is for the "
             "atom in\n"
             "         this molecule that matches the first atom in the "
             "query.\n");
    molBundleClassDoc =
        "A class for storing groups of related molecules.\n\
    Here related means that the molecules have to have the same number of atoms.\n\
\n";
    python::class_<FixedMolSizeMolBundle, python::bases<MolBundle>>(
        "FixedMolSizeMolBundle", molBundleClassDoc.c_str(),
        python::init<>(python::args("self")));

    python::def("MolBundleCanSerialize", MolBundleCanSerialize,
                "Returns True if the MolBundle is serializable "
                "(requires boost serialization");
  };
};
}  // namespace RDKit

void wrap_molbundle() { RDKit::molbundle_wrap::wrap(); }
