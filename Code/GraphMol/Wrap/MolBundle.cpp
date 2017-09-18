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
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDBoost/PySequenceHolder.h>
#include <RDBoost/iterator_next.h>

namespace python = boost::python;

namespace RDKit {
PyObject *convertMatches(MatchVectType &matches);

template <typename T1, T2>
bool HasSubstructMatch(const T1 &mol, const T2 &query,
                       bool recursionPossible = true, bool useChirality = false,
                       bool useQueryQueryMatches = false) {
  NOGIL gil;
  MatchVectType res;
  return SubstructMatch(mol, query, res, recursionPossible, useChirality,
                        useQueryQueryMatches);
}

PyObject *GetSubstructMatch(const ROMol &mol, const ROMol &query,
                            bool useChirality = false,
                            bool useQueryQueryMatches = false) {
  MatchVectType matches;
  {
    NOGIL gil;
    SubstructMatch(mol, query, matches, true, useChirality,
                   useQueryQueryMatches);
  }
  return convertMatches(matches);
}

PyObject *GetSubstructMatches(const ROMol &mol, const ROMol &query,
                              bool uniquify = true, bool useChirality = false,
                              bool useQueryQueryMatches = false,
                              unsigned int maxMatches = 1000) {
  std::vector<MatchVectType> matches;
  int matched;
  {
    NOGIL gil;
    matched = SubstructMatch(mol, query, matches, uniquify, true, useChirality,
                             useQueryQueryMatches, maxMatches);
  }
  PyObject *res = PyTuple_New(matched);
  for (int idx = 0; idx < matched; idx++) {
    PyTuple_SetItem(res, idx, convertMatches(matches[idx]));
  }
  return res;
}
std::string molBundleClassDoc =
    "A class for storing gropus of related molecules.\n\
    Here related means that the molecules have to have the same number of atoms.\n\
\n";
struct molbundle_wrap {
  static void wrap() {
    python::class_<MolBundle, boost::noncopyable>(
        "MolBundle", molBundleClassDoc.c_str(), python::init<>())
        .def("__getitem__", &MolBundle::getMol)
        .def("__len__", &MolBundle::size)
        .def("addMol", &MolBundle::addMol)
        .def("getMol", &MolBundle::getMol)
        .def("size", &MolBundle::size)
#if 0
        .def("reset", &ResonanceMolSupplier::reset,
             "Resets our position in the resonance structure supplier to the "
             "beginning.\n")
        .def("__len__", &ResonanceMolSupplier::length)
        .def("atEnd", &ResonanceMolSupplier::atEnd,
             "Returns whether or not we have hit the end of the resonance "
             "structure supplier.\n")
        .def("GetNumConjGrps", &ResonanceMolSupplier::getNumConjGrps,
             "Returns the number of individual conjugated groups in the "
             "molecule\n")
        .def("GetBondConjGrpIdx",
             (unsigned int (ResonanceMolSupplier::*)(unsigned int)) &
                 ResonanceMolSupplier::getBondConjGrpIdx,
             "Given a bond index, it returns the index of the conjugated group"
             "the bond belongs to, or -1 if it is not conjugated\n")
        .def("GetAtomConjGrpIdx",
             (unsigned int (ResonanceMolSupplier::*)(unsigned int)) &
                 ResonanceMolSupplier::getAtomConjGrpIdx,
             "Given an atom index, it returns the index of the conjugated group"
             "the atom belongs to, or -1 if it is not conjugated\n")
        .def("SetNumThreads",
             (void (ResonanceMolSupplier::*)(unsigned int)) &
                 ResonanceMolSupplier::setNumThreads,
             "Sets the number of threads to be used to enumerate resonance\n"
             "structures (defaults to 1; 0 selects the number of concurrent\n"
             "threads supported by the hardware; negative values are added\n"
             "to the number of concurrent threads supported by the hardware)\n")
        .def("Enumerate", &ResonanceMolSupplier::enumerate,
             "Ask ResonanceMolSupplier to enumerate resonance structures"
             "(automatically done as soon as any attempt to access them is "
             "made)\n")
        .def("GetIsEnumerated", &ResonanceMolSupplier::getIsEnumerated,
             "Returns true if resonance structure enumeration has already "
             "happened\n")
        .def("GetSubstructMatch", GetResonanceSubstructMatch,
             (python::arg("self"), python::arg("query"),
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false),
             "Returns the indices of the molecule's atoms that match a "
             "substructure query,\n"
             "taking into account all resonance structures in "
             "ResonanceMolSupplier.\n\n"
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
        .def("GetSubstructMatches", GetResonanceSubstructMatches,
             (python::arg("self"), python::arg("query"),
              python::arg("uniquify") = false,
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false,
              python::arg("maxMatches") = 1000, python::arg("numThreads") = 1),
             "Returns tuples of the indices of the molecule's atoms that match "
             "a substructure query,\n"
             "taking into account all resonance structures in "
             "ResonanceMolSupplier.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a Molecule.\n"
             "    - uniquify: (optional) determines whether or not the matches "
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
             "    - numThreads: The number of threads to be used (defaults to "
             "1; 0 selects the\n"
             "                  number of concurrent threads supported by the "
             "hardware; negative\n"
             "                  values are added to the number of concurrent "
             "threads supported\n"
             "                  by the hardware).\n\n"
             "  RETURNS: a tuple of tuples of integers\n\n"
             "  NOTE:\n"
             "     - the ordering of the indices corresponds to the atom "
             "ordering\n"
             "         in the query. For example, the first index is for the "
             "atom in\n"
             "         this molecule that matches the first atom in the "
             "query.\n")
#endif

        // substructures
        .def("HasSubstructMatch", HasSubstructMatch,
             (python::arg("self"), python::arg("query"),
              python::arg("recursionPossible") = true,
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false),
             "Queries whether or not the molecule contains a particular "
             "substructure.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a Molecule\n\n"
             "    - recursionPossible: (optional)\n\n"
             "    - useChirality: enables the use of stereochemistry in the "
             "matching\n\n"
             "    - useQueryQueryMatches: use query-query matching logic\n\n"
             "  RETURNS: True or False\n")
        .def("GetSubstructMatch", GetSubstructMatch,
             (python::arg("self"), python::arg("query"),
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false),
             "Returns the indices of the molecule's atoms that match a "
             "substructure query.\n\n"
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

        .def("GetSubstructMatches", GetSubstructMatches,
             (python::arg("self"), python::arg("query"),
              python::arg("uniquify") = true,
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false,
              python::arg("maxMatches") = 1000),
             "Returns tuples of the indices of the molecule's atoms that "
             "match "
             "a substructure query.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a Molecule.\n"
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

        ;
  };
};
}

void wrap_molbundle() { RDKit::molbundle_wrap::wrap(); }
