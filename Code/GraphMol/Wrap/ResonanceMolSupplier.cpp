// $Id$
//
//  Copyright (C) 2015 Paolo Tosco
//
//  Copyright (C) 2003-2010  Greg Landrum and Rational Discovery LLC
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
#include <GraphMol/Resonance.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDBoost/PySequenceHolder.h>
#include <RDBoost/iterator_next.h>

#include "MolSupplier.h"
#include "substructmethods.h"

namespace python = boost::python;

namespace RDKit {

PyObject *GetResonanceSubstructMatches(
    ResonanceMolSupplier &suppl, const ROMol &query, bool uniquify = false,
    bool useChirality = false, bool useQueryQueryMatches = false,
    unsigned int maxMatches = 1000, int numThreads = 1) {
  std::vector<MatchVectType> matches;
  int matched =
      SubstructMatch(suppl, query, matches, uniquify, true, useChirality,
                     useQueryQueryMatches, maxMatches, numThreads);
  PyObject *res = PyTuple_New(matched);
  for (int idx = 0; idx < matched; idx++) {
    PyTuple_SetItem(res, idx, convertMatches(matches[idx]));
  }
  return res;
}

std::string resonanceMolSupplierClassDoc =
    "A class which supplies resonance structures (as mols) from a mol.\n\
\n\
  Usage examples:\n\
\n\
    1) Lazy evaluation: the resonance structures are not constructed\n\
       until we ask for them:\n\n\
       >>> suppl = ResonanceMolSupplier(mol)\n\
       >>> for resMol in suppl:\n\
       ...    resMol.GetNumAtoms()\n\
\n\
    2) Lazy evaluation 2:\n\n\
       >>> suppl = ResonanceMolSupplier(mol)\n\
       >>> resMol1 = next(suppl)\n\
       >>> resMol2 = next(suppl)\n\
       >>> suppl.reset()\n\
       >>> resMol3 = next(suppl)\n\
       # resMol3 and resMol1 are the same: \n\
       >>> MolToSmiles(resMol3)==MolToSmiles(resMol1)\n\
\n\
    3) Random Access:\n\n\
       >>> suppl = ResonanceMolSupplier(mol)\n\
       >>> resMol1 = suppl[0] \n\
       >>> resMol2 = suppl[1] \n\n\
       NOTE: this will generate an IndexError if the supplier doesn't have that many\n\
       molecules.\n\
\n\
    4) Random Access 2: looping over all resonance structures\n\
       >>> suppl = ResonanceMolSupplier(mol)\n\
       >>> nResMols = len(suppl)\n\
       >>> for i in range(nResMols):\n\
       ...   suppl[i].GetNumAtoms()\n\
\n";
struct resmolsup_wrap {
  static void wrap() {
    python::enum_<ResonanceMolSupplier::ResonanceFlags>("ResonanceFlags")
        .value("ALLOW_INCOMPLETE_OCTETS",
               ResonanceMolSupplier::ALLOW_INCOMPLETE_OCTETS)
        .value("ALLOW_CHARGE_SEPARATION",
               ResonanceMolSupplier::ALLOW_CHARGE_SEPARATION)
        .value("KEKULE_ALL", ResonanceMolSupplier::KEKULE_ALL)
        .value("UNCONSTRAINED_CATIONS",
               ResonanceMolSupplier::UNCONSTRAINED_CATIONS)
        .value("UNCONSTRAINED_ANIONS",
               ResonanceMolSupplier::UNCONSTRAINED_ANIONS)
        .export_values();
    python::class_<ResonanceMolSupplier, boost::noncopyable>(
        "ResonanceMolSupplier", resonanceMolSupplierClassDoc.c_str(),
        python::init<ROMol &, unsigned int, unsigned int>(
            (python::arg("mol"), python::arg("flags") = 0,
             python::arg("maxStructs") = 1000)))
        .def(
            "__iter__",
            (ResonanceMolSupplier * (*)(ResonanceMolSupplier *)) & MolSupplIter,
            python::return_internal_reference<1>())
        .def(NEXT_METHOD,
             (ROMol * (*)(ResonanceMolSupplier *)) &
                 MolSupplNextAcceptNullLastMolecule,
             "Returns the next resonance structure in the supplier. Raises "
             "_StopIteration_ on end.\n",
             python::return_value_policy<python::manage_new_object>())
        .def("__getitem__",
             (ROMol * (*)(ResonanceMolSupplier *, int)) & MolSupplGetItem,
             python::return_value_policy<python::manage_new_object>())
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
        .def("GetSubstructMatch",
             (PyObject * (*)(ResonanceMolSupplier & m, const ROMol &query, bool,
                             bool)) GetSubstructMatch,
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
             "query.\n");
  };
};
}  // namespace RDKit

void wrap_resmolsupplier() { RDKit::resmolsup_wrap::wrap(); }
