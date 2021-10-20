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

class PyResonanceMolSupplierCallback
    : public ResonanceMolSupplierCallback,
      public python::wrapper<ResonanceMolSupplierCallback> {
 public:
  PyResonanceMolSupplierCallback() {}
  PyResonanceMolSupplierCallback(const python::object &pyCallbackObject) {
    PyResonanceMolSupplierCallback *pyCallback =
        python::extract<PyResonanceMolSupplierCallback *>(pyCallbackObject);
    *this = *pyCallback;
    d_pyCallbackObject = pyCallbackObject;
    pyCallback->d_cppCallback = this;
  }
  inline unsigned int wrapGetNumConjGrps() const {
    return d_cppCallback->getNumConjGrps();
  }
  inline size_t wrapGetMaxStructures() const {
    return d_cppCallback->getMaxStructures();
  }
  inline size_t wrapGetNumStructures(unsigned int conjGrpIdx) const {
    return d_cppCallback->getNumStructures(conjGrpIdx);
  }
  inline size_t wrapGetNumDiverseStructures(unsigned int conjGrpIdx) const {
    return d_cppCallback->getNumDiverseStructures(conjGrpIdx);
  }
  inline python::object getCallbackOverride() const {
    return get_override("__call__");
  }
  bool operator()() override { return getCallbackOverride()(); }
  python::object getPyCallbackObject() { return d_pyCallbackObject; }

 private:
  PyResonanceMolSupplierCallback *d_cppCallback;
  python::object d_pyCallbackObject;
};

python::object getProgressCallbackHelper(const ResonanceMolSupplier &suppl) {
  PyResonanceMolSupplierCallback *cppCallback =
      dynamic_cast<PyResonanceMolSupplierCallback *>(
          suppl.getProgressCallback());
  python::object res;
  if (cppCallback) {
    res = cppCallback->getPyCallbackObject();
  }
  return res;
};

void setProgressCallbackHelper(ResonanceMolSupplier &suppl,
                               PyObject *callback) {
  PRECONDITION(callback, "callback must not be NULL");
  if (callback == Py_None) {
    suppl.setProgressCallback(nullptr);
    return;
  }
  python::object callbackObject(python::handle<>(python::borrowed(callback)));
  python::extract<PyResonanceMolSupplierCallback *> extractCallback(
      callbackObject);
  if (extractCallback.check()) {
    if (!PyCallable_Check(extractCallback()->getCallbackOverride().ptr())) {
      PyErr_SetString(PyExc_AttributeError,
                      "The __call__ attribute in the "
                      "rdchem.ResonanceMolSupplierCallback subclass "
                      "must exist and be a callable method");
      python::throw_error_already_set();
    } else {
      suppl.setProgressCallback(
          new PyResonanceMolSupplierCallback(callbackObject));
    }
  } else {
    PyErr_SetString(PyExc_TypeError,
                    "Expected an instance of a "
                    "rdchem.ResonanceMolSupplierCallback subclass");
    python::throw_error_already_set();
  }
}

std::string resonanceMolSupplierCallbackClassDoc =
    "Create a derived class from this abstract base class and\n\
    implement the __call__() method.\n\
    The __call__() method is called at each iteration of the\n\
    algorithm, and provides a mechanism to monitor or stop\n\
    its progress.\n\n\
    To have your callback called, pass an instance of your\n\
    derived class to ResonanceMolSupplier.SetProgressCallback()\n";

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

    python::class_<PyResonanceMolSupplierCallback, boost::noncopyable>(
        "ResonanceMolSupplierCallback",
        resonanceMolSupplierCallbackClassDoc.c_str(), python::init<>())
        .def("GetNumConjGrps",
             &PyResonanceMolSupplierCallback::wrapGetNumConjGrps,
             "Returns the number of individual conjugated groups in the "
             "molecule.\n")
        .def("GetMaxStructures",
             &PyResonanceMolSupplierCallback::wrapGetMaxStructures,
             "Get the number of conjugated groups this molecule has.\n")
        .def("GetNumStructures",
             (size_t(PyResonanceMolSupplierCallback::*)(unsigned int)) &
                 PyResonanceMolSupplierCallback::wrapGetNumStructures,
             "Get the number of resonance structures generated so far "
             "for the passed conjugated group index.\n")
        .def("GetNumDiverseStructures",
             (size_t(PyResonanceMolSupplierCallback::*)(unsigned int)) &
                 PyResonanceMolSupplierCallback::wrapGetNumDiverseStructures,
             "Get the number of non-degenrate resonance structures "
             "generated so far for the passed conjugated group index.\n")
        .def("__call__",
             python::pure_virtual(&PyResonanceMolSupplierCallback::operator()),
             "This must be implemented in the derived class. "
             "Return True if the resonance structure generation "
             "should continue; False if the resonance structure "
             "generation should stop.\n");

    python::class_<ResonanceMolSupplier, boost::noncopyable>(
        "ResonanceMolSupplier", resonanceMolSupplierClassDoc.c_str(),
        python::init<ROMol &, unsigned int, unsigned int>(
            (python::arg("mol"), python::arg("flags") = 0,
             python::arg("maxStructs") = 1000)))
        .def(
            "__iter__",
            (ResonanceMolSupplier * (*)(ResonanceMolSupplier *)) & MolSupplIter,
            python::return_internal_reference<1>())
        .def("__next__", (ROMol * (*)(ResonanceMolSupplier *)) & MolSupplNext,
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
             "molecule.\n")
        .def("GetBondConjGrpIdx",
             (unsigned int (ResonanceMolSupplier::*)(unsigned int)) &
                 ResonanceMolSupplier::getBondConjGrpIdx,
             "Given a bond index, it returns the index of the conjugated group"
             "the bond belongs to, or -1 if it is not conjugated.\n")
        .def("GetAtomConjGrpIdx",
             (unsigned int (ResonanceMolSupplier::*)(unsigned int)) &
                 ResonanceMolSupplier::getAtomConjGrpIdx,
             "Given an atom index, it returns the index of the conjugated group"
             "the atom belongs to, or -1 if it is not conjugated.\n")
        .def(
            "SetNumThreads",
            (void (ResonanceMolSupplier::*)(unsigned int)) &
                ResonanceMolSupplier::setNumThreads,
            "Sets the number of threads to be used to enumerate resonance\n"
            "structures (defaults to 1; 0 selects the number of concurrent\n"
            "threads supported by the hardware; negative values are added\n"
            "to the number of concurrent threads supported by the hardware).\n")
        .def("SetProgressCallback", &setProgressCallbackHelper,
             "Pass an instance of a class derived from\n"
             "ResonanceMolSupplierCallback, which must implement the\n"
             "__call__() method.\n")
        .def("GetProgressCallback", &getProgressCallbackHelper,
             "Get the ResonanceMolSupplierCallback subclass instance,\n"
             "or None if none was set.\n")
        .def("WasCanceled", &ResonanceMolSupplier::wasCanceled,
             "Returns True if the resonance structure generation was "
             "canceled.\n")
        .def("Enumerate", &ResonanceMolSupplier::enumerate,
             "Ask ResonanceMolSupplier to enumerate resonance structures"
             "(automatically done as soon as any attempt to access them is "
             "made).\n")
        .def("GetIsEnumerated", &ResonanceMolSupplier::getIsEnumerated,
             "Returns true if resonance structure enumeration has already "
             "happened.\n")
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
