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
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

// ours
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Resonance.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/FileParseException.h>
#include <RDBoost/Wrap_nb.h>

#include "substructmethods.h"

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {

namespace {
template <typename T>
T *MolSupplIter(T *suppl) {
  suppl->reset();
  return suppl;
}

template <typename T>
ROMol *MolSupplNext(T *suppl) {
  ROMol *res = nullptr;
  if (!suppl->atEnd()) {
    try {
      res = suppl->next();
    } catch (const FileParseException &) {
      throw;
    } catch (...) {
      res = nullptr;
    }
  } else {
    throw nb::stop_iteration();
  }
  return res;
}

template <typename T>
ROMol *MolSupplGetItem(T *suppl, int idx) {
  ROMol *res = nullptr;
  if (idx < 0) {
    idx = static_cast<int>(suppl->length()) + idx;
    if (idx < 0) {
      throw nb::index_error("invalid index");
    }
  }
  try {
    res = (*suppl)[static_cast<unsigned int>(idx)];
  } catch (...) {
    if (suppl->atEnd()) {
      throw nb::index_error("invalid index");
    } else {
      res = nullptr;
    }
  }
  return res;
}
}  // namespace

std::vector<std::vector<int>> GetResonanceSubstructMatches(
    ResonanceMolSupplier &suppl, const ROMol &query, bool uniquify = false,
    bool useChirality = false, bool useQueryQueryMatches = false,
    unsigned int maxMatches = 1000, int numThreads = 1) {
  std::vector<MatchVectType> matches;
  int matched =
      SubstructMatch(suppl, query, matches, uniquify, true, useChirality,
                     useQueryQueryMatches, maxMatches, numThreads);
  std::vector<std::vector<int>> res;
  res.reserve(matched);
  for (int idx = 0; idx < matched; idx++) {
    res.push_back(convertMatch(matches[idx]));
  }
  return res;
}

class PyResonanceMolSupplierCallback : public ResonanceMolSupplierCallback {
 public:
  PyResonanceMolSupplierCallback() : d_cppCallback(this) {}
  explicit PyResonanceMolSupplierCallback(nb::object pyCallbackObject)
      : d_cppCallback(this), d_pyCallbackObject(std::move(pyCallbackObject)) {
    auto *pyCallback =
        nb::cast<PyResonanceMolSupplierCallback *>(d_pyCallbackObject);
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
  inline nb::object getCallbackOverride() const {
    if (d_pyCallbackObject.is_valid()) {
      return d_pyCallbackObject.attr("__call__");
    }
    return nb::none();
  }

  bool operator()() override {
    PyGILStateHolder h;
    auto override = getCallbackOverride();
    if (!override.is_valid() || override.is_none()) {
      throw nb::attribute_error(
          "The __call__ attribute in the rdchem.ResonanceMolSupplierCallback "
          "subclass must exist and be a callable method");
    }
    return nb::cast<bool>(override());
  }
  nb::object getPyCallbackObject() { return d_pyCallbackObject; }

 private:
  PyResonanceMolSupplierCallback *d_cppCallback;
  nb::object d_pyCallbackObject;
};

nb::object getProgressCallbackHelper(const ResonanceMolSupplier &suppl) {
  PyResonanceMolSupplierCallback *cppCallback =
      dynamic_cast<PyResonanceMolSupplierCallback *>(
          suppl.getProgressCallback());
  nb::object res = nb::none();
  if (cppCallback) {
    res = cppCallback->getPyCallbackObject();
  }
  return res;
};

void setProgressCallbackHelper(ResonanceMolSupplier &suppl,
                               const nb::object &callback) {
  if (!callback.is_valid() || callback.is_none()) {
    suppl.setProgressCallback(nullptr);
    return;
  }
  try {
    auto *cb = nb::cast<PyResonanceMolSupplierCallback *>(callback);
    auto override = cb->getCallbackOverride();
    if (!override.is_valid() || !PyCallable_Check(override.ptr())) {
      throw nb::attribute_error(
          "The __call__ attribute in the rdchem.ResonanceMolSupplierCallback "
          "subclass must exist and be a callable method");
    }
    suppl.setProgressCallback(new PyResonanceMolSupplierCallback(callback));
  } catch (const nb::cast_error &) {
    throw nb::type_error(
        "Expected an instance of a rdchem.ResonanceMolSupplierCallback "
        "subclass");
  }
}

std::string resonanceMolSupplierCallbackClassDoc =
    R"DOC(Create a derived class from this abstract base class and
          implement the __call__() method.
          The __call__() method is called at each iteration of the
          algorithm, and provides a mechanism to monitor or stop
          its progress.

          To have your callback called, pass an instance of your
          derived class to ResonanceMolSupplier.SetProgressCallback()
)DOC";

std::string resonanceMolSupplierClassDoc =
    R"DOC(A class which supplies resonance structures (as mols) from a mol.

     Usage examples:

          1) Lazy evaluation: the resonance structures are not constructed
                until we ask for them:

                >>> suppl = ResonanceMolSupplier(mol)
                >>> for resMol in suppl:
                ...    resMol.GetNumAtoms()

          2) Lazy evaluation 2:

                >>> suppl = ResonanceMolSupplier(mol)
                >>> resMol1 = next(suppl)
                >>> resMol2 = next(suppl)
                >>> suppl.reset()
                >>> resMol3 = next(suppl)
                # resMol3 and resMol1 are the same:
                >>> MolToSmiles(resMol3)==MolToSmiles(resMol1)

          3) Random Access:

                >>> suppl = ResonanceMolSupplier(mol)
                >>> resMol1 = suppl[0]
                >>> resMol2 = suppl[1]

                NOTE: this will generate an IndexError if the supplier doesn't have that many
                molecules.

          4) Random Access 2: looping over all resonance structures
                >>> suppl = ResonanceMolSupplier(mol)
                >>> nResMols = len(suppl)
                >>> for i in range(nResMols):
                ...   suppl[i].GetNumAtoms()

)DOC";
struct resmolsup_wrap {
  static void wrap(nb::module_ &m) {
    nb::enum_<ResonanceMolSupplier::ResonanceFlags>(m, "ResonanceFlags")
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

    nb::class_<PyResonanceMolSupplierCallback>(
        m, "ResonanceMolSupplierCallback",
        resonanceMolSupplierCallbackClassDoc.c_str())
        .def(nb::init<>())
        .def(
            "GetNumConjGrps",
            &PyResonanceMolSupplierCallback::wrapGetNumConjGrps,
            R"DOC(Returns the number of individual conjugated groups in the molecule.
)DOC")
        .def("GetMaxStructures",
             &PyResonanceMolSupplierCallback::wrapGetMaxStructures,
             R"DOC(Get the number of conjugated groups this molecule has.
)DOC")
        .def(
            "GetNumStructures",
            (size_t (PyResonanceMolSupplierCallback::*)(
                unsigned int))&PyResonanceMolSupplierCallback::
                wrapGetNumStructures,
            "conjGrpIdx"_a,
            R"DOC(Get the number of resonance structures generated so far for the passed conjugated group index.
)DOC")
        .def(
            "GetNumDiverseStructures",
            (size_t (PyResonanceMolSupplierCallback::*)(
                unsigned int))&PyResonanceMolSupplierCallback::
                wrapGetNumDiverseStructures,
            "conjGrpIdx"_a,
            R"DOC(Get the number of non-degenrate resonance structures generated so far for the passed conjugated group index.
)DOC")
        .def(
            "__call__",
            [](PyResonanceMolSupplierCallback &) {
              throw nb::attribute_error(
                  "This must be implemented in the derived class");
              return false;
            },
            R"DOC(This must be implemented in the derived class. Return True if the resonance structure generation should continue; False if the resonance structure generation should stop.
)DOC");

    nb::class_<ResonanceMolSupplier>(m, "ResonanceMolSupplier",
                                     resonanceMolSupplierClassDoc.c_str())
        .def(nb::init<ROMol &, unsigned int, unsigned int>(), "mol"_a,
             "flags"_a = 0, "maxStructs"_a = 1000)
        .def(
            "__iter__",
            (ResonanceMolSupplier * (*)(ResonanceMolSupplier *)) & MolSupplIter,
            nb::rv_policy::reference_internal)
        .def(
            "__next__", (ROMol * (*)(ResonanceMolSupplier *)) & MolSupplNext,
            nb::rv_policy::take_ownership,
            R"DOC(Returns the next resonance structure in the supplier. Raises _StopIteration_ on end.
)DOC")
        .def("__getitem__",
             (ROMol * (*)(ResonanceMolSupplier *, int)) & MolSupplGetItem,
             nb::rv_policy::take_ownership, "idx"_a)
        .def(
            "reset", &ResonanceMolSupplier::reset,
            R"DOC(Resets our position in the resonance structure supplier to the beginning.
)DOC")
        .def("__len__", &ResonanceMolSupplier::length)
        .def(
            "atEnd", &ResonanceMolSupplier::atEnd,
            R"DOC(Returns whether or not we have hit the end of the resonance structure supplier.
)DOC")
        .def(
            "GetNumConjGrps", &ResonanceMolSupplier::getNumConjGrps,
            R"DOC(Returns the number of individual conjugated groups in the molecule.
)DOC")
        .def(
            "GetBondConjGrpIdx",
            (int (ResonanceMolSupplier::*)(
                unsigned int))&ResonanceMolSupplier::getBondConjGrpIdx,
            "bi"_a,
            R"DOC(Given a bond index, it returns the index of the conjugated groupthe bond belongs to, or -1 if it is not conjugated.
)DOC")
        .def(
            "GetAtomConjGrpIdx",
            (int (ResonanceMolSupplier::*)(
                unsigned int))&ResonanceMolSupplier::getAtomConjGrpIdx,
            "ai"_a,
            R"DOC(Given an atom index, it returns the index of the conjugated groupthe atom belongs to, or -1 if it is not conjugated.
)DOC")
        .def("SetNumThreads",
             (void (ResonanceMolSupplier::*)(
                 int))&ResonanceMolSupplier::setNumThreads,
             "numThreads"_a,
             R"DOC(Sets the number of threads to be used to enumerate resonance
structures (defaults to 1; 0 selects the number of concurrent
threads supported by the hardware; negative values are added
to the number of concurrent threads supported by the hardware).
)DOC")
        .def("SetProgressCallback", &setProgressCallbackHelper, "callback"_a,
             R"DOC(Pass an instance of a class derived from
ResonanceMolSupplierCallback, which must implement the
__call__() method.
)DOC")
        .def("GetProgressCallback", &getProgressCallbackHelper,
             R"DOC(Get the ResonanceMolSupplierCallback subclass instance,
or None if none was set.
)DOC")
        .def(
            "WasCanceled", &ResonanceMolSupplier::wasCanceled,
            R"DOC(Returns True if the resonance structure generation was canceled.
)DOC")
        .def(
            "Enumerate", &ResonanceMolSupplier::enumerate,
            R"DOC(Ask ResonanceMolSupplier to enumerate resonance structures(automatically done as soon as any attempt to access them is made).
)DOC")
        .def(
            "GetIsEnumerated", &ResonanceMolSupplier::getIsEnumerated,
            R"DOC(Returns true if resonance structure enumeration has already happened.
)DOC")
        .def(
            "GetSubstructMatch",
            (std::vector<int> (*)(ResonanceMolSupplier &, const ROMol &, bool,
                                  bool))GetSubstructMatch,
            "query"_a, "useChirality"_a = false,
            "useQueryQueryMatches"_a = false,
            R"DOC(Returns the indices of the molecule's atoms that match a substructure query,
taking into account all resonance structures in ResonanceMolSupplier.

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
            "GetSubstructMatches", GetResonanceSubstructMatches, "query"_a,
            "uniquify"_a = false, "useChirality"_a = false,
            "useQueryQueryMatches"_a = false, "maxMatches"_a = 1000,
            "numThreads"_a = 1,
            R"DOC(Returns tuples of the indices of the molecule's atoms that match a substructure query,
taking into account all resonance structures in ResonanceMolSupplier.

  ARGUMENTS:
    - query: a Molecule.
    - uniquify: (optional) determines whether or not the matches are uniquified.
                Defaults to 1.

    - useChirality: enables the use of stereochemistry in the matching

    - useQueryQueryMatches: use query-query matching logic

    - maxMatches: The maximum number of matches that will be returned.
                  In high-symmetry cases with medium-sized molecules, it is
                  very easy to end up with a combinatorial explosion in the
                  number of possible matches. This argument prevents that from
                  having unintended consequences

    - numThreads: The number of threads to be used (defaults to 1; 0 selects the
                  number of concurrent threads supported by the hardware; negative
                  values are added to the number of concurrent threads supported
                  by the hardware).

  RETURNS: a tuple of tuples of integers

  NOTE:
     - the ordering of the indices corresponds to the atom ordering
         in the query. For example, the first index is for the atom in
         this molecule that matches the first atom in the query.
)DOC");
  };
};
}  // namespace RDKit

void wrap_resmolsupplier(nb::module_ &m) { RDKit::resmolsup_wrap::wrap(m); }
