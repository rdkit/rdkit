//
//  Copyright (C) 2020-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/trampoline.h>
#include <nanobind/make_iterator.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolStandardize/MolStandardize.h>
#include <GraphMol/MolStandardize/Tautomer.h>
#include <sstream>

NB_MAKE_OPAQUE(std::vector<RDKit::MolStandardize::TautomerScoringFunctions::SubstructTerm>);

#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

namespace {

std::shared_ptr<RDKit::ROMol> toStd(const RDKit::ROMOL_SPTR &bptr) {
  return {bptr.get(), [b = bptr](RDKit::ROMol *) {}};
}

class PyTautomerEnumeratorResult {
 public:
  PyTautomerEnumeratorResult(const MolStandardize::TautomerEnumeratorResult &tr)
      : d_tr(new MolStandardize::TautomerEnumeratorResult(std::move(tr))) {
    nb::list atList;
    nb::list bndList;
    for (unsigned int i = 0; i < d_tr->modifiedAtoms().size(); ++i) {
      if (d_tr->modifiedAtoms().test(i)) {
        atList.append(i);
      }
    }
    for (unsigned int i = 0; i < d_tr->modifiedBonds().size(); ++i) {
      if (d_tr->modifiedBonds().test(i)) {
        bndList.append(i);
      }
    }
    d_atTuple = nb::tuple(atList);
    d_bndTuple = nb::tuple(bndList);
  }
  inline std::vector<std::shared_ptr<RDKit::ROMol>> tautomers() const {
    std::vector<std::shared_ptr<RDKit::ROMol>> out;
    for (const auto &sptr : d_tr->tautomers()) {
      out.push_back(toStd(sptr));
    }
    return out;
  }
  inline std::vector<std::string> smiles() const {
    return std::vector<std::string>(d_tr->smiles());
  }
  inline const MolStandardize::SmilesTautomerMap &smilesTautomerMap() const {
    return d_tr->smilesTautomerMap();
  }
  inline MolStandardize::TautomerEnumeratorStatus status() const {
    return d_tr->status();
  }
  inline nb::tuple modifiedAtoms() const { return d_atTuple; }
  inline nb::tuple modifiedBonds() const { return d_bndTuple; }
  inline const MolStandardize::TautomerEnumeratorResult::const_iterator begin()
      const {
    return d_tr->begin();
  }
  inline const MolStandardize::TautomerEnumeratorResult::const_iterator end()
      const {
    return d_tr->end();
  }
  inline int size() const { return d_tr->size(); }
  RDKit::ROMol *at(int pos) const {
    if (pos < 0) {
      pos += size();
    }
    if (pos < 0 || pos >= size()) {
      PyErr_SetString(PyExc_IndexError, "index out of bounds");
      throw nb::python_error();
    }
    return new RDKit::ROMol(*d_tr->at(pos));
  }
  const MolStandardize::TautomerEnumeratorResult *get() { return d_tr.get(); }

 private:
  std::shared_ptr<MolStandardize::TautomerEnumeratorResult> d_tr;
  nb::tuple d_atTuple;
  nb::tuple d_bndTuple;
};

struct TautomerEnumeratorCallbackTrampoline
    : MolStandardize::TautomerEnumeratorCallback {
  NB_TRAMPOLINE(MolStandardize::TautomerEnumeratorCallback, 1);

  bool operator()(const ROMol &mol,
                  const MolStandardize::TautomerEnumeratorResult &res) override {
    nb::gil_scoped_acquire guard;
    nb::object self = nb::find(this);
    auto *pyRes = new PyTautomerEnumeratorResult(res);
    nb::object pyResObj = nb::cast(pyRes, nb::rv_policy::take_ownership);
    return nb::cast<bool>(self.attr("__call__")(mol, pyResObj));
  }
};

// Wrapper to hold a Python callback object without needing to copy trampolines
struct PyCallbackWrapper : MolStandardize::TautomerEnumeratorCallback {
  nb::object d_pyObj;
  explicit PyCallbackWrapper(nb::object obj) : d_pyObj(std::move(obj)) {}
  bool operator()(const ROMol &mol,
                  const MolStandardize::TautomerEnumeratorResult &res) override {
    nb::gil_scoped_acquire guard;
    auto *pyRes = new PyTautomerEnumeratorResult(res);
    nb::object pyResObj = nb::cast(pyRes, nb::rv_policy::take_ownership);
    return nb::cast<bool>(d_pyObj(mol, pyResObj));
  }
};

nb::tuple smilesTautomerMapKeysHelper(
    const MolStandardize::SmilesTautomerMap &self) {
  nb::list keys;
  for (const auto &pair : self) {
    keys.append(pair.first);
  }
  return nb::tuple(keys);
}

nb::tuple smilesTautomerMapValuesHelper(
    const MolStandardize::SmilesTautomerMap &self) {
  nb::list values;
  for (const auto &pair : self) {
    auto *t = new MolStandardize::Tautomer(pair.second);
    values.append(nb::cast(t, nb::rv_policy::take_ownership));
  }
  return nb::tuple(values);
}

nb::tuple smilesTautomerMapItemsHelper(
    const MolStandardize::SmilesTautomerMap &self) {
  nb::list items;
  for (const auto &pair : self) {
    auto *t = new MolStandardize::Tautomer(pair.second);
    items.append(nb::make_tuple(
        pair.first, nb::cast(t, nb::rv_policy::take_ownership)));
  }
  return nb::tuple(items);
}

nb::object getCallbackHelper(const MolStandardize::TautomerEnumerator &te) {
  PyCallbackWrapper *cppCallback =
      dynamic_cast<PyCallbackWrapper *>(te.getCallback());
  if (cppCallback) {
    return cppCallback->d_pyObj;
  }
  return nb::none();
}

void setCallbackHelper(MolStandardize::TautomerEnumerator &te,
                       nb::object callback) {
  if (callback.is_none()) {
    te.setCallback(nullptr);
    return;
  }
  if (!nb::isinstance(
          callback,
          nb::type<MolStandardize::TautomerEnumeratorCallback>())) {
    throw nb::type_error(
        "Expected an instance of a "
        "rdMolStandardize.TautomerEnumeratorCallback subclass");
  }
  // Verify that the Python subclass has a properly overridden __call__:
  // - it must be defined in the immediate class dict (not just inherited)
  // - it must be callable
  nb::object cls_dict =
      nb::borrow<nb::object>(PyObject_Type(callback.ptr())).attr("__dict__");
  if (!nb::cast<bool>(cls_dict.attr("__contains__")("__call__"))) {
    throw nb::attribute_error(
        "TautomerEnumeratorCallback subclass must override __call__");
  }
  nb::object call_attr = cls_dict.attr("__getitem__")("__call__");
  nb::object builtins = nb::module_::import_("builtins");
  if (!nb::cast<bool>(builtins.attr("callable")(call_attr))) {
    throw nb::attribute_error(
        "TautomerEnumeratorCallback.__call__ must be callable");
  }
  te.setCallback(new PyCallbackWrapper(callback));
}

class pyobjFunctor {
 public:
  pyobjFunctor(nb::object obj) : dp_obj(std::move(obj)) {}
  ~pyobjFunctor() = default;
  int operator()(const ROMol &m) {
    nb::object result = dp_obj(m);
    try {
      return nb::cast<int>(result);
    } catch (const nb::cast_error &) {
      throw nb::type_error(
          "scoring function must return a numeric value");
    }
  }

 private:
  nb::object dp_obj;
};

ROMol *canonicalizeHelper(const MolStandardize::TautomerEnumerator &self,
                          const ROMol &mol) {
  return self.canonicalize(mol);
}

ROMol *canonicalizeHelper2(const MolStandardize::TautomerEnumerator &self,
                           const ROMol &mol, nb::object scoreFunc) {
  pyobjFunctor ftor(scoreFunc);
  return self.canonicalize(mol, ftor);
}

inline std::vector<ROMOL_SPTR> extractPythonIterable(const nb::object &o) {
  if (!nb::hasattr(o, "__iter__")) {
    throw nb::type_error(
        "the passed object should be an iterable of Mol objects");
  }
  std::vector<ROMOL_SPTR> result;
  for (nb::handle h : o) {
    RDKit::ROMol *molPtr = nullptr;
    try {
      molPtr = nb::cast<RDKit::ROMol *>(h);
    } catch (const nb::cast_error &) {
      throw nb::type_error(
          "the passed object should be an iterable of Mol objects");
    }
    // Copy the molecule so the vector owns its lifetime independently of Python
    result.push_back(ROMOL_SPTR(new RDKit::ROMol(*molPtr)));
  }
  return result;
}

ROMol *pickCanonicalHelper(const MolStandardize::TautomerEnumerator &self,
                           const nb::object &o) {
  try {
    PyTautomerEnumeratorResult *e =
        nb::cast<PyTautomerEnumeratorResult *>(o);
    return self.pickCanonical(*e->get());
  } catch (const nb::cast_error &) {
    return self.pickCanonical(extractPythonIterable(o));
  }
}

ROMol *pickCanonicalHelper2(const MolStandardize::TautomerEnumerator &self,
                            const nb::object &o, nb::object scoreFunc) {
  pyobjFunctor ftor(scoreFunc);
  try {
    PyTautomerEnumeratorResult *e =
        nb::cast<PyTautomerEnumeratorResult *>(o);
    return self.pickCanonical(*e->get(), ftor);
  } catch (const nb::cast_error &) {
    return self.pickCanonical(extractPythonIterable(o), ftor);
  }
}

PyTautomerEnumeratorResult *enumerateHelper(
    const MolStandardize::TautomerEnumerator &self, const ROMol &mol) {
  return new PyTautomerEnumeratorResult(self.enumerate(mol));
}

std::vector<MolStandardize::TautomerScoringFunctions::SubstructTerm>
GetDefaultTautomerSubstructsHelper() {
  std::vector<MolStandardize::TautomerScoringFunctions::SubstructTerm> terms;
  for (auto term : MolStandardize::TautomerScoringFunctions::
           getDefaultTautomerScoreSubstructs()) {
    terms.emplace_back(term);
  }
  return terms;
}

}  // namespace

void wrap_tautomer(nb::module_ &m) {
  nb::enum_<MolStandardize::TautomerEnumeratorStatus>(
      m, "TautomerEnumeratorStatus")
      .value("Completed",
             MolStandardize::TautomerEnumeratorStatus::Completed)
      .value("MaxTautomersReached",
             MolStandardize::TautomerEnumeratorStatus::MaxTautomersReached)
      .value("MaxTransformsReached",
             MolStandardize::TautomerEnumeratorStatus::MaxTransformsReached)
      .value("Canceled",
             MolStandardize::TautomerEnumeratorStatus::Canceled);

  nb::class_<MolStandardize::TautomerEnumeratorCallback,
             TautomerEnumeratorCallbackTrampoline>(
      m, "TautomerEnumeratorCallback",
      R"DOC(Create a derived class from this abstract base class and
implement the __call__() method.
The __call__() method is called in the innermost loop of the
algorithm, and provides a mechanism to monitor or stop
its progress.

To have your callback called, pass an instance of your
derived class to TautomerEnumerator.SetCallback())DOC")
      .def(nb::init<>())
      .def("__call__",
           &MolStandardize::TautomerEnumeratorCallback::operator(), "mol"_a,
           "res"_a,
           "This must be implemented in the derived class. "
           "Return True if the tautomer enumeration should continue; "
           "False if the tautomer enumeration should stop.\n");

  nb::class_<MolStandardize::Tautomer>(
      m, "Tautomer",
      "used to hold the aromatic and kekulized versions of each tautomer")
      .def_prop_ro("tautomer",
                   [](const MolStandardize::Tautomer &self) {
                     return new ROMol(*self.tautomer);
                   },
                   nb::rv_policy::take_ownership,
                   "aromatic version of the tautomer")
      .def_prop_ro("kekulized",
                   [](const MolStandardize::Tautomer &self) {
                     return new ROMol(*self.getKekulized());
                   },
                   nb::rv_policy::take_ownership,
                   "kekulized version of the tautomer");

  nb::class_<MolStandardize::SmilesTautomerMap>(
      m, "SmilesTautomerMap",
      "maps SMILES strings to the respective Tautomer objects")
      .def("keys", &smilesTautomerMapKeysHelper)
      .def("values", &smilesTautomerMapValuesHelper)
      .def("items", &smilesTautomerMapItemsHelper)
      .def("__len__",
           [](const MolStandardize::SmilesTautomerMap &self) {
             return self.size();
           });

  nb::class_<PyTautomerEnumeratorResult>(
      m, "TautomerEnumeratorResult",
      "used to return tautomer enumeration results")
      .def_prop_ro("tautomers", &PyTautomerEnumeratorResult::tautomers,
                   "tautomers generated by the enumerator")
      .def_prop_ro("smiles", &PyTautomerEnumeratorResult::smiles,
                   "SMILES of tautomers generated by the enumerator")
      .def_prop_ro(
          "smilesTautomerMap",
          &PyTautomerEnumeratorResult::smilesTautomerMap,
          nb::rv_policy::reference_internal,
          "dictionary mapping SMILES strings to the respective Tautomer objects")
      .def_prop_ro("status", &PyTautomerEnumeratorResult::status,
                   "whether the enumeration completed or not; see "
                   "TautomerEnumeratorStatus for possible values")
      .def_prop_ro("modifiedAtoms", &PyTautomerEnumeratorResult::modifiedAtoms,
                   "tuple of atom indices modified by the transforms")
      .def_prop_ro("modifiedBonds", &PyTautomerEnumeratorResult::modifiedBonds,
                   "tuple of bond indices modified by the transforms")
      .def("__call__",
           [](PyTautomerEnumeratorResult &self) { return self.tautomers(); },
           "tautomers generated by the enumerator")
      .def("__iter__",
           [](PyTautomerEnumeratorResult &self) {
             auto taus = self.tautomers();
             nb::list result;
             for (const auto &mol : taus) {
               result.append(nb::cast(mol));
             }
             return result.attr("__iter__")();
           })
      .def("__getitem__",
           [](PyTautomerEnumeratorResult &self, int pos) {
             return self.at(pos);
           },
           "pos"_a, nb::rv_policy::take_ownership)
      .def("__len__", &PyTautomerEnumeratorResult::size);

  nb::class_<MolStandardize::TautomerEnumerator>(m, "TautomerEnumerator")
      .def(nb::init<>())
      .def(nb::init<const MolStandardize::CleanupParameters &>(), "params"_a)
      .def(nb::init<const MolStandardize::TautomerEnumerator &>(), "other"_a)
      .def("Enumerate", &enumerateHelper, "mol"_a,
           nb::rv_policy::take_ownership,
           R"DOC(Generates the tautomers for a molecule.

The enumeration rules are inspired by the publication:
M. Sitzmann et al., "Tautomerism in Large Databases.", JCAMD 24:521 (2010)
https://doi.org/10.1007/s10822-010-9346-4

Note: the definitions used here are that the atoms modified during
tautomerization are the atoms at the beginning and end of each tautomer
transform (the H "donor" and H "acceptor" in the transform) and the bonds
modified during transformation are any bonds whose order is changed during
the tautomer transform (these are the bonds between the "donor" and the
"acceptor").)DOC")
      .def("Canonicalize", &canonicalizeHelper, "mol"_a,
           nb::rv_policy::take_ownership,
           R"DOC(Returns the canonical tautomer for a molecule.

The default scoring scheme is inspired by the publication:
M. Sitzmann et al., "Tautomerism in Large Databases.", JCAMD 24:521 (2010)
https://doi.org/10.1007/s10822-010-9346-4

Note that the canonical tautomer is very likely not the most stable tautomer
for any given conditions. The default scoring rules are designed to produce
"reasonable" tautomers, but the primary concern is that the results are
canonical: you always get the same canonical tautomer for a molecule
regardless of what the input tautomer or atom ordering were.)DOC")
      .def("Canonicalize", &canonicalizeHelper2, "mol"_a, "scoreFunc"_a,
           nb::rv_policy::take_ownership,
           "picks the canonical tautomer from an iterable of molecules "
           "using a custom scoring function")
      .def("PickCanonical", &pickCanonicalHelper, "iterable"_a,
           nb::rv_policy::take_ownership,
           "picks the canonical tautomer from an iterable of molecules")
      .def("PickCanonical", &pickCanonicalHelper2, "iterable"_a,
           "scoreFunc"_a, nb::rv_policy::take_ownership,
           "returns the canonical tautomer for a molecule using a custom "
           "scoring function")
      .def_static(
          "ScoreTautomer",
          &MolStandardize::TautomerScoringFunctions::scoreTautomer, "mol"_a,
          "returns the score for a tautomer using the default scoring scheme.")
      .def("SetMaxTautomers",
           &MolStandardize::TautomerEnumerator::setMaxTautomers,
           "maxTautomers"_a,
           "set the maximum number of tautomers to be generated.")
      .def("GetMaxTautomers",
           &MolStandardize::TautomerEnumerator::getMaxTautomers,
           "returns the maximum number of tautomers to be generated.")
      .def("SetMaxTransforms",
           &MolStandardize::TautomerEnumerator::setMaxTransforms,
           "maxTransforms"_a,
           "set the maximum number of transformations to be applied. "
           "This limit is usually hit earlier than the maxTautomers limit "
           "and leads to a more linear scaling of CPU time with increasing "
           "number of tautomeric centers (see Sitzmann et al.).")
      .def("GetMaxTransforms",
           &MolStandardize::TautomerEnumerator::getMaxTransforms,
           "returns the maximum number of transformations to be applied.")
      .def("SetRemoveSp3Stereo",
           &MolStandardize::TautomerEnumerator::setRemoveSp3Stereo,
           "removeSp3Stereo"_a,
           "set to True if you wish stereochemistry information "
           "to be removed from sp3 atoms involved in tautomerism. "
           "This means that S-aminoacids will lose their stereochemistry "
           "after going through tautomer enumeration because of the "
           "amido-imidol tautomerism. This defaults to True in RDKit, "
           "and to False in the workflow described by Sitzmann et al.")
      .def("GetRemoveSp3Stereo",
           &MolStandardize::TautomerEnumerator::getRemoveSp3Stereo,
           "returns whether stereochemistry information will be removed from "
           "sp3 atoms involved in tautomerism.")
      .def("SetRemoveBondStereo",
           &MolStandardize::TautomerEnumerator::setRemoveBondStereo,
           "removeBondStereo"_a,
           "set to True if you wish stereochemistry information "
           "to be removed from double bonds involved in tautomerism. "
           "This means that enols will lose their E/Z stereochemistry "
           "after going through tautomer enumeration because of the "
           "keto-enolic tautomerism. This defaults to True in the "
           "RDKit and also in the workflow described by Sitzmann et al.")
      .def("GetRemoveBondStereo",
           &MolStandardize::TautomerEnumerator::getRemoveBondStereo,
           "returns whether stereochemistry information "
           "will be removed from double bonds involved in tautomerism.")
      .def("SetReassignStereo",
           &MolStandardize::TautomerEnumerator::setReassignStereo,
           "reassignStereo"_a,
           "set to True if you wish AssignStereochemistry to be called "
           "on each tautomer generated by the Enumerate() method. "
           "This defaults to True.")
      .def("GetReassignStereo",
           &MolStandardize::TautomerEnumerator::getReassignStereo,
           "returns whether AssignStereochemistry will be called "
           "on each tautomer generated by the Enumerate() method.")
      .def("SetCallback", &setCallbackHelper, "callback"_a,
           "Pass an instance of a class derived from\n"
           "TautomerEnumeratorCallback, which must implement the\n"
           "__call__() method.")
      .def("GetCallback", &getCallbackHelper,
           "Get the TautomerEnumeratorCallback subclass instance,\n"
           "or None if none was set.")
      .def_prop_ro_static(
          "tautomerScoreVersion",
          [](nb::handle) {
            return MolStandardize::TautomerScoringFunctions::
                tautomerScoringVersion;
          });

  m.def("GetV1TautomerEnumerator", MolStandardize::getV1TautomerEnumerator,
        "return a TautomerEnumerator using v1 of the enumeration rules",
        nb::rv_policy::take_ownership);

  m.def("ScoreRings",
        MolStandardize::TautomerScoringFunctions::scoreRings, "mol"_a,
        "scores the ring system of the tautomer for canonicalization\n"
        "Aromatic rings score 100, all carbon aromatic rings score 250");

  m.def("ScoreHeteroHs",
        MolStandardize::TautomerScoringFunctions::scoreHeteroHs, "mol"_a,
        "scores the number of heteroHs of the tautomer for canonicalization\n"
        "This gives a negative penalty to hydrogens attached to S,P, Se and Te");

  nb::class_<MolStandardize::TautomerScoringFunctions::SubstructTerm>(
      m, "SubstructTerm",
      "Sets the score of this particular tautomer substructure, higher scores "
      "are more preferable\n"
      "Aromatic rings score 100, all carbon aromatic rings score 250")
      .def(nb::init<std::string, std::string, int>(), "name"_a, "smarts"_a,
           "score"_a)
      .def_ro(
          "name",
          &MolStandardize::TautomerScoringFunctions::SubstructTerm::name)
      .def_ro(
          "smarts",
          &MolStandardize::TautomerScoringFunctions::SubstructTerm::smarts)
      .def_ro(
          "score",
          &MolStandardize::TautomerScoringFunctions::SubstructTerm::score);

  nb::bind_vector<
      std::vector<MolStandardize::TautomerScoringFunctions::SubstructTerm>>(
      m, "SubstructTermVector");

  m.def(
      "ScoreSubstructs",
      [](const ROMol &mol,
         nb::object terms) -> int {
        if (terms.is_none()) {
          return MolStandardize::TautomerScoringFunctions::scoreSubstructs(mol);
        }
        return MolStandardize::TautomerScoringFunctions::scoreSubstructs(
            mol,
            nb::cast<std::vector<
                MolStandardize::TautomerScoringFunctions::SubstructTerm> &>(
                terms));
      },
      "mol"_a, "terms"_a = nb::none(), "scores the tautomer substructures");

  m.def("GetDefaultTautomerScoreSubstructs",
        GetDefaultTautomerSubstructsHelper,
        "Return the default tautomer substructure scoring terms");
}
