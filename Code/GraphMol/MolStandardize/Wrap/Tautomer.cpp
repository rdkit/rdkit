//
//  Copyright (C) 2020 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define NO_IMPORT_ARRAY
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <RDBoost/Wrap.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolStandardize/MolStandardize.h>
#include <GraphMol/MolStandardize/Tautomer.h>
#include <sstream>

namespace python = boost::python;
using namespace RDKit;

namespace {
class PyTautomerEnumeratorResult {
 public:
  PyTautomerEnumeratorResult(const MolStandardize::TautomerEnumeratorResult &tr)
      : d_tr(new MolStandardize::TautomerEnumeratorResult(std::move(tr))) {
    python::list atList;
    python::list bndList;
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
    d_atTuple = python::tuple(atList);
    d_bndTuple = python::tuple(bndList);
  }
  inline const std::vector<ROMOL_SPTR> *tautomers() const {
    return new std::vector<ROMOL_SPTR>(d_tr->tautomers());
  }
  inline const std::vector<std::string> *smiles() const {
    return new std::vector<std::string>(d_tr->smiles());
  }
  inline const MolStandardize::SmilesTautomerMap &smilesTautomerMap() const {
    return d_tr->smilesTautomerMap();
  }
  inline MolStandardize::TautomerEnumeratorStatus status() const {
    return d_tr->status();
  }
  inline python::tuple modifiedAtoms() const { return d_atTuple; }
  inline python::tuple modifiedBonds() const { return d_bndTuple; }
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
      python::throw_error_already_set();
      return nullptr;
    }
    return new RDKit::ROMol(*d_tr->at(pos));
  }
  const MolStandardize::TautomerEnumeratorResult *get() { return d_tr.get(); }

 private:
  boost::shared_ptr<MolStandardize::TautomerEnumeratorResult> d_tr;
  python::tuple d_atTuple;
  python::tuple d_bndTuple;
};

class PyTautomerEnumeratorCallback
    : public MolStandardize::TautomerEnumeratorCallback,
      public python::wrapper<MolStandardize::TautomerEnumeratorCallback> {
 public:
  PyTautomerEnumeratorCallback() {}
  PyTautomerEnumeratorCallback(const python::object &pyCallbackObject) {
    PyTautomerEnumeratorCallback *pyCallback =
        python::extract<PyTautomerEnumeratorCallback *>(pyCallbackObject);
    *this = *pyCallback;
    d_pyCallbackObject = pyCallbackObject;
    pyCallback->d_cppCallback = this;
  }
  inline python::object getCallbackOverride() const {
    return get_override("__call__");
  }
  bool operator()(
      const ROMol &mol,
      const MolStandardize::TautomerEnumeratorResult &res) override {
    PyTautomerEnumeratorResult pyRes(res);
    return getCallbackOverride()(boost::ref(mol), boost::ref(pyRes));
  }
  python::object getPyCallbackObject() { return d_pyCallbackObject; }

 private:
  PyTautomerEnumeratorCallback *d_cppCallback;
  python::object d_pyCallbackObject;
};

typedef boost::shared_ptr<MolStandardize::Tautomer> TAUT_SPTR;

ROMol *getTautomerHelper(const TAUT_SPTR &self) {
  return new ROMol(*self->tautomer);
}

ROMol *getKekulizedHelper(const TAUT_SPTR &self) {
  return new ROMol(*self->kekulized);
}

python::tuple smilesTautomerMapKeysHelper(
    const MolStandardize::SmilesTautomerMap &self) {
  python::list keys;
  for (const auto &pair : self) {
    keys.append(pair.first);
  }
  return python::tuple(keys);
}

python::tuple smilesTautomerMapValuesHelper(
    const MolStandardize::SmilesTautomerMap &self) {
  python::list values;
  for (const auto &pair : self) {
    values.append(TAUT_SPTR(new MolStandardize::Tautomer(pair.second)));
  }
  return python::tuple(values);
}

python::tuple smilesTautomerMapItemsHelper(
    const MolStandardize::SmilesTautomerMap &self) {
  python::list items;
  for (const auto &pair : self) {
    items.append(python::make_tuple(
        pair.first, TAUT_SPTR(new MolStandardize::Tautomer(pair.second))));
  }
  return python::tuple(items);
}

python::object getCallbackHelper(const MolStandardize::TautomerEnumerator &te) {
  PyTautomerEnumeratorCallback *cppCallback =
      dynamic_cast<PyTautomerEnumeratorCallback *>(te.getCallback());
  python::object res;
  if (cppCallback) {
    res = cppCallback->getPyCallbackObject();
  }
  return res;
};

void setCallbackHelper(MolStandardize::TautomerEnumerator &te,
                       PyObject *callback) {
  PRECONDITION(callback, "callback must not be NULL");
  if (callback == Py_None) {
    te.setCallback(nullptr);
    return;
  }
  python::object callbackObject(python::handle<>(python::borrowed(callback)));
  python::extract<PyTautomerEnumeratorCallback *> extractCallback(
      callbackObject);
  if (extractCallback.check()) {
    if (!PyCallable_Check(extractCallback()->getCallbackOverride().ptr())) {
      PyErr_SetString(PyExc_AttributeError,
                      "The __call__ attribute in the "
                      "rdMolStandardize.TautomerEnumeratorCallback subclass "
                      "must exist and be a callable method");
      python::throw_error_already_set();
    } else {
      te.setCallback(new PyTautomerEnumeratorCallback(callbackObject));
    }
  } else {
    PyErr_SetString(PyExc_TypeError,
                    "Expected an instance of a "
                    "rdMolStandardize.TautomerEnumeratorCallback subclass");
    python::throw_error_already_set();
  }
}

std::string tautomerEnumeratorCallbackClassDoc =
    R"DOC(Create a derived class from this abstract base class and
    implement the __call__() method.
    The __call__() method is called in the innermost loop of the
    algorithm, and provides a mechanism to monitor or stop
    its progress.

    To have your callback called, pass an instance of your
    derived class to TautomerEnumerator.SetCallback())DOC";

MolStandardize::TautomerEnumerator *EnumeratorFromParams(
    const MolStandardize::CleanupParameters &params) {
  return new MolStandardize::TautomerEnumerator(params);
}

MolStandardize::TautomerEnumerator *createDefaultEnumerator() {
  MolStandardize::CleanupParameters ps;
  return EnumeratorFromParams(ps);
}

MolStandardize::TautomerEnumerator *copyEnumerator(
    const MolStandardize::TautomerEnumerator &other) {
  return new MolStandardize::TautomerEnumerator(other);
}

class pyobjFunctor {
 public:
  pyobjFunctor(python::object obj) : dp_obj(std::move(obj)) {}
  ~pyobjFunctor() = default;
  int operator()(const ROMol &m) {
    return python::extract<int>(dp_obj(boost::ref(m)));
  }

 private:
  python::object dp_obj;
};

ROMol *canonicalizeHelper(const MolStandardize::TautomerEnumerator &self,
                          const ROMol &mol) {
  return self.canonicalize(mol);
}

ROMol *canonicalizeHelper2(const MolStandardize::TautomerEnumerator &self,
                           const ROMol &mol, python::object scoreFunc) {
  pyobjFunctor ftor(scoreFunc);
  return self.canonicalize(mol, ftor);
}

inline std::vector<ROMOL_SPTR> extractPythonIterable(const python::object &o) {
  if (!PyObject_HasAttrString(o.ptr(), "__iter__")) {
    PyErr_SetString(PyExc_TypeError,
                    "the passed object should be an iterable of Mol objects");
    python::throw_error_already_set();
    return std::vector<ROMOL_SPTR>();
  }
  return std::vector<ROMOL_SPTR>(python::stl_input_iterator<ROMOL_SPTR>(o),
                                 python::stl_input_iterator<ROMOL_SPTR>());
}

ROMol *pickCanonicalHelper(const MolStandardize::TautomerEnumerator &self,
                           const python::object &o) {
  python::extract<PyTautomerEnumeratorResult *> e(o);
  if (e.check()) {
    return self.pickCanonical(*e()->get());
  }
  return self.pickCanonical(extractPythonIterable(o));
}

ROMol *pickCanonicalHelper2(const MolStandardize::TautomerEnumerator &self,
                            const python::object &o, python::object scoreFunc) {
  pyobjFunctor ftor(scoreFunc);
  python::extract<PyTautomerEnumeratorResult *> e(o);
  if (e.check()) {
    return self.pickCanonical(*e()->get(), ftor);
  }
  return self.pickCanonical(extractPythonIterable(o), ftor);
}

PyTautomerEnumeratorResult *enumerateHelper(
    const MolStandardize::TautomerEnumerator &self, const ROMol &mol) {
  return new PyTautomerEnumeratorResult(self.enumerate(mol));
}

}  // namespace

struct tautomer_wrapper {
  static void wrap() {
    python::enum_<MolStandardize::TautomerEnumeratorStatus>(
        "TautomerEnumeratorStatus")
        .value("Completed", MolStandardize::TautomerEnumeratorStatus::Completed)
        .value("MaxTautomersReached",
               MolStandardize::TautomerEnumeratorStatus::MaxTautomersReached)
        .value("MaxTransformsReached",
               MolStandardize::TautomerEnumeratorStatus::MaxTransformsReached)
        .value("Canceled", MolStandardize::TautomerEnumeratorStatus::Canceled);

    python::class_<PyTautomerEnumeratorCallback, boost::noncopyable>(
        "TautomerEnumeratorCallback",
        tautomerEnumeratorCallbackClassDoc.c_str(), python::init<>())
        .def("__call__",
             python::pure_virtual(&PyTautomerEnumeratorCallback::operator()),
             "This must be implemented in the derived class. "
             "Return True if the tautomer enumeration should continue; "
             "False if the tautomer enumeration should stop.\n");

    python::class_<MolStandardize::Tautomer, TAUT_SPTR>(
        "Tautomer",
        "used to hold the aromatic and kekulized versions "
        "of each tautomer",
        python::no_init)
        .add_property(
            "tautomer",
            python::make_function(
                &getTautomerHelper,
                python::return_value_policy<python::manage_new_object>()),
            "aromatic version of the tautomer")
        .add_property(
            "kekulized",
            python::make_function(
                &getKekulizedHelper,
                python::return_value_policy<python::manage_new_object>()),
            "kekulized version of the tautomer");

    python::class_<MolStandardize::SmilesTautomerMap, boost::noncopyable>(
        "SmilesTautomerMap",
        "maps SMILES strings to the respective Tautomer objects",
        python::no_init)
        .def(python::map_indexing_suite<MolStandardize::SmilesTautomerMap,
                                        true>())
        .def("keys", &smilesTautomerMapKeysHelper)
        .def("values", &smilesTautomerMapValuesHelper)
        .def("items", &smilesTautomerMapItemsHelper);

    python::class_<PyTautomerEnumeratorResult, boost::noncopyable>(
        "TautomerEnumeratorResult",
        "used to return tautomer enumeration results", python::no_init)
        .add_property(
            "tautomers",
            python::make_function(
                &PyTautomerEnumeratorResult::tautomers,
                python::return_value_policy<python::manage_new_object>()),
            "tautomers generated by the enumerator")
        .add_property(
            "smiles",
            python::make_function(
                &PyTautomerEnumeratorResult::smiles,
                python::return_value_policy<python::manage_new_object>()),
            "SMILES of tautomers generated by the enumerator")
        .add_property("smilesTautomerMap",
                      python::make_function(
                          &PyTautomerEnumeratorResult::smilesTautomerMap,
                          python::return_value_policy<
                              python::reference_existing_object>()),
                      "dictionary mapping SMILES strings to the respective "
                      "Tautomer objects")
        .add_property("status", &PyTautomerEnumeratorResult::status,
                      "whether the enumeration completed or not; see "
                      "TautomerEnumeratorStatus for possible values")
        .add_property(
            "modifiedAtoms",
            python::make_function(&PyTautomerEnumeratorResult::modifiedAtoms),
            "tuple of atom indices modified by the transforms")
        .add_property(
            "modifiedBonds",
            python::make_function(&PyTautomerEnumeratorResult::modifiedBonds),
            "tuple of bond indices modified by the transforms")
        .def("__call__", &PyTautomerEnumeratorResult::tautomers,
             python::return_value_policy<python::manage_new_object>(),
             "tautomers generated by the enumerator")
        .def("__iter__", python::range(&PyTautomerEnumeratorResult::begin,
                                       &PyTautomerEnumeratorResult::end))
        .def("__getitem__", &PyTautomerEnumeratorResult::at,
             python::return_value_policy<python::manage_new_object>())
        .def("__len__", &PyTautomerEnumeratorResult::size);

    python::class_<MolStandardize::TautomerEnumerator, boost::noncopyable>(
        "TautomerEnumerator", python::no_init)
        .def("__init__", python::make_constructor(createDefaultEnumerator))
        .def("__init__", python::make_constructor(EnumeratorFromParams))
        .def("__init__", python::make_constructor(copyEnumerator))
        .def("Enumerate", &enumerateHelper,
             (python::arg("self"), python::arg("mol")),
             python::return_value_policy<python::manage_new_object>(),
             R"DOC(Generates the tautomers for a molecule.
             
  The enumeration rules are inspired by the publication:
  M. Sitzmann et al., “Tautomerism in Large Databases.”, JCAMD 24:521 (2010)
  https://doi.org/10.1007/s10822-010-9346-4
  
  Note: the definitions used here are that the atoms modified during
  tautomerization are the atoms at the beginning and end of each tautomer
  transform (the H "donor" and H "acceptor" in the transform) and the bonds
  modified during transformation are any bonds whose order is changed during
  the tautomer transform (these are the bonds between the "donor" and the
  "acceptor").)DOC")
        .def("Canonicalize", &canonicalizeHelper,
             (python::arg("self"), python::arg("mol")),
             python::return_value_policy<python::manage_new_object>(),
             R"DOC(Returns the canonical tautomer for a molecule.

  The default scoring scheme is inspired by the publication:
  M. Sitzmann et al., “Tautomerism in Large Databases.”, JCAMD 24:521 (2010)
  https://doi.org/10.1007/s10822-010-9346-4

  Note that the canonical tautomer is very likely not the most stable tautomer
  for any given conditions. The default scoring rules are designed to produce
  "reasonable" tautomers, but the primary concern is that the results are
  canonical: you always get the same canonical tautomer for a molecule
  regardless of what the input tautomer or atom ordering were.)DOC")
        .def(
            "Canonicalize", &canonicalizeHelper2,
            (python::arg("self"), python::arg("mol"), python::arg("scoreFunc")),
            python::return_value_policy<python::manage_new_object>(),
            "picks the canonical tautomer from an iterable of molecules "
            "using a custom scoring function")
        .def("PickCanonical", &pickCanonicalHelper,
             (python::arg("self"), python::arg("iterable")),
             python::return_value_policy<python::manage_new_object>(),
             "picks the canonical tautomer from an iterable of molecules")
        .def("PickCanonical", &pickCanonicalHelper2,
             (python::arg("self"), python::arg("iterable"),
              python::arg("scoreFunc")),
             python::return_value_policy<python::manage_new_object>(),
             "returns the canonical tautomer for a molecule using a custom "
             "scoring function")
        .def("ScoreTautomer",
             &MolStandardize::TautomerScoringFunctions::scoreTautomer,
             (python::arg("mol")),
             "returns the score for a tautomer using the default scoring "
             "scheme.")
        .staticmethod("ScoreTautomer")
        .def("SetMaxTautomers",
             &MolStandardize::TautomerEnumerator::setMaxTautomers,
             (python::arg("self"), python::arg("maxTautomers")),
             "set the maximum number of tautomers to be generated.")
        .def("GetMaxTautomers",
             &MolStandardize::TautomerEnumerator::getMaxTautomers,
             (python::arg("self")),
             "returns the maximum number of tautomers to be generated.")
        .def("SetMaxTransforms",
             &MolStandardize::TautomerEnumerator::setMaxTransforms,
             (python::arg("self"), python::arg("maxTransforms")),
             "set the maximum number of transformations to be applied. "
             "This limit is usually hit earlier than the maxTautomers limit "
             "and leads to a more linear scaling of CPU time with increasing "
             "number of tautomeric centers (see Sitzmann et al.).")
        .def("GetMaxTransforms",
             &MolStandardize::TautomerEnumerator::getMaxTransforms,
             (python::arg("self")),
             "returns the maximum number of transformations to be applied.")
        .def("SetRemoveSp3Stereo",
             &MolStandardize::TautomerEnumerator::setRemoveSp3Stereo,
             (python::arg("self"), python::arg("removeSp3Stereo")),
             "set to True if you wish stereochemistry information "
             "to be removed from sp3 atoms involved in tautomerism. "
             "This means that S-aminoacids will lose their stereochemistry "
             "after going through tautomer enumeration because of the "
             "amido-imidol tautomerism. This defaults to True in RDKit, "
             "and to False in the workflow described by Sitzmann et al.")
        .def("GetRemoveSp3Stereo",
             &MolStandardize::TautomerEnumerator::getRemoveSp3Stereo,
             (python::arg("self")),
             "returns whether stereochemistry information will be removed from "
             "sp3 atoms involved in tautomerism.")
        .def("SetRemoveBondStereo",
             &MolStandardize::TautomerEnumerator::setRemoveBondStereo,
             (python::arg("self"), python::arg("removeBondStereo")),
             "set to True if you wish stereochemistry information "
             "to be removed from double bonds involved in tautomerism. "
             "This means that enols will lose their E/Z stereochemistry "
             "after going through tautomer enumeration because of the "
             "keto-enolic tautomerism. This defaults to True in the "
             "RDKit and also in the workflow described by Sitzmann et al.")
        .def("GetRemoveBondStereo",
             &MolStandardize::TautomerEnumerator::getRemoveBondStereo,
             (python::arg("self")),
             "returns whether stereochemistry information "
             "will be removed from double bonds involved in tautomerism.")
        .def("SetReassignStereo",
             &MolStandardize::TautomerEnumerator::setReassignStereo,
             (python::arg("self"), python::arg("reassignStereo")),
             "set to True if you wish AssignStereochemistry to be called "
             "on each tautomer generated by the Enumerate() method. "
             "This defaults to True.")
        .def("GetReassignStereo",
             &MolStandardize::TautomerEnumerator::getReassignStereo,
             (python::arg("self")),
             "returns whether AssignStereochemistry will be called "
             "on each tautomer generated by the Enumerate() method.")
        .def("SetCallback", &setCallbackHelper,
             "Pass an instance of a class derived from\n"
             "TautomerEnumeratorCallback, which must implement the\n"
             "__call__() method.")
        .def("GetCallback", &getCallbackHelper,
             "Get the TautomerEnumeratorCallback subclass instance,\n"
             "or None if none was set.")
        .def_readonly(
            "tautomerScoreVersion",
            MolStandardize::TautomerScoringFunctions::tautomerScoringVersion);
    python::def("GetV1TautomerEnumerator",
                MolStandardize::getV1TautomerEnumerator,
                "return a TautomerEnumerator using v1 of the enumeration rules",
                python::return_value_policy<python::manage_new_object>());
  }
};

void wrap_tautomer() { tautomer_wrapper::wrap(); }
