//
//  Copyright (C) 2020 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
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
      : d_tr(std::move(tr)) {
    python::list atList;
    python::list bndList;
    for (unsigned int i = 0; i < tr.modifiedAtoms.size(); ++i) {
      if (d_tr.modifiedAtoms.test(i)) {
        atList.append(i);
      }
    }
    for (unsigned int i = 0; i < tr.modifiedBonds.size(); ++i) {
      if (d_tr.modifiedBonds.test(i)) {
        bndList.append(i);
      }
    }
    d_atTuple = python::tuple(atList);
    d_bndTuple = python::tuple(bndList);
  }
  inline const std::vector<ROMOL_SPTR> &tautomers() { return d_tr.tautomers; }
  inline MolStandardize::TautomerEnumeratorStatus status() {
    return d_tr.status;
  }
  inline python::tuple modifiedAtoms() { return d_atTuple; }
  inline python::tuple modifiedBonds() { return d_bndTuple; }

 private:
  const MolStandardize::TautomerEnumeratorResult d_tr;
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
  bool operator()(const ROMol &mol,
                  const MolStandardize::TautomerEnumeratorResult &res) {
    PyTautomerEnumeratorResult pyRes(res);
    return getCallbackOverride()(boost::ref(mol), boost::ref(pyRes));
  }
  python::object getPyCallbackObject() { return d_pyCallbackObject; }

 private:
  PyTautomerEnumeratorCallback *d_cppCallback;
  python::object d_pyCallbackObject;
};

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

ROMol *canonicalizeHelper(const MolStandardize::TautomerEnumerator &self,
                          const ROMol &mol) {
  return self.canonicalize(mol);
}

class pyobjFunctor {
 public:
  pyobjFunctor(python::object obj) : dp_obj(std::move(obj)) {}
  ~pyobjFunctor() {}
  int operator()(const ROMol &m) {
    return python::extract<int>(dp_obj(boost::ref(m)));
  }

 private:
  python::object dp_obj;
};

ROMol *canonicalizeHelper2(const MolStandardize::TautomerEnumerator &self,
                           const ROMol &mol, python::object scoreFunc) {
  pyobjFunctor ftor(scoreFunc);
  return self.canonicalize(mol, ftor);
}

PyTautomerEnumeratorResult *enumerateHelper(
    const MolStandardize::TautomerEnumerator &self, const ROMol &mol,
    bool reassignStereo = true) {
  return new PyTautomerEnumeratorResult(self.enumerate(mol, reassignStereo));
}

}  // namespace

struct tautomer_wrapper {
  static void wrap() {
    python::enum_<MolStandardize::TautomerEnumeratorStatus>(
        "TautomerEnumeratorStatus")
        .value("Completed", MolStandardize::Completed)
        .value("MaxTautomersReached", MolStandardize::MaxTautomersReached)
        .value("MaxTransformsReached", MolStandardize::MaxTransformsReached)
        .value("Canceled", MolStandardize::Canceled);

    python::class_<PyTautomerEnumeratorCallback, boost::noncopyable>(
        "TautomerEnumeratorCallback",
        tautomerEnumeratorCallbackClassDoc.c_str(), python::init<>())
        .def("__call__",
             python::pure_virtual(&PyTautomerEnumeratorCallback::operator()),
             "This must be implemented in the derived class. "
             "Return True if the tautomer enumeration should continue; "
             "False if the tautomer enumeration should stop.\n");

    python::class_<PyTautomerEnumeratorResult, boost::noncopyable>(
        "TautomerEnumeratorResult",
        "used to return tautomer enumeration results", python::no_init)
        .add_property(
            "tautomers",
            python::make_function(&PyTautomerEnumeratorResult::tautomers,
                                  python::return_internal_reference<>()),
            "tautomers generated by the enumerator")
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
            "tuple of bond indices modified by the transforms");

    python::class_<MolStandardize::TautomerEnumerator, boost::noncopyable>(
        "TautomerEnumerator", python::no_init)
        .def("__init__", python::make_constructor(createDefaultEnumerator))
        .def("__init__", python::make_constructor(EnumeratorFromParams))
        .def("Enumerate", &enumerateHelper,
             (python::arg("self"), python::arg("mol"),
              python::arg("reassignStereo") = true),
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
  "acceptor").
  The reassignStereo parameter determines whether AssignStereochemistry is
  called on all generated tautomers or not (defaults to True).)DOC")
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
             "False in the workflow described by Sitzmann et al.")
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
        .def("SetCallback", &setCallbackHelper,
             "Pass an instance of a class derived from\n"
             "TautomerEnumeratorCallback, which must implement the\n"
             "__call__() method.\n")
        .def("GetCallback", &getCallbackHelper,
             "Get the TautomerEnumeratorCallback subclass instance,\n"
             "or None if none was set.\n")
        .def_readonly(
            "tautomerScoreVersion",
            MolStandardize::TautomerScoringFunctions::tautomerScoringVersion);
  }
};

void wrap_tautomer() { tautomer_wrapper::wrap(); }
