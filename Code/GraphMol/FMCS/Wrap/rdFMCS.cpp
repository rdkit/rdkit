//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/python.h>
#include <GraphMol/ROMol.h>
#include <RDBoost/Wrap.h>
#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/FMCS/RingMatchTableSet.h>

#define COMPARE_FUNC_NAME "compare"
#define CALLBACK_FUNC_NAME "callback"

namespace python = boost::python;

namespace RDKit {

struct PyMCSAtomCompare {
  inline bool checkAtomRingMatch(const MCSAtomCompareParameters& p,
                                 const ROMol& mol1, unsigned int atom1,
                                 const ROMol& mol2, unsigned int atom2) const {
    return RDKit::checkAtomRingMatch(p, mol1, atom1, mol2, atom2);
  }
  inline bool checkAtomCharge(const MCSAtomCompareParameters& p,
                              const ROMol& mol1, unsigned int atom1,
                              const ROMol& mol2, unsigned int atom2) const {
    return RDKit::checkAtomCharge(p, mol1, atom1, mol2, atom2);
  }
  inline bool checkAtomChirality(const MCSAtomCompareParameters& p,
                                 const ROMol& mol1, unsigned int atom1,
                                 const ROMol& mol2, unsigned int atom2) const {
    return RDKit::checkAtomChirality(p, mol1, atom1, mol2, atom2);
  }
  virtual bool compare(const MCSAtomCompareParameters&,
                       const ROMol&, unsigned int, const ROMol&, unsigned int) {
    PyErr_SetString(PyExc_AttributeError,
      "The " COMPARE_FUNC_NAME "() method "
      "must be overridden in the rdFMCS.MCSAtomCompare subclass");
    python::throw_error_already_set();
    return false;
  }
};

struct PyMCSBondCompare {
  inline bool checkBondStereo(const MCSBondCompareParameters& p,
                              const ROMol& mol1, unsigned int bond1,
                              const ROMol& mol2, unsigned int bond2) const {
    return RDKit::checkBondStereo(p, mol1, bond1, mol2, bond2);
  }
  inline bool checkBondRingMatch(const MCSBondCompareParameters& p,
                                 const ROMol& mol1, unsigned int bond1,
                                 const ROMol& mol2, unsigned int bond2) {
    if (!ringMatchTablesMols.count(&mol1)) {
      ringMatchTables.init(&mol1);
      ringMatchTables.computeRingMatchTable(&mol1, &mol1, *mcsParameters);
      ringMatchTablesMols.insert(&mol1);
    }
    if (!ringMatchTablesMols.count(&mol2)) {
      ringMatchTables.computeRingMatchTable(&mol1, &mol2, *mcsParameters);
      ringMatchTables.addTargetBondRingsIndeces(&mol2);
      ringMatchTablesMols.insert(&mol2);
    }
    return RDKit::checkBondRingMatch(p, mol1, bond1, mol2, bond2, &ringMatchTables);
  }
  virtual bool compare(const MCSBondCompareParameters&,
                       const ROMol&, unsigned int, const ROMol&, unsigned int) {
    PyErr_SetString(PyExc_AttributeError,
      "The " COMPARE_FUNC_NAME "() method "
      "must be overridden in the rdFMCS.MCSBondCompare subclass");
    python::throw_error_already_set();
    return false;
  };
  const MCSParameters *mcsParameters;
  std::set<const ROMol *> ringMatchTablesMols;
  FMCS::RingMatchTableSet ringMatchTables;
};

struct PyCompareFunctionUserData {
  const MCSParameters *mcsParameters;
  std::set<const ROMol *> *ringMatchTablesMols;
  FMCS::RingMatchTableSet *ringMatchTables;
  python::object pyAtomComp;
  python::object pyBondComp;
};

struct PyProgressCallbackUserData {
  const MCSProgressData *mcsProgressData;
  python::object pyMCSProgress;
  python::object pyAtomComp;
  python::object pyBondComp;
};

struct PyMCSProgress {
  virtual bool callback(const MCSProgressData&,
                        const MCSParameters&) {
    PyErr_SetString(PyExc_AttributeError,
      "The " CALLBACK_FUNC_NAME "() method "
      "must be overridden in the rdFMCS.MCSProgress subclass");
    python::throw_error_already_set();
    return false;
  }
};

class PyMCSProgressData {
public:
  PyMCSProgressData() :
    pd(new MCSProgressData()),
    pcud(new PyProgressCallbackUserData()) {
      pcud->mcsProgressData = pd.get();
    }
  PyMCSProgressData(const MCSProgressData &other) :
    PyMCSProgressData() {
      *pd = other;
    }
  unsigned int getNumAtoms() const {
    return pd->NumAtoms;
  }
  unsigned int getNumBonds() const {
    return pd->NumBonds;
  }
  unsigned int getSeedProcessed() const {
    return pd->SeedProcessed;
  }
private:
  std::unique_ptr<MCSProgressData> pd;
  std::unique_ptr<PyProgressCallbackUserData> pcud;
};

class PyMCSParameters {
public:
  PyMCSParameters() :
    p(new MCSParameters()),
    cfud(new PyCompareFunctionUserData()),
    pcud(new PyProgressCallbackUserData()) {
      cfud->mcsParameters = p.get();
      pcud->mcsProgressData = nullptr;
    }
  PyMCSParameters(const MCSParameters &other,
                  const PyProgressCallbackUserData &pcudOther) :
    PyMCSParameters() {
      *p = other;
      pcud->pyMCSProgress = pcudOther.pyMCSProgress;
      cfud->pyAtomComp = pcudOther.pyAtomComp;
      cfud->pyBondComp = pcudOther.pyBondComp;
    }
  const MCSParameters *get() const {
    return p.get();
  }
  bool getMaximizeBonds() const {
    return p->MaximizeBonds;
  }
  void setMaximizeBonds(bool value) {
    p->MaximizeBonds = value;
  }
  double getThreshold() const {
    return p->Threshold;
  }
  void setThreshold(double value) {
    p->Threshold = value;
  }
  unsigned int getTimeout() const {
    return p->Timeout;
  }
  void setTimeout(unsigned int value) {
    p->Timeout = value;
  }
  bool getVerbose() const {
    return p->Verbose;
  }
  void setVerbose(bool value) {
    p->Verbose = value;
  }
  const MCSAtomCompareParameters& getAtomCompareParameters() const {
    return p->AtomCompareParameters;
  }
  void setAtomCompareParameters(const MCSAtomCompareParameters &value) {
    p->AtomCompareParameters = value;
  }
  const MCSBondCompareParameters& getBondCompareParameters() const {
    return p->BondCompareParameters;
  }
  void setBondCompareParameters(const MCSBondCompareParameters &value) {
    p->BondCompareParameters = value;
  }
  std::string getInitialSeed() const {
    return p->InitialSeed;
  }
  void setInitialSeed(const std::string &value) {
    p->InitialSeed = value;
  }
  void setMCSAtomTyper(PyObject *atomComp) {
    PRECONDITION(atomComp, "atomComp must not be NULL");
    python::object atomCompObject(python::handle<>(python::borrowed(atomComp)));
    python::extract<AtomComparator> extractAtomComparator(atomCompObject);
    if (extractAtomComparator.check()) {
      AtomComparator ac(extractAtomComparator());
      p->setMCSAtomTyperFromEnum(ac);
    }
    else {
      python::extract<PyMCSAtomCompare *> extractPyMCSAtomCompare(atomCompObject);
      if (extractPyMCSAtomCompare.check()) {
        PyObject *callable = PyObject_GetAttrString(atomCompObject.ptr(), COMPARE_FUNC_NAME);
        if (!callable) {
          // should never get here
          PyErr_SetString(PyExc_AttributeError,
            "The " COMPARE_FUNC_NAME "() method must be defined "
            "in the rdFMCS.MCSAtomCompare subclass");
          python::throw_error_already_set();
        }
        if (!PyCallable_Check(callable)) {
          PyErr_SetString(PyExc_TypeError,
            "The " COMPARE_FUNC_NAME " attribute in the "
            "rdFMCS.MCSAtomCompare subclass is not a callable method");
          python::throw_error_already_set();
        }
        p->CompareFunctionsUserData = cfud.get();
        p->AtomTyper = MCSAtomComparePyFunc;
        cfud->pyAtomComp = atomCompObject;
        cfud->mcsParameters = p.get();
      }
      else {
        PyErr_SetString(PyExc_TypeError,
          "expected an instance of a rdFMCS.MCSAtomCompare subclass "
          "or a member of the AtomCompare class");
        python::throw_error_already_set();
      }
    }
  }
  python::object getMCSAtomTyper() {
    static const std::map<RDKit::MCSAtomCompareFunction,
                          RDKit::AtomComparator> atomTyperToComp = {
      { MCSAtomCompareAny, AtomCompareAny },
      { MCSAtomCompareElements, AtomCompareElements },
      { MCSAtomCompareIsotopes, AtomCompareIsotopes },
      { MCSAtomCompareAnyHeavyAtom, AtomCompareAnyHeavyAtom }
    };
    if (!cfud->pyAtomComp.is_none()) {
      return cfud->pyAtomComp;
    }
    python::object res;
    try {
      res = python::object(atomTyperToComp.at(p->AtomTyper));
    }
    catch(const std::out_of_range&) {
      PyErr_SetString(PyExc_TypeError,
        "Unknown AtomTyper");
      python::throw_error_already_set();
    }
    return res;
  }
  void setMCSBondTyper(PyObject *bondComp) {
    PRECONDITION(bondComp, "bondComp must not be NULL");
    python::object bondCompObject(python::handle<>(python::borrowed(bondComp)));
    python::extract<BondComparator> extractBondComparator(bondCompObject);
    if (extractBondComparator.check()) {
      BondComparator bc(extractBondComparator());
      p->setMCSBondTyperFromEnum(bc);
    }
    else {
      python::extract<PyMCSBondCompare *> extractPyMCSBondCompare(bondCompObject);
      if (extractPyMCSBondCompare.check()) {
        PyObject *callable = PyObject_GetAttrString(bondCompObject.ptr(), COMPARE_FUNC_NAME);
        if (!callable) {
          // should never get here
          PyErr_SetString(PyExc_AttributeError,
            "The " COMPARE_FUNC_NAME "() method must be defined "
            "in the rdFMCS.MCSBondCompare subclass");
          python::throw_error_already_set();
        }
        if (!PyCallable_Check(callable)) {
          PyErr_SetString(PyExc_TypeError,
            "The " COMPARE_FUNC_NAME " attribute in the "
            "rdFMCS.MCSBondCompare subclass is not a callable method");
          python::throw_error_already_set();
        }
        p->CompareFunctionsUserData = cfud.get();
        p->BondTyper = MCSBondComparePyFunc;
        cfud->pyBondComp = bondCompObject;
        PyMCSBondCompare *bc = extractPyMCSBondCompare();
        bc->mcsParameters = p.get();
        cfud->mcsParameters = p.get();
        cfud->ringMatchTablesMols = &bc->ringMatchTablesMols;
        cfud->ringMatchTables = &bc->ringMatchTables;
      }
      else {
        PyErr_SetString(PyExc_TypeError,
          "expected an instance of a rdFMCS.MCSBondCompare subclass "
          "or a member of the BondCompare class");
        python::throw_error_already_set();
      }
    }
  }
  python::object getMCSBondTyper() {
    static const std::map<RDKit::MCSBondCompareFunction,
                          RDKit::BondComparator> bondTyperToComp = {
      { MCSBondCompareAny, BondCompareAny },
      { MCSBondCompareOrder, BondCompareOrder },
      { MCSBondCompareOrderExact, BondCompareOrderExact }
    };
    if (!cfud->pyBondComp.is_none()) {
      return cfud->pyBondComp;
    }
    python::object res;
    try {
      res = python::object(bondTyperToComp.at(p->BondTyper));
    }
    catch(const std::out_of_range&) {
      PyErr_SetString(PyExc_TypeError,
        "Unknown BondTyper");
      python::throw_error_already_set();
    }
    return res;
  }
  void setMCSProgressCallback(PyObject *progress) {
    PRECONDITION(progress, "progress must not be NULL");
    python::object progressObject(python::handle<>(python::borrowed(progress)));
    python::extract<PyMCSProgress *> extractMCSProgress(progressObject);
    if (extractMCSProgress.check()) {
      p->ProgressCallbackUserData = pcud.get();
      p->ProgressCallback = MCSProgressCallbackPyFunc;
      pcud->pyMCSProgress = progressObject;
      pcud->pyAtomComp = cfud->pyAtomComp;
      pcud->pyBondComp = cfud->pyBondComp;
    }
    else {
      PyErr_SetString(PyExc_TypeError,
        "expected an instance of a rdFMCS.MCSProgress subclass");
      python::throw_error_already_set();
    }
  }
  python::object getMCSProgressCallback() {
    if (!pcud->pyMCSProgress.is_none()) {
      return pcud->pyMCSProgress;
    }
    return python::object();
  }
  void clearRingTableCache() {
    if (cfud->ringMatchTablesMols) {
      cfud->ringMatchTablesMols->clear();
    }
    if (cfud->ringMatchTables) {
      cfud->ringMatchTables->clear();
    }
  }
private:
  static bool MCSAtomComparePyFunc(const MCSAtomCompareParameters& p,
                                   const ROMol& mol1, unsigned int atom1,
                                   const ROMol& mol2, unsigned int atom2,
                                   void* userData) {
    PRECONDITION(userData, "userData must not be NULL");
    PyCompareFunctionUserData *cfud = static_cast<PyCompareFunctionUserData *>(userData);
    bool res = false;
    {
      PyGILStateHolder h;
      res = python::call_method<bool>(cfud->pyAtomComp.ptr(), COMPARE_FUNC_NAME,
                                      boost::ref(p), boost::ref(mol1),
                                      atom1, boost::ref(mol2), atom2);
    }
    return res;
  }
  static bool MCSBondComparePyFunc(const MCSBondCompareParameters& p,
                                   const ROMol& mol1, unsigned int bond1,
                                   const ROMol& mol2, unsigned int bond2,
                                   void* userData) {
    PRECONDITION(userData, "userData must not be NULL");
    PyCompareFunctionUserData *cfud = static_cast<PyCompareFunctionUserData *>(userData);
    bool res = false;
    {
      PyGILStateHolder h;
      res = python::call_method<bool>(cfud->pyBondComp.ptr(), COMPARE_FUNC_NAME,
                                      boost::ref(p), boost::ref(mol1),
                                      bond1, boost::ref(mol2), bond2);
    }
    return res;
  }
  static bool MCSProgressCallbackPyFunc(const MCSProgressData& stat,
                                        const MCSParameters& params, void* userData) {
    PRECONDITION(userData, "userData must not be NULL");
    PyProgressCallbackUserData *pcud = static_cast<PyProgressCallbackUserData *>(userData);
    bool res = false;
    {
      PyMCSParameters ps(params, *pcud);
      PyMCSProgressData pd(stat);
      PyGILStateHolder h;
      res = python::call_method<bool>(pcud->pyMCSProgress.ptr(), CALLBACK_FUNC_NAME,
                                      boost::ref(pd), boost::ref(ps));
    }
    return res;
  }
  std::unique_ptr<MCSParameters> p;
  std::unique_ptr<PyCompareFunctionUserData> cfud;
  std::unique_ptr<PyProgressCallbackUserData> pcud;
};

MCSResult *FindMCSWrapper(python::object mols, bool maximizeBonds,
                          double threshold, unsigned timeout, bool verbose,
                          bool matchValences, bool ringMatchesRingOnly,
                          bool completeRingsOnly, bool matchChiralTag,
                          AtomComparator atomComp, BondComparator bondComp,
                          RingComparator ringComp, std::string seedSmarts) {
  std::vector<ROMOL_SPTR> ms;
  unsigned int nElems = python::extract<unsigned int>(mols.attr("__len__")());
  ms.resize(nElems);
  for (unsigned int i = 0; i < nElems; ++i) {
    if (!mols[i]) {
      throw_value_error("molecule is None");
    }
    ms[i] = python::extract<ROMOL_SPTR>(mols[i]);
  }
  MCSParameters p;
  p.Threshold = threshold;
  p.MaximizeBonds = maximizeBonds;
  p.Timeout = timeout;
  p.Verbose = verbose;
  p.InitialSeed = seedSmarts;
  p.AtomCompareParameters.MatchValences = matchValences;
  p.AtomCompareParameters.MatchChiralTag = matchChiralTag;
  p.AtomCompareParameters.RingMatchesRingOnly = ringMatchesRingOnly;
  p.setMCSAtomTyperFromEnum(atomComp);
  p.setMCSBondTyperFromEnum(bondComp);
  p.BondCompareParameters.RingMatchesRingOnly = ringMatchesRingOnly;
  p.BondCompareParameters.CompleteRingsOnly = completeRingsOnly;
  p.BondCompareParameters.MatchFusedRings = (ringComp != IgnoreRingFusion);
  p.BondCompareParameters.MatchFusedRingsStrict = (ringComp == StrictRingFusion);

  MCSResult *res = nullptr;
  {
    NOGIL gil;
    res = new MCSResult(findMCS(ms, &p));
  }
  return res;
}

MCSResult *FindMCSWrapper2(python::object mols, PyMCSParameters &pyMcsParams) {
  std::vector<ROMOL_SPTR> ms;
  unsigned int nElems = python::extract<unsigned int>(mols.attr("__len__")());
  ms.resize(nElems);
  for (unsigned int i = 0; i < nElems; ++i) {
    if (!mols[i]) {
      throw_value_error("molecule is None");
    }
    ms[i] = python::extract<ROMOL_SPTR>(mols[i]);
  }

  MCSResult *res = nullptr;
  pyMcsParams.clearRingTableCache();
  {
    NOGIL gil;
    res = new MCSResult(findMCS(ms, pyMcsParams.get()));
  }
  return res;
}
}  // namespace RDKit
namespace {
struct mcsresult_wrapper {
  static void wrap() {
    python::class_<RDKit::MCSResult>("MCSResult", "used to return MCS results",
                                     python::no_init)
        .def_readonly("numAtoms", &RDKit::MCSResult::NumAtoms,
                      "number of atoms in MCS")
        .def_readonly("numBonds", &RDKit::MCSResult::NumBonds,
                      "number of bonds in MCS")
        .def_readonly("queryMol", &RDKit::MCSResult::QueryMol,
                      "query molecule for the MCS")
        .def_readonly("smartsString", &RDKit::MCSResult::SmartsString,
                      "SMARTS string for the MCS")
        .def_readonly("canceled", &RDKit::MCSResult::Canceled,
                      "if True, the MCS calculation did not finish");
  }
};
}  // namespace

BOOST_PYTHON_MODULE(rdFMCS) {
  python::scope().attr("__doc__") =
      "Module containing a C++ implementation of the FMCS algorithm";
  mcsresult_wrapper::wrap();

  python::enum_<RDKit::AtomComparator>("AtomCompare")
      .value("CompareAny", RDKit::AtomCompareAny)
      .value("CompareElements", RDKit::AtomCompareElements)
      .value("CompareIsotopes", RDKit::AtomCompareIsotopes)
      .value("CompareAnyHeavyAtom", RDKit::AtomCompareAnyHeavyAtom);
  python::enum_<RDKit::BondComparator>("BondCompare")
      .value("CompareAny", RDKit::BondCompareAny)
      .value("CompareOrder", RDKit::BondCompareOrder)
      .value("CompareOrderExact", RDKit::BondCompareOrderExact);
  python::enum_<RDKit::RingComparator>("RingCompare")
      .value("IgnoreRingFusion", RDKit::IgnoreRingFusion)
      .value("PermissiveRingFusion", RDKit::PermissiveRingFusion)
      .value("StrictRingFusion", RDKit::StrictRingFusion);

  std::string docString = "Find the MCS for a set of molecules";
  python::def(
      "FindMCS", RDKit::FindMCSWrapper,
      (python::arg("mols"), python::arg("maximizeBonds") = true,
       python::arg("threshold") = 1.0, python::arg("timeout") = 3600,
       python::arg("verbose") = false, python::arg("matchValences") = false,
       python::arg("ringMatchesRingOnly") = false,
       python::arg("completeRingsOnly") = false,
       python::arg("matchChiralTag") = false,
       python::arg("atomCompare") = RDKit::AtomCompareElements,
       python::arg("bondCompare") = RDKit::BondCompareOrder,
       python::arg("ringCompare") = RDKit::IgnoreRingFusion,
       python::arg("seedSmarts") = ""),
      python::return_value_policy<python::manage_new_object>(),
      docString.c_str());

  python::class_<RDKit::PyMCSParameters, boost::noncopyable>(
      "MCSParameters", "Parameters controlling how the MCS is constructed")
      .add_property("MaximizeBonds",
                    &RDKit::PyMCSParameters::getMaximizeBonds,
                    &RDKit::PyMCSParameters::setMaximizeBonds,
                    "toggles maximizing the number of bonds (instead of the "
                    "number of atoms)")
      .add_property("Threshold",
                    &RDKit::PyMCSParameters::getThreshold,
                    &RDKit::PyMCSParameters::setThreshold,
                    "fraction of the dataset that must contain the MCS")
      .add_property("Timeout",
                    &RDKit::PyMCSParameters::getTimeout,
                    &RDKit::PyMCSParameters::setTimeout,
                    "timeout (in seconds) for the calculation")
      .add_property("Verbose",
                    &RDKit::PyMCSParameters::getVerbose,
                    &RDKit::PyMCSParameters::setVerbose,
                    "toggles verbose mode")
      .add_property("AtomCompareParameters", python::make_function(
                    &RDKit::PyMCSParameters::getAtomCompareParameters,
                    python::return_internal_reference<>()),
                    &RDKit::PyMCSParameters::setAtomCompareParameters,
                    "parameters for comparing atoms")
      .add_property("BondCompareParameters", python::make_function(
                    &RDKit::PyMCSParameters::getBondCompareParameters,
                    python::return_internal_reference<>()),
                    &RDKit::PyMCSParameters::setBondCompareParameters,
                    "parameters for comparing bonds")
      .add_property("AtomTyper",
                    &RDKit::PyMCSParameters::getMCSAtomTyper,
                    &RDKit::PyMCSParameters::setMCSAtomTyper,
                    "atom typer to be used. Must be one of the "
                    "members of the rdFMCS.AtomCompare class or "
                    "an instance of a user-defined subclass of "
                    "rdFMCS.MCSAtomCompare")
      .add_property("BondTyper",
                    &RDKit::PyMCSParameters::getMCSBondTyper,
                    &RDKit::PyMCSParameters::setMCSBondTyper,
                    "bond typer to be used. Must be one of the "
                    "members of the rdFMCS.BondCompare class or "
                    "an instance of a user-defined subclass of "
                    "rdFMCS.MCSBondCompare")
      .add_property("ProgressCallback",
                    &RDKit::PyMCSParameters::getMCSProgressCallback,
                    &RDKit::PyMCSParameters::setMCSProgressCallback,
                    "progress callback class. Must be a "
                    "user-defined subclass of rdFMCS.Progress")
      .add_property("InitialSeed",
                    &RDKit::PyMCSParameters::getInitialSeed,
                    &RDKit::PyMCSParameters::setInitialSeed,
                    "SMILES string to be used as the seed of the MCS")
      .def("SetAtomTyper", &RDKit::PyMCSParameters::setMCSAtomTyper,
           (python::arg("self"), python::arg("comparator")),
           "DEPRECATED: please use the AtomTyper property instead. "
           "Sets the atom typer to be used. The argument must be one "
           "of the members of the rdFMCS.MCSAtomCompare class or an "
           "instance of a user-defined subclass of rdFMCS.MCSAtomCompare")
      .def("SetBondTyper", &RDKit::PyMCSParameters::setMCSBondTyper,
           (python::arg("self"), python::arg("comparator")),
           "DEPRECATED: please use the BondTyper property instead. "
           "Sets the bond typer to be used. The argument must be one "
           "of the members of the rdFMCS.MCSBondCompare class or an "
           "instance of a user-defined subclass of rdFMCS.MCSBondCompare");

  python::class_<RDKit::MCSAtomCompareParameters, boost::noncopyable>(
      "MCSAtomCompareParameters",
      "Parameters controlling how atom-atom matching is done")
      .def_readwrite("MatchValences",
                     &RDKit::MCSAtomCompareParameters::MatchValences,
                     "include atom valences in the match")
      .def_readwrite("MatchChiralTag",
                     &RDKit::MCSAtomCompareParameters::MatchChiralTag,
                     "include atom chirality in the match")
      .def_readwrite("MatchFormalCharge",
                     &RDKit::MCSAtomCompareParameters::MatchFormalCharge,
                     "include formal charge in the match")
      .def_readwrite("RingMatchesRingOnly",
                     &RDKit::MCSAtomCompareParameters::RingMatchesRingOnly,
                     "ring atoms are only allowed to match other ring atoms")
      .def_readwrite("MatchIsotope",
                     &RDKit::MCSAtomCompareParameters::MatchIsotope,
                     "use isotope atom queries in MCSResults");

  python::class_<RDKit::MCSBondCompareParameters, boost::noncopyable>(
      "MCSBondCompareParameters",
      "Parameters controlling how bond-bond matching is done")
      .def_readwrite("RingMatchesRingOnly",
                     &RDKit::MCSBondCompareParameters::RingMatchesRingOnly,
                     "ring bonds are only allowed to match other ring bonds")
      .def_readwrite("CompleteRingsOnly",
                     &RDKit::MCSBondCompareParameters::CompleteRingsOnly,
                     "results cannot include partial rings")
      .def_readwrite("MatchFusedRings",
                     &RDKit::MCSBondCompareParameters::MatchFusedRings,
                     "enforce check on ring fusion, i.e. alpha-methylnaphthalene "
                     "won't match beta-methylnaphtalene, but decalin "
                     "will match cyclodecane unless MatchFusedRingsStrict is True")
      .def_readwrite("MatchFusedRingsStrict",
                     &RDKit::MCSBondCompareParameters::MatchFusedRingsStrict,
                     "only enforced if MatchFusedRings is True; the ring fusion "
                     "must be the same in both query and target, i.e. decalin "
                     "won't match cyclodecane")
      .def_readwrite("MatchStereo",
                     &RDKit::MCSBondCompareParameters::MatchStereo,
                     "include bond stereo in the comparison");

  python::class_<RDKit::PyMCSProgress, boost::noncopyable>(
      "MCSProgress", "Base class. Subclass and override "
      "MCSProgress." CALLBACK_FUNC_NAME "callback() "
      "to define a custom callback function")
      .def(CALLBACK_FUNC_NAME, &RDKit::PyMCSProgress::callback,
           (python::arg("self"), python::arg("stat"),
              python::arg("parameters")),
              "override to implement a custom progress callback");

  python::class_<RDKit::PyMCSProgressData, boost::noncopyable>(
      "MCSProgressData", "Information about the MCS progress")
      .add_property("numAtoms", &RDKit::PyMCSProgressData::getNumAtoms,
                    "number of atoms in MCS")
      .add_property("numBonds", &RDKit::PyMCSProgressData::getNumBonds,
                    "number of bonds in MCS")
      .add_property("seedProcessed", &RDKit::PyMCSProgressData::getSeedProcessed,
                    "number of processed seeds");

  python::class_<RDKit::PyMCSAtomCompare, boost::noncopyable>(
      "MCSAtomCompare", "Base class. Subclass and override "
      "MCSAtomCompare." COMPARE_FUNC_NAME "() to define custom "
      "atom compare functions, then set MCSParameters.AtomTyper "
      "to an instance of the subclass")
      .def("CheckAtomRingMatch", &RDKit::PyMCSAtomCompare::checkAtomRingMatch,
           (python::arg("self"), python::arg("parameters"),
              python::arg("mol1"), python::arg("atom1"),
              python::arg("mol2"), python::arg("atom2")),
              "Return True if both atoms are, or are not, in a ring")
      .def("CheckAtomCharge", &RDKit::PyMCSAtomCompare::checkAtomCharge,
           (python::arg("self"), python::arg("parameters"),
              python::arg("mol1"), python::arg("atom1"),
              python::arg("mol2"), python::arg("atom2")),
              "Return True if both atoms have the same formal charge")
      .def("CheckAtomChirality", &RDKit::PyMCSAtomCompare::checkAtomChirality,
           (python::arg("self"), python::arg("parameters"),
              python::arg("mol1"), python::arg("atom1"),
              python::arg("mol2"), python::arg("atom2")),
              "Return True if both atoms have, or have not, a chiral tag")
      .def(COMPARE_FUNC_NAME, &RDKit::PyMCSAtomCompare::compare,
           (python::arg("self"), python::arg("parameters"),
              python::arg("mol1"), python::arg("atom1"),
              python::arg("mol2"), python::arg("atom2")),
              "override to implement custom atom comparison");

  python::class_<RDKit::PyMCSBondCompare, boost::noncopyable>(
      "MCSBondCompare", "Base class. Subclass and override "
      "MCSBondCompare." COMPARE_FUNC_NAME "() to define custom "
      "bond compare functions, then set MCSParameters.BondTyper "
      "to an instance of the subclass")
      .def("CheckBondStereo", &RDKit::PyMCSBondCompare::checkBondStereo,
           (python::arg("self"), python::arg("parameters"),
              python::arg("mol1"), python::arg("bond1"),
              python::arg("mol2"), python::arg("bond2")),
              "Return True if both bonds have, or have not, a stereo descriptor")
      .def("CheckBondRingMatch", &RDKit::PyMCSBondCompare::checkBondRingMatch,
           (python::arg("self"), python::arg("parameters"),
              python::arg("mol1"), python::arg("bond1"),
              python::arg("mol2"), python::arg("bond2")),
              "Return True if both bonds are, or are not, part of a ring")
      .def(COMPARE_FUNC_NAME, &RDKit::PyMCSBondCompare::compare,
           (python::arg("self"), python::arg("parameters"),
              python::arg("mol1"), python::arg("bond1"),
              python::arg("mol2"), python::arg("bond2")),
              "override to implement custom bond comparison");

  python::def("FindMCS", RDKit::FindMCSWrapper2,
              (python::arg("mols"), python::arg("parameters")),
              python::return_value_policy<python::manage_new_object>(),
              docString.c_str());
}
