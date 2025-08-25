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
#include <GraphMol/Wrap/substructmethods.h>
#include <RDBoost/Wrap.h>
#include <GraphMol/FMCS/FMCS.h>

#define COMPARE_FUNC_NAME "__call__"
#define CALLBACK_FUNC_NAME "__call__"

namespace python = boost::python;

namespace RDKit {

struct PyMCSWrapper : public boost::python::wrapper<PyMCSWrapper> {
 public:
  PyMCSWrapper() {}
  PyMCSWrapper(PyObject *obj) {
    PRECONDITION(obj, "PyObject* must not be NULL");
    d_pyObject.reset(
        new python::object(python::handle<>(python::borrowed(obj))));
  }
  virtual ~PyMCSWrapper() {}
  virtual const char *subclassName() const {
    throw std::invalid_argument(
        "subclassName() must be overridden in the "
        "derived class");
  }
  const python::object &pyObject() const { return *d_pyObject; }

 protected:
  void failedToExtractPyObject() const {
    std::stringstream ss;
    ss << "Failed to extract object from " << subclassName() << " subclass";
    PyErr_SetString(PyExc_RuntimeError, ss.str().c_str());
    python::throw_error_already_set();
  }
  void errorNotOverridden() const {
    std::stringstream ss;
    ss << "The " COMPARE_FUNC_NAME "() method must be overridden in the rdFMCS."
       << subclassName() << " subclass";
    PyErr_SetString(PyExc_AttributeError, ss.str().c_str());
    python::throw_error_already_set();
  }
  void errorNotDefined() const {
    // should never happen as the method is virtual but not pure in the C++
    // class
    std::stringstream ss;
    ss << "The " CALLBACK_FUNC_NAME
          "() method must be defined "
          "in the "
       << subclassName() << " subclass";
    PyErr_SetString(PyExc_AttributeError, ss.str().c_str());
    python::throw_error_already_set();
  }
  void errorNotCallable() const {
    std::stringstream ss;
    ss << "The " COMPARE_FUNC_NAME " attribute in the " << subclassName()
       << " subclass is not a callable method";
    PyErr_SetString(PyExc_TypeError, ss.str().c_str());
    python::throw_error_already_set();
  }
  virtual bool hasPythonOverride(const char *attrName) const {
    auto obj = get_override(attrName);
    return PyCallable_Check(obj.ptr());
  }
  void extractPyMCSWrapper() {
    d_pyObjectExtractor.reset(new python::extract<PyMCSWrapper *>(*d_pyObject));
    if (d_pyObjectExtractor->check()) {
      PyObject *callable =
          PyObject_GetAttrString(d_pyObject->ptr(), CALLBACK_FUNC_NAME);
      if (!callable) {
        errorNotDefined();
      }
      if (!PyCallable_Check(callable)) {
        errorNotCallable();
      }
      if (!pyObjectExtract()->hasPythonOverride(CALLBACK_FUNC_NAME)) {
        errorNotOverridden();
      }
    } else {
      std::stringstream ss;
      ss << "expected an instance of the rdFMCS." << subclassName()
         << " subclass";
      PyErr_SetString(PyExc_TypeError, ss.str().c_str());
      python::throw_error_already_set();
    }
  }
  PyMCSWrapper *pyObjectExtract() const { return (*d_pyObjectExtractor)(); }

 private:
  std::unique_ptr<python::object> d_pyObject;
  std::unique_ptr<python::extract<PyMCSWrapper *>> d_pyObjectExtractor;
};

struct PyMCSAtomCompare : public PyMCSWrapper {
 public:
  PyMCSAtomCompare() {}
  PyMCSAtomCompare(PyObject *obj) : PyMCSWrapper(obj) {}
  ~PyMCSAtomCompare() override {}
  bool extractAtomComparator(AtomComparator &ac) {
    bool res = false;
    python::extract<AtomComparator> predefinedAtomComparator(pyObject());
    if (predefinedAtomComparator.check()) {
      ac = predefinedAtomComparator();
      res = true;
    } else {
      extractPyMCSWrapper();
    }
    return res;
  }
  PyMCSAtomCompare *extractPyObject() const {
    auto res = dynamic_cast<PyMCSAtomCompare *>(pyObjectExtract());
    if (!res) {
      failedToExtractPyObject();
    }
    return res;
  }
  inline bool checkAtomRingMatch(const MCSAtomCompareParameters &p,
                                 const ROMol &mol1, unsigned int atom1,
                                 const ROMol &mol2, unsigned int atom2) const {
    return RDKit::checkAtomRingMatch(p, mol1, atom1, mol2, atom2);
  }
  inline bool checkAtomCharge(const MCSAtomCompareParameters &p,
                              const ROMol &mol1, unsigned int atom1,
                              const ROMol &mol2, unsigned int atom2) const {
    return RDKit::checkAtomCharge(p, mol1, atom1, mol2, atom2);
  }
  inline bool checkAtomChirality(const MCSAtomCompareParameters &p,
                                 const ROMol &mol1, unsigned int atom1,
                                 const ROMol &mol2, unsigned int atom2) const {
    return RDKit::checkAtomChirality(p, mol1, atom1, mol2, atom2);
  }
  inline const char *subclassName() const override { return "MCSAtomCompare"; }
  virtual bool operator()(const MCSAtomCompareParameters &, const ROMol &,
                          unsigned int, const ROMol &, unsigned int) const {
    errorNotOverridden();
    return false;
  }
};

struct PyMCSBondCompare : public PyMCSWrapper {
  PyMCSBondCompare() {}
  PyMCSBondCompare(PyObject *obj) : PyMCSWrapper(obj) {}
  ~PyMCSBondCompare() override {}
  bool extractBondComparator(BondComparator &bc) {
    bool res = false;
    python::extract<BondComparator> predefinedBondComparator(pyObject());
    if (predefinedBondComparator.check()) {
      bc = predefinedBondComparator();
      res = true;
    } else {
      extractPyMCSWrapper();
    }
    return res;
  }
  PyMCSBondCompare *extractPyObject() const {
    auto res = dynamic_cast<PyMCSBondCompare *>(pyObjectExtract());
    if (!res) {
      failedToExtractPyObject();
    }
    return res;
  }
  inline bool checkBondStereo(const MCSBondCompareParameters &p,
                              const ROMol &mol1, unsigned int bond1,
                              const ROMol &mol2, unsigned int bond2) const {
    return RDKit::checkBondStereo(p, mol1, bond1, mol2, bond2);
  }
  inline bool checkBondRingMatch(const MCSBondCompareParameters &p,
                                 const ROMol &mol1, unsigned int bond1,
                                 const ROMol &mol2, unsigned int bond2) {
    return RDKit::checkBondRingMatch(p, mol1, bond1, mol2, bond2);
  }
  inline const char *subclassName() const override { return "MCSBondCompare"; }
  virtual bool operator()(const MCSBondCompareParameters &, const ROMol &,
                          unsigned int, const ROMol &, unsigned int) const {
    errorNotOverridden();
    return false;
  }
  const MCSParameters *mcsParameters;
};

struct PyAtomBondCompData {
  python::object pyAtomComp;
  python::object pyBondComp;
  MCSBondCompareFunction standardBondTyperFunc = nullptr;
};

struct PyBaseUserData {
  PyAtomBondCompData pyAtomBondCompData;
};

struct PyCompareFunctionUserData : public PyBaseUserData {
  const MCSParameters *mcsParameters;
};

struct PyProgressCallbackUserData : public PyBaseUserData {
  const MCSProgressData *mcsProgressData;
  python::object pyMCSProgress;
};

struct PyMCSFinalMatchCheckFunctionUserData : public PyBaseUserData {
  python::object pyMCSFinalMatchCheck;
};

struct PyMCSAcceptanceFunctionUserData : public PyBaseUserData {
  python::object pyMCSAcceptance;
};

struct PyMCSProgress : public PyMCSWrapper {
  PyMCSProgress() {}
  PyMCSProgress(PyObject *obj) : PyMCSWrapper(obj) { extractPyMCSWrapper(); }
  ~PyMCSProgress() override {}
  PyMCSProgress *extractPyObject() const {
    auto res = dynamic_cast<PyMCSProgress *>(pyObjectExtract());
    if (!res) {
      failedToExtractPyObject();
    }
    return res;
  }
  inline const char *subclassName() const override { return "MCSProgress"; }
  virtual bool operator()(const MCSProgressData &,
                          const MCSParameters &) const {
    errorNotOverridden();
    return false;
  }
};

class PyMCSProgressData {
 public:
  PyMCSProgressData()
      : pd(new MCSProgressData()), pcud(new PyProgressCallbackUserData()) {
    pcud->mcsProgressData = pd.get();
  }
  PyMCSProgressData(const MCSProgressData &other) : PyMCSProgressData() {
    *pd = other;
  }
  unsigned int getNumAtoms() const { return pd->NumAtoms; }
  unsigned int getNumBonds() const { return pd->NumBonds; }
  unsigned int getSeedProcessed() const { return pd->SeedProcessed; }

 private:
  std::unique_ptr<MCSProgressData> pd;
  std::unique_ptr<PyProgressCallbackUserData> pcud;
};

struct PyMCSFinalMatchCheck : public PyMCSWrapper {
  PyMCSFinalMatchCheck() {}
  PyMCSFinalMatchCheck(PyObject *obj) : PyMCSWrapper(obj) {
    extractPyMCSWrapper();
  }
  ~PyMCSFinalMatchCheck() override {}
  PyMCSFinalMatchCheck *extractPyObject() const {
    auto res = dynamic_cast<PyMCSFinalMatchCheck *>(pyObjectExtract());
    if (!res) {
      failedToExtractPyObject();
    }
    return res;
  }
  inline const char *subclassName() const override {
    return "MCSFinalMatchCheck";
  }
  virtual bool operator()() const {
    errorNotOverridden();
    return false;
  }
};

struct PyMCSAcceptance : public PyMCSWrapper {
  PyMCSAcceptance() {}
  PyMCSAcceptance(PyObject *obj) : PyMCSWrapper(obj) { extractPyMCSWrapper(); }
  ~PyMCSAcceptance() override {}
  PyMCSAcceptance *extractPyObject() const {
    auto res = dynamic_cast<PyMCSAcceptance *>(pyObjectExtract());
    if (!res) {
      failedToExtractPyObject();
    }
    return res;
  }
  virtual bool operator()() const {
    errorNotOverridden();
    return false;
  }
};

class PyMCSParameters : public boost::noncopyable {
 public:
  PyMCSParameters() : p(new MCSParameters()) {
    cfud.mcsParameters = p.get();
    pcud.mcsProgressData = nullptr;
  }
  PyMCSParameters(const MCSParameters &other) : PyMCSParameters() {
    *p = other;
  }
  PyMCSParameters(const MCSParameters &other,
                  const PyProgressCallbackUserData &pcudOther)
      : PyMCSParameters(other) {
    pcud.pyMCSProgress = pcudOther.pyMCSProgress;
    cfud.pyAtomBondCompData = pcudOther.pyAtomBondCompData;
  }
  PyMCSParameters(const MCSParameters &other,
                  const PyMCSFinalMatchCheckFunctionUserData &fmudOther)
      : PyMCSParameters(other) {
    fmud.pyMCSFinalMatchCheck = fmudOther.pyMCSFinalMatchCheck;
    cfud.pyAtomBondCompData = fmudOther.pyAtomBondCompData;
  }
  PyMCSParameters(const MCSParameters &other,
                  const PyMCSAcceptanceFunctionUserData &afudOther)
      : PyMCSParameters(other) {
    afud.pyMCSAcceptance = afudOther.pyMCSAcceptance;
    cfud.pyAtomBondCompData = afudOther.pyAtomBondCompData;
  }
  const MCSParameters *get() const { return p.get(); }
  bool getMaximizeBonds() const { return p->MaximizeBonds; }
  void setMaximizeBonds(bool value) { p->MaximizeBonds = value; }
  double getThreshold() const { return p->Threshold; }
  void setThreshold(double value) { p->Threshold = value; }
  unsigned int getTimeout() const { return p->Timeout; }
  void setTimeout(unsigned int value) { p->Timeout = value; }
  bool getVerbose() const { return p->Verbose; }
  void setVerbose(bool value) { p->Verbose = value; }
  bool getUseCliqueDection() const { return p->useCliqueDetection; }
  void setUseCliqueDetection(bool value) { p->useCliqueDetection = value; }
  void setMinMCSSSize(unsigned int value) { p->MinMCSSSize = value; }
  unsigned int getMinMCSSSize() { return p->MinMCSSSize; }
  const MCSAtomCompareParameters &getAtomCompareParameters() const {
    return p->AtomCompareParameters;
  }
  void setAtomCompareParameters(const MCSAtomCompareParameters &value) {
    p->AtomCompareParameters = value;
  }
  const MCSBondCompareParameters &getBondCompareParameters() const {
    return p->BondCompareParameters;
  }
  void setBondCompareParameters(const MCSBondCompareParameters &value) {
    p->BondCompareParameters = value;
  }
  std::string getInitialSeed() const { return p->InitialSeed; }
  void setInitialSeed(const std::string &value) { p->InitialSeed = value; }
  bool getStoreAll() const { return p->StoreAll; }
  void setStoreAll(bool value) { p->StoreAll = value; }
  void setMCSAtomTyper(PyObject *atomComp) {
    PyMCSAtomCompare pyMCSAtomCompare(atomComp);
    AtomComparator ac;
    if (pyMCSAtomCompare.extractAtomComparator(ac)) {
      p->setMCSAtomTyperFromEnum(ac);
    } else {
      p->CompareFunctionsUserData = &cfud;
      p->AtomTyper = MCSAtomComparePyFunc;
      cfud.pyAtomBondCompData.pyAtomComp = pyMCSAtomCompare.pyObject();
      cfud.mcsParameters = p.get();
    }
  }
  python::object getMCSAtomTyper() const {
    static const std::map<RDKit::MCSAtomCompareFunction, RDKit::AtomComparator>
        atomTyperToComp = {
            {MCSAtomCompareAny, AtomCompareAny},
            {MCSAtomCompareElements, AtomCompareElements},
            {MCSAtomCompareIsotopes, AtomCompareIsotopes},
            {MCSAtomCompareAnyHeavyAtom, AtomCompareAnyHeavyAtom}};
    if (!cfud.pyAtomBondCompData.pyAtomComp.is_none()) {
      return cfud.pyAtomBondCompData.pyAtomComp;
    }
    python::object res;
    try {
      res = python::object(atomTyperToComp.at(p->AtomTyper));
    } catch (const std::out_of_range &) {
      PyErr_SetString(PyExc_TypeError, "Unknown AtomTyper");
      python::throw_error_already_set();
    }
    return res;
  }
  void setMCSBondTyper(PyObject *bondComp) {
    PyMCSBondCompare pyMCSBondCompare(bondComp);
    BondComparator bc;
    if (pyMCSBondCompare.extractBondComparator(bc)) {
      p->setMCSBondTyperFromEnum(bc);
    } else {
      p->CompareFunctionsUserData = &cfud;
      p->BondTyper = MCSBondComparePyFunc;
      cfud.pyAtomBondCompData.pyBondComp = pyMCSBondCompare.pyObject();
      auto bc = pyMCSBondCompare.extractPyObject();
      bc->mcsParameters = p.get();
      cfud.mcsParameters = p.get();
    }
  }
  python::object getMCSBondTyper() const {
    static const std::map<RDKit::MCSBondCompareFunction, RDKit::BondComparator>
        bondTyperToComp = {{MCSBondCompareAny, BondCompareAny},
                           {MCSBondCompareOrder, BondCompareOrder},
                           {MCSBondCompareOrderExact, BondCompareOrderExact}};
    if (!cfud.pyAtomBondCompData.pyBondComp.is_none()) {
      return cfud.pyAtomBondCompData.pyBondComp;
    }
    python::object res;
    try {
      res = python::object(bondTyperToComp.at(p->BondTyper));
    } catch (const std::out_of_range &) {
      PyErr_SetString(PyExc_TypeError, "Unknown BondTyper");
      python::throw_error_already_set();
    }
    return res;
  }
  void setMCSProgressCallback(PyObject *progress) {
    PyMCSProgress pyMCSProgress(progress);
    p->ProgressCallbackUserData = &pcud;
    p->ProgressCallback = MCSProgressCallbackPyFunc;
    pcud.pyMCSProgress = pyMCSProgress.pyObject();
    pcud.pyAtomBondCompData = cfud.pyAtomBondCompData;
  }
  python::object getMCSProgressCallback() {
    if (!pcud.pyMCSProgress.is_none()) {
      return pcud.pyMCSProgress;
    }
    return python::object();
  }
  void setFinalMatchCheck(PyObject *finalMatchCheck) {
    PyMCSFinalMatchCheck pyMCSFinalMatchCheck(finalMatchCheck);
    p->FinalMatchChecker = MCSFinalMatchCheckPyFunc;
    p->FinalMatchCheckerUserData = &fmud;
    fmud.pyMCSFinalMatchCheck = pyMCSFinalMatchCheck.pyObject();
    fmud.pyAtomBondCompData = cfud.pyAtomBondCompData;
  }
  python::object getFinalMatchCheck() const {
    if (!fmud.pyMCSFinalMatchCheck.is_none()) {
      return fmud.pyMCSFinalMatchCheck;
    }
    return python::object();
  }
  void setShouldAcceptMCS(PyObject *mcsAcceptance) {
    PyMCSAcceptance pyMCSAcceptance(mcsAcceptance);
    p->ShouldAcceptMCS = MCSAcceptancePyFunc;
    p->ShouldAcceptMCSUserData = &afud;
    afud.pyMCSAcceptance = pyMCSAcceptance.pyObject();
    afud.pyAtomBondCompData = cfud.pyAtomBondCompData;
  }
  python::object getShouldAcceptMCS() const {
    if (!afud.pyMCSAcceptance.is_none()) {
      return afud.pyMCSAcceptance;
    }
    return python::object();
  }

 private:
  static bool MCSAtomComparePyFunc(const MCSAtomCompareParameters &p,
                                   const ROMol &mol1, unsigned int atom1,
                                   const ROMol &mol2, unsigned int atom2,
                                   void *userData) {
    PRECONDITION(userData, "userData must not be NULL");
    PyCompareFunctionUserData *cfud =
        static_cast<PyCompareFunctionUserData *>(userData);
    CHECK_INVARIANT(cfud, "");
    bool res = false;
    {
      PyGILStateHolder h;
      res = python::call_method<bool>(
          cfud->pyAtomBondCompData.pyAtomComp.ptr(), COMPARE_FUNC_NAME,
          boost::ref(p), boost::ref(mol1), atom1, boost::ref(mol2), atom2);
    }
    return res;
  }
  static bool MCSBondComparePyFunc(const MCSBondCompareParameters &p,
                                   const ROMol &mol1, unsigned int bond1,
                                   const ROMol &mol2, unsigned int bond2,
                                   void *userData) {
    PRECONDITION(userData, "userData must not be NULL");
    PyCompareFunctionUserData *cfud =
        static_cast<PyCompareFunctionUserData *>(userData);
    CHECK_INVARIANT(cfud, "");
    bool res = false;
    {
      PyGILStateHolder h;
      res = python::call_method<bool>(
          cfud->pyAtomBondCompData.pyBondComp.ptr(), COMPARE_FUNC_NAME,
          boost::ref(p), boost::ref(mol1), bond1, boost::ref(mol2), bond2);
    }
    return res;
  }
  static bool MCSProgressCallbackPyFunc(const MCSProgressData &stat,
                                        const MCSParameters &params,
                                        void *userData) {
    PRECONDITION(userData, "userData must not be NULL");
    PyProgressCallbackUserData *pcud =
        static_cast<PyProgressCallbackUserData *>(userData);
    bool res = false;
    CHECK_INVARIANT(pcud, "");
    MCSParameters paramsCopy(params);
    if (pcud->pyAtomBondCompData.standardBondTyperFunc) {
      paramsCopy.BondTyper = pcud->pyAtomBondCompData.standardBondTyperFunc;
    }
    {
      PyGILStateHolder h;
      PyMCSParameters ps(paramsCopy, *pcud);
      PyMCSProgressData pd(stat);
      res = python::call_method<bool>(pcud->pyMCSProgress.ptr(),
                                      CALLBACK_FUNC_NAME, boost::ref(pd),
                                      boost::ref(ps));
    }
    return res;
  }
  static bool MCSFinalMatchCheckPyFunc(
      const std::uint32_t c1[], const std::uint32_t c2[], const ROMol &mol1,
      const FMCS::Graph &query, const ROMol &mol2, const FMCS::Graph &target,
      const MCSParameters *params) {
    PRECONDITION(params, "params must not be NULL");
    PyMCSFinalMatchCheckFunctionUserData *fmud =
        static_cast<PyMCSFinalMatchCheckFunctionUserData *>(
            params->FinalMatchCheckerUserData);
    CHECK_INVARIANT(fmud, "");
    bool res = false;
    {
      PyGILStateHolder h;
      PyMCSParameters ps(*params, *fmud);
      auto numMcsAtoms = boost::num_vertices(query);
      PyObject *pyAtomIdxMatch = PyTuple_New(numMcsAtoms);
      for (unsigned int i = 0; i < numMcsAtoms; ++i) {
        PyObject *pyPair = PyTuple_New(2);
        PyTuple_SetItem(pyPair, 0,
                        PyInt_FromLong(query[c1[boost::vertex(i, query)]]));
        PyTuple_SetItem(pyPair, 1,
                        PyInt_FromLong(target[c2[boost::vertex(i, query)]]));
        PyTuple_SetItem(pyAtomIdxMatch, i, pyPair);
      }
      auto numMcsBonds = boost::num_edges(query);
      auto queryBondIt = boost::edges(query).first;
      PyObject *pyBondIdxMatch = PyTuple_New(numMcsBonds);
      for (unsigned int i = 0; i < numMcsBonds; ++i, ++queryBondIt) {
        const auto queryBond = mol1.getBondBetweenAtoms(
            query[c1[boost::source(*queryBondIt, query)]],
            query[c1[boost::target(*queryBondIt, query)]]);
        CHECK_INVARIANT(queryBond, "");
        const auto targetBond = mol2.getBondBetweenAtoms(
            target[c2[boost::source(*queryBondIt, query)]],
            target[c2[boost::target(*queryBondIt, query)]]);
        CHECK_INVARIANT(targetBond, "");
        PyObject *pyPair = PyTuple_New(2);
        PyTuple_SetItem(pyPair, 0, PyInt_FromLong(queryBond->getIdx()));
        PyTuple_SetItem(pyPair, 1, PyInt_FromLong(targetBond->getIdx()));
        PyTuple_SetItem(pyBondIdxMatch, i, pyPair);
      }
      res = python::call_method<bool>(
          fmud->pyMCSFinalMatchCheck.ptr(), CALLBACK_FUNC_NAME,
          boost::ref(mol1), boost::ref(mol2), python::handle<>(pyAtomIdxMatch),
          python::handle<>(pyBondIdxMatch), boost::ref(ps));
    }
    return res;
  }
  static bool MCSAcceptancePyFunc(
      const ROMol &query, const ROMol &target,
      const std::vector<std::pair<int, int>> &atomIdxMatch,
      const std::vector<std::pair<int, int>> &bondIdxMatch,
      const MCSParameters *params) {
    PRECONDITION(params, "params must not be NULL");
    PyMCSAcceptanceFunctionUserData *afud =
        static_cast<PyMCSAcceptanceFunctionUserData *>(
            params->ShouldAcceptMCSUserData);
    CHECK_INVARIANT(afud, "");
    bool res = false;
    {
      PyGILStateHolder h;
      PyMCSParameters ps(*params, *afud);
      auto pyAtomIdxMatch = convertMatchesToTupleOfPairs(atomIdxMatch);
      auto pyBondIdxMatch = convertMatchesToTupleOfPairs(bondIdxMatch);
      res = python::call_method<bool>(
          afud->pyMCSAcceptance.ptr(), CALLBACK_FUNC_NAME, boost::ref(query),
          boost::ref(target), python::handle<>(pyAtomIdxMatch),
          python::handle<>(pyBondIdxMatch), boost::ref(ps));
    }
    return res;
  }
  std::unique_ptr<MCSParameters> p;
  PyCompareFunctionUserData cfud;
  PyProgressCallbackUserData pcud;
  PyMCSFinalMatchCheckFunctionUserData fmud;
  PyMCSAcceptanceFunctionUserData afud;
};

MCSResult *FindMCSWrapper(python::object mols, bool maximizeBonds,
                          double threshold, unsigned int timeout, bool verbose,
                          bool useCliqueDetection, unsigned int minMCSSSize,
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
  p.useCliqueDetection = useCliqueDetection;
  p.MinMCSSSize = minMCSSSize;
  p.InitialSeed = seedSmarts;
  p.AtomCompareParameters.MatchValences = matchValences;
  p.AtomCompareParameters.MatchChiralTag = matchChiralTag;
  p.AtomCompareParameters.RingMatchesRingOnly = ringMatchesRingOnly;
  p.setMCSAtomTyperFromEnum(atomComp);
  p.setMCSBondTyperFromEnum(bondComp);
  p.BondCompareParameters.RingMatchesRingOnly = ringMatchesRingOnly;
  p.BondCompareParameters.CompleteRingsOnly = completeRingsOnly;
  p.BondCompareParameters.MatchFusedRings = (ringComp != IgnoreRingFusion);
  p.BondCompareParameters.MatchFusedRingsStrict =
      (ringComp == StrictRingFusion);

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
  {
    NOGIL gil;
    res = new MCSResult(findMCS(ms, pyMcsParams.get()));
  }
  return res;
}
}  // namespace RDKit
namespace {
python::object degenerateSmartsQueryMolDictHelper(
    const RDKit::MCSResult &self) {
  python::dict res;
  for (const auto &pair : self.DegenerateSmartsQueryMolDict) {
    res[pair.first] = pair.second;
  }
  return res;
}
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
                      "if True, the MCS calculation did not finish")
        .add_property(
            "degenerateSmartsQueryMolDict", degenerateSmartsQueryMolDictHelper,
            "Dictionary collecting all degenerate (SMARTS, queryMol) pairs "
            "(empty if MCSParameters.StoreAll is False)");
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
       python::arg("verbose") = false,
       python::arg("useCliqueDetection") = false,
       python::arg("minMCSSSize") = 1, python::arg("matchValences") = false,
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
      .add_property("MaximizeBonds", &RDKit::PyMCSParameters::getMaximizeBonds,
                    &RDKit::PyMCSParameters::setMaximizeBonds,
                    "toggles maximizing the number of bonds (instead of the "
                    "number of atoms)")
      .add_property("Threshold", &RDKit::PyMCSParameters::getThreshold,
                    &RDKit::PyMCSParameters::setThreshold,
                    "fraction of the dataset that must contain the MCS")
      .add_property("Timeout", &RDKit::PyMCSParameters::getTimeout,
                    &RDKit::PyMCSParameters::setTimeout,
                    "timeout (in seconds) for the calculation")
      .add_property("Verbose", &RDKit::PyMCSParameters::getVerbose,
                    &RDKit::PyMCSParameters::setVerbose, "toggles verbose mode")
      .add_property("UseCliqueDetection",
                    &RDKit::PyMCSParameters::getUseCliqueDection,
                    &RDKit::PyMCSParameters::setUseCliqueDetection,
                    "toggles useCliqueDetection mode")
      .add_property("AtomCompareParameters",
                    python::make_function(
                        &RDKit::PyMCSParameters::getAtomCompareParameters,
                        python::return_internal_reference<>()),
                    &RDKit::PyMCSParameters::setAtomCompareParameters,
                    "parameters for comparing atoms")
      .add_property("BondCompareParameters",
                    python::make_function(
                        &RDKit::PyMCSParameters::getBondCompareParameters,
                        python::return_internal_reference<>()),
                    &RDKit::PyMCSParameters::setBondCompareParameters,
                    "parameters for comparing bonds")
      .add_property("AtomTyper", &RDKit::PyMCSParameters::getMCSAtomTyper,
                    &RDKit::PyMCSParameters::setMCSAtomTyper,
                    "atom typer to be used. Must be one of the "
                    "members of the rdFMCS.AtomCompare class or "
                    "an instance of a user-defined subclass of "
                    "rdFMCS.MCSAtomCompare")
      .add_property("BondTyper", &RDKit::PyMCSParameters::getMCSBondTyper,
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
      .add_property("FinalMatchChecker",
                    &RDKit::PyMCSParameters::getFinalMatchCheck,
                    &RDKit::PyMCSParameters::setFinalMatchCheck,
                    "seed final match checker callback class. Must be a "
                    "user-defined subclass of rdFMCS.MCSFinalMatchCheck")
      .add_property("ShouldAcceptMCS",
                    &RDKit::PyMCSParameters::getShouldAcceptMCS,
                    &RDKit::PyMCSParameters::setShouldAcceptMCS,
                    "MCS acceptance callback class. Must be a "
                    "user-defined subclass of rdFMCS.MCSAcceptance")
      .add_property("InitialSeed", &RDKit::PyMCSParameters::getInitialSeed,
                    &RDKit::PyMCSParameters::setInitialSeed,
                    "SMILES string to be used as the seed of the MCS")
      .add_property("StoreAll", &RDKit::PyMCSParameters::getStoreAll,
                    &RDKit::PyMCSParameters::setStoreAll,
                    "toggles storage of degenerate MCSs")
      .add_property("MinMCSSSize", &RDKit::PyMCSParameters::getMinMCSSSize,
                    &RDKit::PyMCSParameters::setMinMCSSSize,
                    "Minimum number of atoms in MCSS.");

  python::class_<RDKit::MCSAtomCompareParameters, boost::noncopyable>(
      "MCSAtomCompareParameters",
      "Parameters controlling how atom-atom matching is done")
      .def_readwrite("MatchValences",
                     &RDKit::MCSAtomCompareParameters::MatchValences,
                     "include atom valences in the match")
      .def_readwrite("MatchChiralTag",
                     &RDKit::MCSAtomCompareParameters::MatchChiralTag,
                     "include atom chirality in the match")
      .def_readwrite("MaxDistance",
                     &RDKit::MCSAtomCompareParameters::MaxDistance,
                     "Require atoms to be within this many angstroms in 3D")
      .def_readwrite("MatchFormalCharge",
                     &RDKit::MCSAtomCompareParameters::MatchFormalCharge,
                     "include formal charge in the match")
      .def_readwrite("RingMatchesRingOnly",
                     &RDKit::MCSAtomCompareParameters::RingMatchesRingOnly,
                     "ring atoms are only allowed to match other ring atoms")
      .def_readwrite("CompleteRingsOnly",
                     &RDKit::MCSAtomCompareParameters::CompleteRingsOnly,
                     "results cannot include lone ring atoms")
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
      .def_readwrite(
          "MatchFusedRings", &RDKit::MCSBondCompareParameters::MatchFusedRings,
          "enforce check on ring fusion, i.e. alpha-methylnaphthalene "
          "won't match beta-methylnaphtalene, but decalin "
          "will match cyclodecane unless MatchFusedRingsStrict is True")
      .def_readwrite(
          "MatchFusedRingsStrict",
          &RDKit::MCSBondCompareParameters::MatchFusedRingsStrict,
          "only enforced if MatchFusedRings is True; the ring fusion "
          "must be the same in both query and target, i.e. decalin "
          "won't match cyclodecane")
      .def_readwrite("MatchStereo",
                     &RDKit::MCSBondCompareParameters::MatchStereo,
                     "include bond stereo in the comparison");

  python::class_<RDKit::PyMCSProgress, boost::noncopyable>(
      "MCSProgress",
      "Base class. Subclass and override "
      "MCSProgress." CALLBACK_FUNC_NAME
      "() "
      "to define a custom callback function")
      .def(
          CALLBACK_FUNC_NAME, &RDKit::PyMCSProgress::operator(),
          (python::arg("self"), python::arg("stat"), python::arg("parameters")),
          "override to implement a custom progress callback");

  python::class_<RDKit::PyMCSProgressData, boost::noncopyable>(
      "MCSProgressData", "Information about the MCS progress")
      .add_property("numAtoms", &RDKit::PyMCSProgressData::getNumAtoms,
                    "number of atoms in MCS")
      .add_property("numBonds", &RDKit::PyMCSProgressData::getNumBonds,
                    "number of bonds in MCS")
      .add_property("seedProcessed",
                    &RDKit::PyMCSProgressData::getSeedProcessed,
                    "number of processed seeds");

  python::class_<RDKit::PyMCSAtomCompare, boost::noncopyable>(
      "MCSAtomCompare",
      "Base class. Subclass and override "
      "MCSAtomCompare." COMPARE_FUNC_NAME
      "() to define custom "
      "atom compare functions, then set MCSParameters.AtomTyper "
      "to an instance of the subclass")
      .def("CheckAtomRingMatch", &RDKit::PyMCSAtomCompare::checkAtomRingMatch,
           (python::arg("self"), python::arg("parameters"), python::arg("mol1"),
            python::arg("atom1"), python::arg("mol2"), python::arg("atom2")),
           "Return True if both atoms are, or are not, in a ring")
      .def("CheckAtomCharge", &RDKit::PyMCSAtomCompare::checkAtomCharge,
           (python::arg("self"), python::arg("parameters"), python::arg("mol1"),
            python::arg("atom1"), python::arg("mol2"), python::arg("atom2")),
           "Return True if both atoms have the same formal charge")
      .def("CheckAtomChirality", &RDKit::PyMCSAtomCompare::checkAtomChirality,
           (python::arg("self"), python::arg("parameters"), python::arg("mol1"),
            python::arg("atom1"), python::arg("mol2"), python::arg("atom2")),
           "Return True if both atoms have, or have not, a chiral tag")
      .def(COMPARE_FUNC_NAME, &RDKit::PyMCSAtomCompare::operator(),
           (python::arg("self"), python::arg("parameters"), python::arg("mol1"),
            python::arg("atom1"), python::arg("mol2"), python::arg("atom2")),
           "override to implement custom atom comparison");

  python::class_<RDKit::PyMCSBondCompare, boost::noncopyable>(
      "MCSBondCompare",
      "Base class. Subclass and override "
      "MCSBondCompare." COMPARE_FUNC_NAME
      "() to define custom "
      "bond compare functions, then set MCSParameters.BondTyper "
      "to an instance of the subclass")
      .def("CheckBondStereo", &RDKit::PyMCSBondCompare::checkBondStereo,
           (python::arg("self"), python::arg("parameters"), python::arg("mol1"),
            python::arg("bond1"), python::arg("mol2"), python::arg("bond2")),
           "Return True if both bonds have, or have not, a stereo descriptor")
      .def("CheckBondRingMatch", &RDKit::PyMCSBondCompare::checkBondRingMatch,
           (python::arg("self"), python::arg("parameters"), python::arg("mol1"),
            python::arg("bond1"), python::arg("mol2"), python::arg("bond2")),
           "Return True if both bonds are, or are not, part of a ring")
      .def(COMPARE_FUNC_NAME, &RDKit::PyMCSBondCompare::operator(),
           (python::arg("self"), python::arg("parameters"), python::arg("mol1"),
            python::arg("bond1"), python::arg("mol2"), python::arg("bond2")),
           "override to implement custom bond comparison");

  python::class_<RDKit::PyMCSFinalMatchCheck, boost::noncopyable>(
      "MCSFinalMatchCheck",
      "Base class. Subclass and override "
      "MCSFinalMatchCheck." CALLBACK_FUNC_NAME
      "() "
      "to define a custom boolean callback function. "
      "Returning True will cause the growing seed to be accepted, "
      "False to be rejected")
      .def(CALLBACK_FUNC_NAME, &RDKit::PyMCSFinalMatchCheck::operator(),
           (python::arg("self")),
           "override to implement a custom seed final match checker callback");

  python::class_<RDKit::PyMCSAcceptance, boost::noncopyable>(
      "MCSAcceptance",
      "Base class. Subclass and override "
      "MCSAcceptance." CALLBACK_FUNC_NAME
      "() "
      "to define a custom boolean callback function. "
      "Returning True will cause the MCS candidate to be accepted, "
      "False to be rejected")
      .def(CALLBACK_FUNC_NAME, &RDKit::PyMCSAcceptance::operator(),
           (python::arg("self")),
           "override to implement a custom MCS acceptance callback");

  python::def("FindMCS", RDKit::FindMCSWrapper2,
              (python::arg("mols"), python::arg("parameters")),
              python::return_value_policy<python::manage_new_object>(),
              docString.c_str());
}
