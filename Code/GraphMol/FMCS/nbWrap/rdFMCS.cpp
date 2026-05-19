//
//  Copyright (C) 2014-2026 Novartis Institutes for BioMedical Research
//  and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <GraphMol/ROMol.h>
#include <RDBoost/Wrap_nb.h>
#include <GraphMol/FMCS/FMCS.h>

#define COMPARE_FUNC_NAME "__call__"
#define CALLBACK_FUNC_NAME "__call__"

namespace nb = nanobind;
using namespace nb::literals;

namespace {
std::shared_ptr<RDKit::ROMol> toStd(const RDKit::ROMOL_SPTR &bptr) {
  return {bptr.get(), [b = bptr](RDKit::ROMol *) {}};
}
}  // namespace

namespace RDKit {

// Helper: build atom index match tuple from FinalMatchCheck c-arrays
static nb::tuple buildAtomIdxMatchTuple(const std::uint32_t c1[],
                                        const std::uint32_t c2[],
                                        const FMCS::Graph &query,
                                        const FMCS::Graph &target) {
  auto numMcsAtoms = boost::num_vertices(query);
  nb::list atomList;
  for (unsigned int i = 0; i < numMcsAtoms; ++i) {
    atomList.append(nb::make_tuple(
        (long)query[c1[boost::vertex(i, query)]],
        (long)target[c2[boost::vertex(i, query)]]));
  }
  return nb::tuple(atomList);
}

// Helper: build bond index match tuple from FinalMatchCheck c-arrays
static nb::tuple buildBondIdxMatchTuple(const std::uint32_t c1[],
                                        const std::uint32_t c2[],
                                        const ROMol &mol1,
                                        const FMCS::Graph &query,
                                        const ROMol &mol2,
                                        const FMCS::Graph &target) {
  auto numMcsBonds = boost::num_edges(query);
  auto queryBondIt = boost::edges(query).first;
  nb::list bondList;
  for (unsigned int i = 0; i < numMcsBonds; ++i, ++queryBondIt) {
    const auto queryBond = mol1.getBondBetweenAtoms(
        query[c1[boost::source(*queryBondIt, query)]],
        query[c1[boost::target(*queryBondIt, query)]]);
    CHECK_INVARIANT(queryBond, "");
    const auto targetBond = mol2.getBondBetweenAtoms(
        target[c2[boost::source(*queryBondIt, query)]],
        target[c2[boost::target(*queryBondIt, query)]]);
    CHECK_INVARIANT(targetBond, "");
    bondList.append(nb::make_tuple((long)queryBond->getIdx(),
                                   (long)targetBond->getIdx()));
  }
  return nb::tuple(bondList);
}

// Helper: convert a vector of (int,int) pairs to a Python tuple of 2-tuples
static nb::tuple convertMatchesToTupleOfPairs(
    const std::vector<std::pair<int, int>> &matches) {
  nb::list res;
  for (const auto &pair : matches) {
    res.append(nb::make_tuple(pair.first, pair.second));
  }
  return nb::tuple(res);
}

// -------------------------------------------------------------------
// PyMCSAtomCompare: wraps a Python object that is either an AtomComparator
// enum value or a user-defined subclass implementing __call__().
// -------------------------------------------------------------------
struct PyMCSAtomCompare {
  PyMCSAtomCompare() = default;
  explicit PyMCSAtomCompare(nb::object obj) : d_pyObject(std::move(obj)) {}

  // Returns true if the object is a predefined AtomComparator enum value.
  bool extractAtomComparator(AtomComparator &ac) const {
    try {
      ac = nb::cast<AtomComparator>(d_pyObject);
      return true;
    } catch (const nb::cast_error &) {
      return false;
    }
  }

  const nb::object &pyObject() const { return d_pyObject; }

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

 private:
  nb::object d_pyObject;
};

// -------------------------------------------------------------------
// PyMCSBondCompare: same pattern as PyMCSAtomCompare but for bonds.
// -------------------------------------------------------------------
struct PyMCSBondCompare {
  PyMCSBondCompare() = default;
  explicit PyMCSBondCompare(nb::object obj) : d_pyObject(std::move(obj)) {}

  bool extractBondComparator(BondComparator &bc) const {
    try {
      bc = nb::cast<BondComparator>(d_pyObject);
      return true;
    } catch (const nb::cast_error &) {
      return false;
    }
  }

  const nb::object &pyObject() const { return d_pyObject; }

  inline bool checkBondStereo(const MCSBondCompareParameters &p,
                              const ROMol &mol1, unsigned int bond1,
                              const ROMol &mol2, unsigned int bond2) const {
    return RDKit::checkBondStereo(p, mol1, bond1, mol2, bond2);
  }
  inline bool checkBondRingMatch(const MCSBondCompareParameters &p,
                                 const ROMol &mol1, unsigned int bond1,
                                 const ROMol &mol2, unsigned int bond2) const {
    return RDKit::checkBondRingMatch(p, mol1, bond1, mol2, bond2);
  }

 private:
  nb::object d_pyObject;
};

// -------------------------------------------------------------------
// User-data structs stored as pointers in MCSParameters callbacks.
// These outlive each call since they are owned by PyMCSParameters.
// -------------------------------------------------------------------
struct PyAtomBondCompData {
  nb::object pyAtomComp = nb::none();
  nb::object pyBondComp = nb::none();
  MCSBondCompareFunction standardBondTyperFunc = nullptr;
};

struct PyBaseUserData {
  PyAtomBondCompData pyAtomBondCompData;
};

struct PyCompareFunctionUserData : public PyBaseUserData {
  const MCSParameters *mcsParameters = nullptr;
};

struct PyProgressCallbackUserData : public PyBaseUserData {
  const MCSProgressData *mcsProgressData = nullptr;
  nb::object pyMCSProgress = nb::none();
};

struct PyMCSFinalMatchCheckFunctionUserData : public PyBaseUserData {
  nb::object pyMCSFinalMatchCheck = nb::none();
};

struct PyMCSAcceptanceFunctionUserData : public PyBaseUserData {
  nb::object pyMCSAcceptance = nb::none();
};

// -------------------------------------------------------------------
// PyMCSProgressData: read-only view of MCSProgressData for Python callbacks.
// -------------------------------------------------------------------
class PyMCSProgressData {
 public:
  PyMCSProgressData()
      : pd(new MCSProgressData()), pcud(new PyProgressCallbackUserData()) {
    pcud->mcsProgressData = pd.get();
  }
  explicit PyMCSProgressData(const MCSProgressData &other)
      : PyMCSProgressData() {
    *pd = other;
  }
  unsigned int getNumAtoms() const { return pd->NumAtoms; }
  unsigned int getNumBonds() const { return pd->NumBonds; }
  unsigned int getSeedProcessed() const { return pd->SeedProcessed; }

 private:
  std::unique_ptr<MCSProgressData> pd;
  std::unique_ptr<PyProgressCallbackUserData> pcud;
};

// -------------------------------------------------------------------
// PyMCSParameters: Python-facing wrapper around MCSParameters.
// Manages all Python callback objects and wires them into C++ callbacks.
// -------------------------------------------------------------------
class PyMCSParameters {
 public:
  PyMCSParameters() : p(new MCSParameters()) {
    cfud.mcsParameters = p.get();
    pcud.mcsProgressData = nullptr;
  }
  explicit PyMCSParameters(const MCSParameters &other) : PyMCSParameters() {
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

  PyMCSParameters(const PyMCSParameters &) = delete;
  PyMCSParameters &operator=(const PyMCSParameters &) = delete;

  const MCSParameters *get() const { return p.get(); }

  bool getMaximizeBonds() const { return p->MaximizeBonds; }
  void setMaximizeBonds(bool value) { p->MaximizeBonds = value; }
  double getThreshold() const { return p->Threshold; }
  void setThreshold(double value) { p->Threshold = value; }
  unsigned int getTimeout() const { return p->Timeout; }
  void setTimeout(unsigned int value) { p->Timeout = value; }
  bool getVerbose() const { return p->Verbose; }
  void setVerbose(bool value) { p->Verbose = value; }
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

  void setMCSAtomTyper(nb::object atomComp) {
    PyMCSAtomCompare pyMCSAtomCompare(atomComp);
    AtomComparator ac;
    if (pyMCSAtomCompare.extractAtomComparator(ac)) {
      p->setMCSAtomTyperFromEnum(ac);
    } else {
      // Require that the Python class itself (not an inherited base) defines a
      // callable __call__ method.
      PyObject *callInDict =
          PyDict_GetItemString(Py_TYPE(atomComp.ptr())->tp_dict, COMPARE_FUNC_NAME);
      if (!callInDict || !PyCallable_Check(callInDict)) {
        throw nb::type_error(
            "AtomTyper must be an AtomCompare enum value or an instance of a "
            "rdFMCS.MCSAtomCompare subclass with a callable __call__() method");
      }
      p->CompareFunctionsUserData = &cfud;
      p->AtomTyper = MCSAtomComparePyFunc;
      cfud.pyAtomBondCompData.pyAtomComp = pyMCSAtomCompare.pyObject();
      cfud.mcsParameters = p.get();
    }
  }
  nb::object getMCSAtomTyper() const {
    static const std::map<RDKit::MCSAtomCompareFunction, RDKit::AtomComparator>
        atomTyperToComp = {
            {MCSAtomCompareAny, AtomCompareAny},
            {MCSAtomCompareElements, AtomCompareElements},
            {MCSAtomCompareIsotopes, AtomCompareIsotopes},
            {MCSAtomCompareAnyHeavyAtom, AtomCompareAnyHeavyAtom}};
    if (!cfud.pyAtomBondCompData.pyAtomComp.is_none()) {
      return cfud.pyAtomBondCompData.pyAtomComp;
    }
    try {
      return nb::cast(atomTyperToComp.at(p->AtomTyper));
    } catch (const std::out_of_range &) {
      PyErr_SetString(PyExc_TypeError, "Unknown AtomTyper");
      throw nb::python_error();
    }
  }
  void setMCSBondTyper(nb::object bondComp) {
    PyMCSBondCompare pyMCSBondCompare(bondComp);
    BondComparator bc;
    if (pyMCSBondCompare.extractBondComparator(bc)) {
      p->setMCSBondTyperFromEnum(bc);
    } else {
      // Require that the Python class itself (not an inherited base) defines a
      // callable __call__ method.
      PyObject *callInDict =
          PyDict_GetItemString(Py_TYPE(bondComp.ptr())->tp_dict, COMPARE_FUNC_NAME);
      if (!callInDict || !PyCallable_Check(callInDict)) {
        throw nb::type_error(
            "BondTyper must be a BondCompare enum value or an instance of a "
            "rdFMCS.MCSBondCompare subclass with a callable __call__() method");
      }
      p->CompareFunctionsUserData = &cfud;
      p->BondTyper = MCSBondComparePyFunc;
      cfud.pyAtomBondCompData.pyBondComp = pyMCSBondCompare.pyObject();
      cfud.mcsParameters = p.get();
    }
  }
  nb::object getMCSBondTyper() const {
    static const std::map<RDKit::MCSBondCompareFunction, RDKit::BondComparator>
        bondTyperToComp = {{MCSBondCompareAny, BondCompareAny},
                           {MCSBondCompareOrder, BondCompareOrder},
                           {MCSBondCompareOrderExact, BondCompareOrderExact}};
    if (!cfud.pyAtomBondCompData.pyBondComp.is_none()) {
      return cfud.pyAtomBondCompData.pyBondComp;
    }
    try {
      return nb::cast(bondTyperToComp.at(p->BondTyper));
    } catch (const std::out_of_range &) {
      PyErr_SetString(PyExc_TypeError, "Unknown BondTyper");
      throw nb::python_error();
    }
  }
  void setMCSProgressCallback(nb::object progress) {
    PyObject *callInDict = PyDict_GetItemString(
        Py_TYPE(progress.ptr())->tp_dict, CALLBACK_FUNC_NAME);
    if (!callInDict || !PyCallable_Check(callInDict)) {
      throw nb::type_error("The __call__() method must be overridden in "
                           "the rdFMCS.MCSProgress subclass");
    }
    p->ProgressCallbackUserData = &pcud;
    p->ProgressCallback = MCSProgressCallbackPyFunc;
    pcud.pyMCSProgress = std::move(progress);
    pcud.pyAtomBondCompData = cfud.pyAtomBondCompData;
  }
  nb::object getMCSProgressCallback() const {
    if (!pcud.pyMCSProgress.is_none()) {
      return pcud.pyMCSProgress;
    }
    return nb::none();
  }
  void setFinalMatchCheck(nb::object finalMatchCheck) {
    PyObject *callInDict = PyDict_GetItemString(
        Py_TYPE(finalMatchCheck.ptr())->tp_dict, CALLBACK_FUNC_NAME);
    if (!callInDict || !PyCallable_Check(callInDict)) {
      throw nb::type_error("The __call__() method must be overridden in "
                           "the rdFMCS.MCSFinalMatchCheck subclass");
    }
    p->FinalMatchChecker = MCSFinalMatchCheckPyFunc;
    p->FinalMatchCheckerUserData = &fmud;
    fmud.pyMCSFinalMatchCheck = std::move(finalMatchCheck);
    fmud.pyAtomBondCompData = cfud.pyAtomBondCompData;
  }
  nb::object getFinalMatchCheck() const {
    if (!fmud.pyMCSFinalMatchCheck.is_none()) {
      return fmud.pyMCSFinalMatchCheck;
    }
    return nb::none();
  }
  void setShouldAcceptMCS(nb::object mcsAcceptance) {
    PyObject *callInDict = PyDict_GetItemString(
        Py_TYPE(mcsAcceptance.ptr())->tp_dict, CALLBACK_FUNC_NAME);
    if (!callInDict || !PyCallable_Check(callInDict)) {
      throw nb::type_error("The __call__() method must be overridden in "
                           "the rdFMCS.MCSAcceptance subclass");
    }
    p->ShouldAcceptMCS = MCSAcceptancePyFunc;
    p->ShouldAcceptMCSUserData = &afud;
    afud.pyMCSAcceptance = std::move(mcsAcceptance);
    afud.pyAtomBondCompData = cfud.pyAtomBondCompData;
  }
  nb::object getShouldAcceptMCS() const {
    if (!afud.pyMCSAcceptance.is_none()) {
      return afud.pyMCSAcceptance;
    }
    return nb::none();
  }

 private:
  static bool MCSAtomComparePyFunc(const MCSAtomCompareParameters &p,
                                   const ROMol &mol1, unsigned int atom1,
                                   const ROMol &mol2, unsigned int atom2,
                                   void *userData) {
    PRECONDITION(userData, "userData must not be NULL");
    auto *cfud = static_cast<PyCompareFunctionUserData *>(userData);
    CHECK_INVARIANT(cfud, "");
    bool res = false;
    {
      PyGILStateHolder h;
      res = nb::cast<bool>(cfud->pyAtomBondCompData.pyAtomComp.attr(
          COMPARE_FUNC_NAME)(nb::cast(&p, nb::rv_policy::reference),
                             nb::cast(&mol1, nb::rv_policy::reference), atom1,
                             nb::cast(&mol2, nb::rv_policy::reference), atom2));
    }
    return res;
  }
  static bool MCSBondComparePyFunc(const MCSBondCompareParameters &p,
                                   const ROMol &mol1, unsigned int bond1,
                                   const ROMol &mol2, unsigned int bond2,
                                   void *userData) {
    PRECONDITION(userData, "userData must not be NULL");
    auto *cfud = static_cast<PyCompareFunctionUserData *>(userData);
    CHECK_INVARIANT(cfud, "");
    bool res = false;
    {
      PyGILStateHolder h;
      res = nb::cast<bool>(cfud->pyAtomBondCompData.pyBondComp.attr(
          COMPARE_FUNC_NAME)(nb::cast(&p, nb::rv_policy::reference),
                             nb::cast(&mol1, nb::rv_policy::reference), bond1,
                             nb::cast(&mol2, nb::rv_policy::reference), bond2));
    }
    return res;
  }
  static bool MCSProgressCallbackPyFunc(const MCSProgressData &stat,
                                        const MCSParameters &params,
                                        void *userData) {
    PRECONDITION(userData, "userData must not be NULL");
    auto *pcud = static_cast<PyProgressCallbackUserData *>(userData);
    CHECK_INVARIANT(pcud, "");
    bool res = false;
    MCSParameters paramsCopy(params);
    if (pcud->pyAtomBondCompData.standardBondTyperFunc) {
      paramsCopy.BondTyper = pcud->pyAtomBondCompData.standardBondTyperFunc;
    }
    {
      PyGILStateHolder h;
      PyMCSParameters ps(paramsCopy, *pcud);
      PyMCSProgressData pd(stat);
      // Use rv_policy::reference because ps and pd are stack-allocated and
      // not copy-constructible (they own unique_ptrs). The Python call is
      // synchronous so the references remain valid for its duration.
      res = nb::cast<bool>(pcud->pyMCSProgress.attr(CALLBACK_FUNC_NAME)(
          nb::cast(&pd, nb::rv_policy::reference),
          nb::cast(&ps, nb::rv_policy::reference)));
    }
    return res;
  }
  static bool MCSFinalMatchCheckPyFunc(
      const std::uint32_t c1[], const std::uint32_t c2[], const ROMol &mol1,
      const FMCS::Graph &query, const ROMol &mol2, const FMCS::Graph &target,
      const MCSParameters *params) {
    PRECONDITION(params, "params must not be NULL");
    auto *fmud = static_cast<PyMCSFinalMatchCheckFunctionUserData *>(
        params->FinalMatchCheckerUserData);
    CHECK_INVARIANT(fmud, "");
    bool res = false;
    {
      PyGILStateHolder h;
      PyMCSParameters ps(*params, *fmud);
      auto pyAtomIdxMatch = buildAtomIdxMatchTuple(c1, c2, query, target);
      auto pyBondIdxMatch =
          buildBondIdxMatchTuple(c1, c2, mol1, query, mol2, target);
      res = nb::cast<bool>(fmud->pyMCSFinalMatchCheck.attr(CALLBACK_FUNC_NAME)(
          nb::cast(&mol1, nb::rv_policy::reference),
          nb::cast(&mol2, nb::rv_policy::reference), pyAtomIdxMatch,
          pyBondIdxMatch, nb::cast(&ps, nb::rv_policy::reference)));
    }
    return res;
  }
  static bool MCSAcceptancePyFunc(
      const ROMol &query, const ROMol &target,
      const std::vector<std::pair<int, int>> &atomIdxMatch,
      const std::vector<std::pair<int, int>> &bondIdxMatch,
      const MCSParameters *params) {
    PRECONDITION(params, "params must not be NULL");
    auto *afud = static_cast<PyMCSAcceptanceFunctionUserData *>(
        params->ShouldAcceptMCSUserData);
    CHECK_INVARIANT(afud, "");
    bool res = false;
    {
      PyGILStateHolder h;
      PyMCSParameters ps(*params, *afud);
      auto pyAtomIdxMatch = convertMatchesToTupleOfPairs(atomIdxMatch);
      auto pyBondIdxMatch = convertMatchesToTupleOfPairs(bondIdxMatch);
      res = nb::cast<bool>(afud->pyMCSAcceptance.attr(CALLBACK_FUNC_NAME)(
          nb::cast(&query, nb::rv_policy::reference),
          nb::cast(&target, nb::rv_policy::reference), pyAtomIdxMatch,
          pyBondIdxMatch, nb::cast(&ps, nb::rv_policy::reference)));
    }
    return res;
  }

  std::unique_ptr<MCSParameters> p;
  PyCompareFunctionUserData cfud;
  PyProgressCallbackUserData pcud;
  PyMCSFinalMatchCheckFunctionUserData fmud;
  PyMCSAcceptanceFunctionUserData afud;
};

// -------------------------------------------------------------------
// FindMCS wrappers
// -------------------------------------------------------------------
MCSResult *FindMCSWrapper(nb::object mols, bool maximizeBonds, double threshold,
                          unsigned int timeout, bool verbose,
                          bool matchValences, bool ringMatchesRingOnly,
                          bool completeRingsOnly, bool matchChiralTag,
                          AtomComparator atomComp, BondComparator bondComp,
                          RingComparator ringComp, std::string seedSmarts) {
  std::vector<ROMOL_SPTR> ms;
  unsigned int nElems = nb::len(mols);
  ms.resize(nElems);
  for (unsigned int i = 0; i < nElems; ++i) {
    nb::object mol_obj = mols[i];
    if (mol_obj.is_none()) {
      throw nb::value_error("molecule is None");
    }
    ROMOL_SPTR react(nb::cast<ROMol *>(mol_obj), [](ROMol *) {});
    ms[i] = react;
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
  p.BondCompareParameters.MatchFusedRingsStrict =
      (ringComp == StrictRingFusion);

  MCSResult *res = nullptr;
  {
    NOGIL gil;
    res = new MCSResult(findMCS(ms, &p));
  }
  return res;
}

MCSResult *FindMCSWrapper2(nb::object mols, PyMCSParameters &pyMcsParams) {
  std::vector<ROMOL_SPTR> ms;
  unsigned int nElems = nb::len(mols);
  ms.resize(nElems);
  for (unsigned int i = 0; i < nElems; ++i) {
    nb::object mol_obj = mols[i];
    if (mol_obj.is_none()) {
      throw nb::value_error("molecule is None");
    }
    ROMOL_SPTR react(nb::cast<ROMol *>(mol_obj), [](ROMol *) {});
    ms[i] = react;
  }

  MCSResult *res = nullptr;
  {
    NOGIL gil;
    res = new MCSResult(findMCS(ms, pyMcsParams.get()));
  }
  return res;
}
}  // namespace RDKit

NB_MODULE(rdFMCS, m) {
  m.doc() = "Module containing a C++ implementation of the FMCS algorithm";

  // MCSResult
  nb::class_<RDKit::MCSResult>(m, "MCSResult", "used to return MCS results")
      .def_ro("numAtoms", &RDKit::MCSResult::NumAtoms,
              "number of atoms in MCS")
      .def_ro("numBonds", &RDKit::MCSResult::NumBonds,
              "number of bonds in MCS")
      .def_prop_ro(
          "queryMol",
          [](const RDKit::MCSResult &self) { return toStd(self.QueryMol); },
          "query molecule for the MCS")
      .def_ro("smartsString", &RDKit::MCSResult::SmartsString,
              "SMARTS string for the MCS")
      .def_ro("canceled", &RDKit::MCSResult::Canceled,
              "if True, the MCS calculation did not finish")
      .def_prop_ro(
          "degenerateSmartsQueryMolDict",
          [](const RDKit::MCSResult &self) {
            nb::dict res;
            for (const auto &pair : self.DegenerateSmartsQueryMolDict) {
              res[nb::cast(pair.first)] = nb::cast(toStd(pair.second));
            }
            return res;
          },
          "Dictionary collecting all degenerate (SMARTS, queryMol) pairs "
          "(empty if MCSParameters.StoreAll is False)");

  // AtomCompare enum
  nb::enum_<RDKit::AtomComparator>(m, "AtomCompare")
      .value("CompareAny", RDKit::AtomCompareAny)
      .value("CompareElements", RDKit::AtomCompareElements)
      .value("CompareIsotopes", RDKit::AtomCompareIsotopes)
      .value("CompareAnyHeavyAtom", RDKit::AtomCompareAnyHeavyAtom);

  // BondCompare enum
  nb::enum_<RDKit::BondComparator>(m, "BondCompare")
      .value("CompareAny", RDKit::BondCompareAny)
      .value("CompareOrder", RDKit::BondCompareOrder)
      .value("CompareOrderExact", RDKit::BondCompareOrderExact);

  // RingCompare enum
  nb::enum_<RDKit::RingComparator>(m, "RingCompare")
      .value("IgnoreRingFusion", RDKit::IgnoreRingFusion)
      .value("PermissiveRingFusion", RDKit::PermissiveRingFusion)
      .value("StrictRingFusion", RDKit::StrictRingFusion);

  // FindMCS (simple parameter version)
  m.def("FindMCS", RDKit::FindMCSWrapper,
        "mols"_a, "maximizeBonds"_a = true, "threshold"_a = 1.0,
        "timeout"_a = 3600, "verbose"_a = false, "matchValences"_a = false,
        "ringMatchesRingOnly"_a = false, "completeRingsOnly"_a = false,
        "matchChiralTag"_a = false,
        "atomCompare"_a = RDKit::AtomCompareElements,
        "bondCompare"_a = RDKit::BondCompareOrder,
        "ringCompare"_a = RDKit::IgnoreRingFusion, "seedSmarts"_a = "",
        "Find the MCS for a set of molecules",
        nb::rv_policy::take_ownership);

  // MCSParameters
  nb::class_<RDKit::PyMCSParameters>(
      m, "MCSParameters",
      "Parameters controlling how the MCS is constructed")
      .def(nb::init<>())
      .def_prop_rw(
          "MaximizeBonds", &RDKit::PyMCSParameters::getMaximizeBonds,
          &RDKit::PyMCSParameters::setMaximizeBonds,
          "toggles maximizing the number of bonds (instead of the "
          "number of atoms)")
      .def_prop_rw("Threshold", &RDKit::PyMCSParameters::getThreshold,
                   &RDKit::PyMCSParameters::setThreshold,
                   "fraction of the dataset that must contain the MCS")
      .def_prop_rw("Timeout", &RDKit::PyMCSParameters::getTimeout,
                   &RDKit::PyMCSParameters::setTimeout,
                   "timeout (in seconds) for the calculation")
      .def_prop_rw("Verbose", &RDKit::PyMCSParameters::getVerbose,
                   &RDKit::PyMCSParameters::setVerbose, "toggles verbose mode")
      .def_prop_rw(
          "AtomCompareParameters",
          [](RDKit::PyMCSParameters &self) -> RDKit::MCSAtomCompareParameters & {
            return const_cast<RDKit::MCSAtomCompareParameters &>(
                self.getAtomCompareParameters());
          },
          &RDKit::PyMCSParameters::setAtomCompareParameters,
          "parameters for comparing atoms")
      .def_prop_rw(
          "BondCompareParameters",
          [](RDKit::PyMCSParameters &self) -> RDKit::MCSBondCompareParameters & {
            return const_cast<RDKit::MCSBondCompareParameters &>(
                self.getBondCompareParameters());
          },
          &RDKit::PyMCSParameters::setBondCompareParameters,
          "parameters for comparing bonds")
      .def_prop_rw("AtomTyper", &RDKit::PyMCSParameters::getMCSAtomTyper,
                   &RDKit::PyMCSParameters::setMCSAtomTyper,
                   R"DOC(atom typer to be used. Must be one of the
members of the rdFMCS.AtomCompare class or
an instance of a user-defined subclass of
rdFMCS.MCSAtomCompare)DOC")
      .def_prop_rw("BondTyper", &RDKit::PyMCSParameters::getMCSBondTyper,
                   &RDKit::PyMCSParameters::setMCSBondTyper,
                   R"DOC(bond typer to be used. Must be one of the
members of the rdFMCS.BondCompare class or
an instance of a user-defined subclass of
rdFMCS.MCSBondCompare)DOC")
      .def_prop_rw("ProgressCallback",
                   &RDKit::PyMCSParameters::getMCSProgressCallback,
                   &RDKit::PyMCSParameters::setMCSProgressCallback,
                   R"DOC(progress callback class. Must be a
user-defined subclass of rdFMCS.MCSProgress)DOC")
      .def_prop_rw("FinalMatchChecker",
                   &RDKit::PyMCSParameters::getFinalMatchCheck,
                   &RDKit::PyMCSParameters::setFinalMatchCheck,
                   R"DOC(seed final match checker callback class. Must be a
user-defined subclass of rdFMCS.MCSFinalMatchCheck)DOC")
      .def_prop_rw("ShouldAcceptMCS",
                   &RDKit::PyMCSParameters::getShouldAcceptMCS,
                   &RDKit::PyMCSParameters::setShouldAcceptMCS,
                   R"DOC(MCS acceptance callback class. Must be a
user-defined subclass of rdFMCS.MCSAcceptance)DOC")
      .def_prop_rw("InitialSeed", &RDKit::PyMCSParameters::getInitialSeed,
                   &RDKit::PyMCSParameters::setInitialSeed,
                   "SMILES string to be used as the seed of the MCS")
      .def_prop_rw("StoreAll", &RDKit::PyMCSParameters::getStoreAll,
                   &RDKit::PyMCSParameters::setStoreAll,
                   "toggles storage of degenerate MCSs")
      .def("__setattr__", &safeSetattr);

  // MCSAtomCompareParameters
  nb::class_<RDKit::MCSAtomCompareParameters>(
      m, "MCSAtomCompareParameters",
      "Parameters controlling how atom-atom matching is done")
      .def(nb::init<>())
      .def_rw("MatchValences",
              &RDKit::MCSAtomCompareParameters::MatchValences,
              "include atom valences in the match")
      .def_rw("MatchChiralTag",
              &RDKit::MCSAtomCompareParameters::MatchChiralTag,
              "include atom chirality in the match")
      .def_rw("MaxDistance", &RDKit::MCSAtomCompareParameters::MaxDistance,
              "Require atoms to be within this many angstroms in 3D")
      .def_rw("MatchFormalCharge",
              &RDKit::MCSAtomCompareParameters::MatchFormalCharge,
              "include formal charge in the match")
      .def_rw("RingMatchesRingOnly",
              &RDKit::MCSAtomCompareParameters::RingMatchesRingOnly,
              "ring atoms are only allowed to match other ring atoms")
      .def_rw("CompleteRingsOnly",
              &RDKit::MCSAtomCompareParameters::CompleteRingsOnly,
              "results cannot include lone ring atoms")
      .def_rw("MatchIsotope",
              &RDKit::MCSAtomCompareParameters::MatchIsotope,
              "use isotope atom queries in MCSResults")
      .def("__setattr__", &safeSetattr);

  // MCSBondCompareParameters
  nb::class_<RDKit::MCSBondCompareParameters>(
      m, "MCSBondCompareParameters",
      "Parameters controlling how bond-bond matching is done")
      .def(nb::init<>())
      .def_rw("RingMatchesRingOnly",
              &RDKit::MCSBondCompareParameters::RingMatchesRingOnly,
              "ring bonds are only allowed to match other ring bonds")
      .def_rw("CompleteRingsOnly",
              &RDKit::MCSBondCompareParameters::CompleteRingsOnly,
              "results cannot include partial rings")
      .def_rw("MatchFusedRings",
              &RDKit::MCSBondCompareParameters::MatchFusedRings,
              R"DOC(enforce check on ring fusion, i.e. alpha-methylnaphthalene
won't match beta-methylnaphtalene, but decalin
will match cyclodecane unless MatchFusedRingsStrict is True)DOC")
      .def_rw("MatchFusedRingsStrict",
              &RDKit::MCSBondCompareParameters::MatchFusedRingsStrict,
              R"DOC(only enforced if MatchFusedRings is True; the ring fusion
must be the same in both query and target, i.e. decalin
won't match cyclodecane)DOC")
      .def_rw("MatchStereo", &RDKit::MCSBondCompareParameters::MatchStereo,
              "include bond stereo in the comparison")
      .def("__setattr__", &safeSetattr);

  // MCSProgressData
  nb::class_<RDKit::PyMCSProgressData>(m, "MCSProgressData",
                                        "Information about the MCS progress")
      .def(nb::init<>())
      .def_prop_ro("numAtoms", &RDKit::PyMCSProgressData::getNumAtoms,
                   "number of atoms in MCS")
      .def_prop_ro("numBonds", &RDKit::PyMCSProgressData::getNumBonds,
                   "number of bonds in MCS")
      .def_prop_ro("seedProcessed",
                   &RDKit::PyMCSProgressData::getSeedProcessed,
                   "number of processed seeds");

  // MCSAtomCompare base class — users subclass and override __call__
  nb::class_<RDKit::PyMCSAtomCompare>(
      m, "MCSAtomCompare",
      R"DOC(Base class. Subclass and override
MCSAtomCompare.__call__() to define custom
atom compare functions, then set MCSParameters.AtomTyper
to an instance of the subclass)DOC")
      .def(nb::init<>())
      .def("CheckAtomRingMatch", &RDKit::PyMCSAtomCompare::checkAtomRingMatch,
           "parameters"_a, "mol1"_a, "atom1"_a, "mol2"_a, "atom2"_a,
           "Return True if both atoms are, or are not, in a ring")
      .def("CheckAtomCharge", &RDKit::PyMCSAtomCompare::checkAtomCharge,
           "parameters"_a, "mol1"_a, "atom1"_a, "mol2"_a, "atom2"_a,
           "Return True if both atoms have the same formal charge")
      .def("CheckAtomChirality",
           &RDKit::PyMCSAtomCompare::checkAtomChirality,
           "parameters"_a, "mol1"_a, "atom1"_a, "mol2"_a, "atom2"_a,
           "Return True if both atoms have, or have not, a chiral tag")
      .def(COMPARE_FUNC_NAME,
           [](RDKit::PyMCSAtomCompare &,
              const RDKit::MCSAtomCompareParameters &, const RDKit::ROMol &,
              unsigned int, const RDKit::ROMol &, unsigned int) -> bool {
             PyErr_SetString(PyExc_AttributeError,
                             "The __call__() method must be overridden in the "
                             "rdFMCS.MCSAtomCompare subclass");
             throw nb::python_error();
           },
           "parameters"_a, "mol1"_a, "atom1"_a, "mol2"_a, "atom2"_a,
           "override to implement custom atom comparison");

  // MCSBondCompare base class
  nb::class_<RDKit::PyMCSBondCompare>(
      m, "MCSBondCompare",
      R"DOC(Base class. Subclass and override
MCSBondCompare.__call__() to define custom
bond compare functions, then set MCSParameters.BondTyper
to an instance of the subclass)DOC")
      .def(nb::init<>())
      .def("CheckBondStereo", &RDKit::PyMCSBondCompare::checkBondStereo,
           "parameters"_a, "mol1"_a, "bond1"_a, "mol2"_a, "bond2"_a,
           "Return True if both bonds have, or have not, a stereo descriptor")
      .def("CheckBondRingMatch", &RDKit::PyMCSBondCompare::checkBondRingMatch,
           "parameters"_a, "mol1"_a, "bond1"_a, "mol2"_a, "bond2"_a,
           "Return True if both bonds are, or are not, part of a ring")
      .def(COMPARE_FUNC_NAME,
           [](RDKit::PyMCSBondCompare &,
              const RDKit::MCSBondCompareParameters &, const RDKit::ROMol &,
              unsigned int, const RDKit::ROMol &, unsigned int) -> bool {
             PyErr_SetString(PyExc_AttributeError,
                             "The __call__() method must be overridden in the "
                             "rdFMCS.MCSBondCompare subclass");
             throw nb::python_error();
           },
           "parameters"_a, "mol1"_a, "bond1"_a, "mol2"_a, "bond2"_a,
           "override to implement custom bond comparison");

  // MCSProgress base class — users subclass and override __call__(stat, params)
  // The callback is stored as a Python object in PyMCSParameters and invoked
  // via nb::object::attr(CALLBACK_FUNC_NAME)(); no C++ base class dispatch.
  struct MCSProgress {};
  nb::class_<MCSProgress>(
      m, "MCSProgress",
      R"DOC(Base class. Subclass and override
MCSProgress.__call__()
to define a custom callback function)DOC")
      .def(nb::init<>())
      .def(CALLBACK_FUNC_NAME,
           [](MCSProgress &, const RDKit::PyMCSProgressData &,
              const RDKit::PyMCSParameters &) -> bool {
             PyErr_SetString(PyExc_AttributeError,
                             "The __call__() method must be overridden in the "
                             "rdFMCS.MCSProgress subclass");
             throw nb::python_error();
           },
           "stat"_a, "parameters"_a,
           "override to implement a custom progress callback");

  // MCSFinalMatchCheck base class
  struct MCSFinalMatchCheck {};
  nb::class_<MCSFinalMatchCheck>(
      m, "MCSFinalMatchCheck",
      R"DOC(Base class. Subclass and override
MCSFinalMatchCheck.__call__()
to define a custom boolean callback function.
Returning True will cause the growing seed to be accepted,
False to be rejected)DOC")
      .def(nb::init<>())
      .def(CALLBACK_FUNC_NAME,
           [](MCSFinalMatchCheck &) -> bool {
             PyErr_SetString(PyExc_AttributeError,
                             "The __call__() method must be overridden in the "
                             "rdFMCS.MCSFinalMatchCheck subclass");
             throw nb::python_error();
           },
           "override to implement a custom seed final match checker callback");

  // MCSAcceptance base class
  struct MCSAcceptance {};
  nb::class_<MCSAcceptance>(
      m, "MCSAcceptance",
      R"DOC(Base class. Subclass and override
MCSAcceptance.__call__()
to define a custom boolean callback function.
Returning True will cause the MCS candidate to be accepted,
False to be rejected)DOC")
      .def(nb::init<>())
      .def(CALLBACK_FUNC_NAME,
           [](MCSAcceptance &) -> bool {
             PyErr_SetString(PyExc_AttributeError,
                             "The __call__() method must be overridden in the "
                             "rdFMCS.MCSAcceptance subclass");
             throw nb::python_error();
           },
           "override to implement a custom MCS acceptance callback");

  // FindMCS (MCSParameters version)
  m.def("FindMCS", RDKit::FindMCSWrapper2, "mols"_a, "parameters"_a,
        "Find the MCS for a set of molecules",
        nb::rv_policy::take_ownership);
}
