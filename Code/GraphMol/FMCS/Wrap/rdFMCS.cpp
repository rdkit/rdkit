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

namespace python = boost::python;

namespace RDKit {

void SetMCSAtomTyper(MCSParameters &p, AtomComparator atomComp) {
  switch (atomComp) {
    case AtomCompareAny:
      p.AtomTyper = MCSAtomCompareAny;
      break;
    case AtomCompareElements:
      p.AtomTyper = MCSAtomCompareElements;
      break;
    case AtomCompareIsotopes:
      p.AtomTyper = MCSAtomCompareIsotopes;
      break;
    case AtomCompareAnyHeavyAtom:
      p.AtomTyper = MCSAtomCompareAnyHeavyAtom;
      break;
  }
}

void SetMCSBondTyper(MCSParameters &p, BondComparator bondComp) {
  switch (bondComp) {
    case BondCompareAny:
      p.BondTyper = MCSBondCompareAny;
      break;
    case BondCompareOrder:
      p.BondTyper = MCSBondCompareOrder;
      break;
    case BondCompareOrderExact:
      p.BondTyper = MCSBondCompareOrderExact;
      break;
  }
}

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
    if (!mols[i]) throw_value_error("molecule is None");
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
  SetMCSAtomTyper(p, atomComp);
  SetMCSBondTyper(p, bondComp);
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

MCSResult *FindMCSWrapper2(python::object mols, const MCSParameters &params) {
  std::vector<ROMOL_SPTR> ms;
  unsigned int nElems = python::extract<unsigned int>(mols.attr("__len__")());
  ms.resize(nElems);
  for (unsigned int i = 0; i < nElems; ++i) {
    if (!mols[i]) throw_value_error("molecule is None");
    ms[i] = python::extract<ROMOL_SPTR>(mols[i]);
  }

  MCSResult *res = nullptr;
  {
    NOGIL gil;
    res = new MCSResult(findMCS(ms, &params));
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

  python::class_<RDKit::MCSParameters, boost::noncopyable>(
      "MCSParameters", "Parameters controlling how the MCS is constructed")
      .def_readwrite("MaximizeBonds", &RDKit::MCSParameters::MaximizeBonds,
                     "toggles maximizing the number of bonds (instead of the "
                     "number of atoms)")
      .def_readwrite("Threshold", &RDKit::MCSParameters::Threshold,
                     "fraction of the dataset that must contain the MCS")
      .def_readwrite("Timeout", &RDKit::MCSParameters::Timeout,
                     "timeout (in seconds) for the calculation")
      .def_readwrite("Verbose", &RDKit::MCSParameters::Verbose,
                     "toggles verbose mode")
      .def_readwrite("AtomCompareParameters",
                     &RDKit::MCSParameters::AtomCompareParameters,
                     "parameters for comparing atoms")
      .def_readwrite("BondCompareParameters",
                     &RDKit::MCSParameters::BondCompareParameters,
                     "parameters for comparing bonds")
      // haven't been able to get these properly working
      // .def_readwrite("AtomTyper", &RDKit::MCSParameters::AtomTyper,
      //                "function for comparing atoms")
      // .def_readwrite("BondTyper", &RDKit::MCSParameters::BondTyper,
      //                "function for comparing bonds")
      .def_readwrite("InitialSeed", &RDKit::MCSParameters::InitialSeed,
                     "SMILES string to be used as the seed of the MCS")
      .def("SetAtomTyper", RDKit::SetMCSAtomTyper,
           (python::arg("self"), python::arg("comparator")),
           "sets the atom typer to be used. The argument should be one of the "
           "members of the rdFMCS.AtomCompare class.")
      .def("SetBondTyper", RDKit::SetMCSBondTyper,
           (python::arg("self"), python::arg("comparator")),
           "sets the bond typer to be used. The argument should be one of the "
           "members of the rdFMCS.BondCompare class.");

  ;
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
                     "ring atoms are only allowed to match other ring atoms");
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

  docString = "Find the MCS for a set of molecules";
  python::def("FindMCS", RDKit::FindMCSWrapper2,
              (python::arg("mols"), python::arg("parameters")),
              python::return_value_policy<python::manage_new_object>(),
              docString.c_str());
}
