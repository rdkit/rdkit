//
//  Copyright (C) 2015-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <RDBoost/boost_shared_ptr.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/MMPA/MMPA.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

nb::tuple buildFragmentResults(
    const std::vector<std::pair<RDKit::ROMOL_SPTR, RDKit::ROMOL_SPTR>> &tres,
    bool resultsAsMols) {
  nb::list pyres;
  for (const auto &pr : tres) {
    nb::list lres;
    if (resultsAsMols) {
      if (pr.first) {
        lres.append(pr.first);
      } else {
        lres.append(nb::none());
      }
      lres.append(pr.second);
    } else {
      if (pr.first) {
        lres.append(RDKit::MolToSmiles(*(pr.first), true));
      } else {
        lres.append(std::string(""));
      }
      lres.append(RDKit::MolToSmiles(*(pr.second), true));
    }
    pyres.append(nb::tuple(lres));
  }
  return nb::tuple(pyres);
}

nb::tuple fragmentMolHelper(const RDKit::ROMol &mol, unsigned int maxCuts,
                            unsigned int maxCutBonds,
                            const std::string &pattern, bool resultsAsMols) {
  std::vector<std::pair<RDKit::ROMOL_SPTR, RDKit::ROMOL_SPTR>> tres;
  bool ok = RDKit::MMPA::fragmentMol(mol, tres, maxCuts, maxCutBonds, pattern);
  if (!ok) {
    return nb::tuple();
  }
  return buildFragmentResults(tres, resultsAsMols);
}

nb::tuple fragmentMolHelper2(const RDKit::ROMol &mol, unsigned int minCuts,
                             unsigned int maxCuts, unsigned int maxCutBonds,
                             const std::string &pattern, bool resultsAsMols) {
  std::vector<std::pair<RDKit::ROMOL_SPTR, RDKit::ROMOL_SPTR>> tres;
  bool ok = RDKit::MMPA::fragmentMol(mol, tres, minCuts, maxCuts, maxCutBonds,
                                     pattern);
  if (!ok) {
    return nb::tuple();
  }
  return buildFragmentResults(tres, resultsAsMols);
}

nb::tuple fragmentMolHelper3(const RDKit::ROMol &mol,
                             const std::vector<unsigned int> &bondsToCut,
                             unsigned int minCuts, unsigned int maxCuts,
                             bool resultsAsMols) {
  if (bondsToCut.empty()) {
    throw nb::value_error("bondsToCut must be non-empty");
  }
  std::vector<std::pair<RDKit::ROMOL_SPTR, RDKit::ROMOL_SPTR>> tres;
  bool ok = RDKit::MMPA::fragmentMol(mol, tres, bondsToCut, minCuts, maxCuts);
  if (!ok) {
    return nb::tuple();
  }
  return buildFragmentResults(tres, resultsAsMols);
}

}  // namespace

NB_MODULE(rdMMPA, m) {
  m.doc() = "Module containing a C++ implementation of code for doing MMPA";

  m.def("FragmentMol", fragmentMolHelper,
        "mol"_a, "maxCuts"_a = 3, "maxCutBonds"_a = 20,
        "pattern"_a = "[#6+0;!$(*=,#[!#6])]!@!=!#[*]",
        "resultsAsMols"_a = true,
        R"DOC(Does the fragmentation necessary for an MMPA analysis)DOC");

  m.def("FragmentMol", fragmentMolHelper2,
        "mol"_a, "minCuts"_a, "maxCuts"_a, "maxCutBonds"_a,
        "pattern"_a = "[#6+0;!$(*=,#[!#6])]!@!=!#[*]",
        "resultsAsMols"_a = true,
        R"DOC(Does the fragmentation necessary for an MMPA analysis)DOC");

  m.def("FragmentMol", fragmentMolHelper3,
        "mol"_a, "bondsToCut"_a, "minCuts"_a = 1, "maxCuts"_a = 3,
        "resultsAsMols"_a = true,
        R"DOC(Does the fragmentation necessary for an MMPA analysis)DOC");
}
