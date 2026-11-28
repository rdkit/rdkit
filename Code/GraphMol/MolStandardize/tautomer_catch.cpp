//
//  Copyright (C) 2026 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/MolStandardize/MolStandardize.h>
#include <GraphMol/MolStandardize/Tautomer.h>

using namespace RDKit;

TEST_CASE("tautomer scoring prefilter is independent of aromaticity",
          "[tautomers]") {
  namespace TS = MolStandardize::TautomerScoringFunctions;
  MolStandardize::TautomerEnumerator te;
  // Tautomerization can change whether an atom is aromatic. The prefilter
  // must retain any term that can match a generated tautomer.
  const std::vector<std::string> inputs = {
      "N=C1NC(=O)CS1",
      "N=C1NC(=O)[CH](S1)C2=CC=CC=C2",
      "N=C1C=CC(C(=O)O)C=N1",
      "Nc1nc(O)cs1",
      "N=c1[nH]c(O)cs1",
      "OC1=NC2=C(C=C(C=C2C(=N1)O)[N+]([O-])=O)[N+]([O-])=O"};
  for (const auto &smi : inputs) {
    CAPTURE(smi);
    std::unique_ptr<RWMol> mol{SmilesToMol(smi)};
    REQUIRE(mol);
    const auto tautomers = te.enumerate(*mol);
    REQUIRE(tautomers.status() ==
            MolStandardize::TautomerEnumeratorStatus::Completed);
    REQUIRE(tautomers.size() > 1);
    const auto scorer = TS::makeOptimizedScorer(*mol);
    for (const auto &taut : tautomers) {
      CAPTURE(MolToSmiles(*taut));
      CHECK(scorer(*taut) == TS::scoreTautomer(*taut));
    }
    std::unique_ptr<ROMol> full{te.canonicalize(*mol, TS::scoreTautomer)};
    std::unique_ptr<ROMol> optimized{te.canonicalize(*mol)};
    REQUIRE(full);
    REQUIRE(optimized);
    CHECK(MolToSmiles(*optimized) == MolToSmiles(*full));
    te.canonicalizeInPlace(*mol);
    CHECK(MolToSmiles(*mol) == MolToSmiles(*full));
  }
}

TEST_CASE("exclude tautomer regions") {
  SECTION("basics") {
    MolStandardize::CleanupParameters params;
    MolStandardize::TautomerEnumerator te(params);
    auto m = "CCCCC=O"_smiles;
    REQUIRE(m);
    {  // baseline
      auto tauts = te.enumerate(*m);
      CHECK(tauts.size() == 2);
    }
    {  // blocking
      std::vector<std::vector<unsigned int>> protectedAtomsVec = {
          {3, 4, 5}, {3}, {4}, {5}, {4, 5}, {3, 5}};
      for (const auto &protectedAtoms : protectedAtomsVec) {
        ROMol mcopy(*m);
        for (auto i : protectedAtoms) {
          mcopy.getAtomWithIdx(i)->setProp(common_properties::_protected, 1);
        }
        auto tauts = te.enumerate(mcopy);
        CHECK(tauts.size() == 1);
      }
    }
    {  // blocking non-participating atoms doesn't affect anything
      ROMol mcopy(*m);
      for (auto i : {0, 1, 2}) {
        mcopy.getAtomWithIdx(i)->setProp(common_properties::_protected, 1);
      }
      auto tauts = te.enumerate(mcopy);
      CHECK(tauts.size() == 2);
    }
  }
  SECTION("canonicalization") {
    MolStandardize::CleanupParameters params;
    MolStandardize::TautomerEnumerator te(params);
    auto m = "CCCC=CO"_smiles;
    REQUIRE(m);
    {  // baseline
      std::unique_ptr<ROMol> taut{te.canonicalize(*m)};
      CHECK(taut);
      CHECK(MolToSmiles(*taut) == "CCCCC=O");
    }
    {  // blocking participating atoms prevents the SMILES from changing
      ROMol mcopy(*m);
      for (auto i : {3, 4, 5}) {
        mcopy.getAtomWithIdx(i)->setProp(common_properties::_protected, 1);
      }
      std::unique_ptr<ROMol> taut{te.canonicalize(mcopy)};
      CHECK(taut);
      CHECK(MolToSmiles(*taut) == "CCCC=CO");
    }
  }
}