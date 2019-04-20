//
//  Copyright (C) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MolStandardize/MolStandardize.h>
#include <GraphMol/MolStandardize/Fragment.h>
#include <GraphMol/MolStandardize/Charge.h>

#include <iostream>
#include <fstream>

using namespace RDKit;

TEST_CASE("SKIP_IF_ALL_MATCH") {
  auto m = "[Na+].[Cl-]"_smiles;
  REQUIRE(m);

  SECTION("default") {
    MolStandardize::FragmentRemover fragRemover;
    std::unique_ptr<ROMol> outm(fragRemover.remove(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "[Na+]");
  }
  SECTION("don't remove all") {
    MolStandardize::FragmentRemover fragRemover("", true, true);
    std::unique_ptr<ROMol> outm(fragRemover.remove(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "[Cl-].[Na+]");
  }
  SECTION("feel free to remove everything") {
    MolStandardize::FragmentRemover fragRemover("", false, false);
    std::unique_ptr<ROMol> outm(fragRemover.remove(*m));
    REQUIRE(outm);
    CHECK(outm->getNumAtoms() == 0);
  }
  SECTION("don't remove all 2") {
    MolStandardize::FragmentRemover fragRemover("", true, true);
    auto m = "[Na+].[Cl-].[Na+].[Cl-]"_smiles;
    REQUIRE(m);
    std::unique_ptr<ROMol> outm(fragRemover.remove(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "[Cl-].[Cl-].[Na+].[Na+]");
  }
}

TEST_CASE("symmetry in the uncharger") {
  SECTION("case 1") {
    auto m = "C[N+](C)(C)CC(C(=O)[O-])CC(=O)[O-]"_smiles;
    REQUIRE(m);
    {
      bool canonicalOrdering = false;
      MolStandardize::Uncharger uncharger(canonicalOrdering);
      std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
      REQUIRE(outm);
      CHECK(MolToSmiles(*outm) == "C[N+](C)(C)CC(CC(=O)[O-])C(=O)O");
    }
    {
      bool canonicalOrdering = true;
      MolStandardize::Uncharger uncharger(canonicalOrdering);
      std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
      REQUIRE(outm);
      CHECK(MolToSmiles(*outm) == "C[N+](C)(C)CC(CC(=O)O)C(=O)[O-]");
    }
    {
      MolStandardize::CleanupParameters params;
      std::unique_ptr<ROMol> outm(MolStandardize::chargeParent(*m, params));
      REQUIRE(outm);
      CHECK(MolToSmiles(*outm) == "C[N+](C)(C)CC(CC(=O)O)C(=O)[O-]");
    }
    {
      MolStandardize::CleanupParameters params;
      params.doCanonical = false;
      std::unique_ptr<ROMol> outm(MolStandardize::chargeParent(*m, params));
      REQUIRE(outm);
      CHECK(MolToSmiles(*outm) == "C[N+](C)(C)CC(CC(=O)[O-])C(=O)O");
    }
  }
}

TEST_CASE("uncharger bug with duplicates") {
  SECTION("case 1") {
    auto m = "[NH3+]CC([O-])C[O-]"_smiles;
    REQUIRE(m);
    MolStandardize::Uncharger uncharger;
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "NCC(O)CO");
  }
  SECTION("case 2") {
    auto m = "CC([O-])C[O-].[Na+]"_smiles;
    REQUIRE(m);
    MolStandardize::Uncharger uncharger;
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "CC([O-])CO.[Na+]");
  }
  SECTION("acids + others 1, github #2392") {
    auto m = "C[N+](C)(C)CC(C[O-])CC(=O)[O-]"_smiles;
    REQUIRE(m);
    bool doCanonical = false;
    MolStandardize::Uncharger uncharger(doCanonical);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "C[N+](C)(C)CC(CO)CC(=O)[O-]");
  }
  SECTION("acids + others 2, github #2392") {
    auto m = "C[N+](C)(C)CC(CC(=O)[O-])C[O-]"_smiles;
    REQUIRE(m);
    bool doCanonical = false;
    MolStandardize::Uncharger uncharger(doCanonical);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "C[N+](C)(C)CC(CO)CC(=O)[O-]");
  }
}

TEST_CASE(
    "github #2411: MolStandardize: FragmentRemover should not sanitize "
    "fragments ") {
  SECTION("demo") {
    std::string smi = "CN(C)(C)C.Cl";
    bool debugParse = false;
    bool sanitize = false;
    std::unique_ptr<ROMol> m(SmilesToMol(smi, debugParse, sanitize));
    REQUIRE(m);

    MolStandardize::FragmentRemover fragRemover;
    std::unique_ptr<ROMol> outm(fragRemover.remove(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "CN(C)(C)C");
  }
}