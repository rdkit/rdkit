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
#include <GraphMol/MolStandardize/Normalize.h>
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

TEST_CASE("symmetry in the uncharger", "uncharger") {
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

TEST_CASE("uncharger bug with duplicates", "uncharger") {
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
    "fragments") {
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

TEST_CASE(
    "github #2452: incorrectly removing charge from boron anions"
    "fragments,uncharger") {
  SECTION("demo") {
    auto m = "C[B-](C)(C)C"_smiles;
    REQUIRE(m);
    bool canonicalOrdering = true;

    MolStandardize::Uncharger uncharger(canonicalOrdering);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(1)->getFormalCharge() == -1);
    CHECK(MolToSmiles(*outm) == "C[B-](C)(C)C");
  }
  SECTION("should be removed") {
    auto m = "C[BH-](C)(C)"_smiles;
    REQUIRE(m);
    bool canonicalOrdering = true;

    MolStandardize::Uncharger uncharger(canonicalOrdering);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(1)->getFormalCharge() == 0);
    CHECK(MolToSmiles(*outm) == "CB(C)C");
  }
}

TEST_CASE("github #2602: Uncharger ignores dications", "uncharger") {
  SECTION("demo") {
    auto m = "[O-]CCC[O-].[Ca+2]"_smiles;
    REQUIRE(m);
    bool canonicalOrdering = true;
    MolStandardize::Uncharger uncharger(canonicalOrdering);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(5)->getFormalCharge() == 2);
    CHECK(outm->getAtomWithIdx(0)->getFormalCharge() == -1);
    CHECK(outm->getAtomWithIdx(4)->getFormalCharge() == -1);
    CHECK(MolToSmiles(*outm) == "[Ca+2].[O-]CCC[O-]");
  }
}

TEST_CASE(
    "github #2605: Uncharger incorrectly neutralizes cations when "
    "non-neutralizable anions are present.",
    "uncharger") {
  SECTION("demo") {
    auto m = "F[B-](F)(F)F.[NH3+]CCC"_smiles;
    REQUIRE(m);
    bool canonicalOrdering = true;
    MolStandardize::Uncharger uncharger(canonicalOrdering);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(1)->getFormalCharge() == -1);
    CHECK(outm->getAtomWithIdx(5)->getFormalCharge() == 1);
    CHECK(MolToSmiles(*outm) == "CCC[NH3+].F[B-](F)(F)F");
  }
  SECTION("multiple positively charged sites") {
    auto m = "F[B-](F)(F)F.[NH3+]CC=C[NH3+]"_smiles;
    REQUIRE(m);
    bool canonicalOrdering = true;
    MolStandardize::Uncharger uncharger(canonicalOrdering);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(1)->getFormalCharge() == -1);
    CHECK(outm->getAtomWithIdx(5)->getFormalCharge() == 0);
    CHECK(outm->getAtomWithIdx(9)->getFormalCharge() == 1);
    CHECK(MolToSmiles(*outm) == "F[B-](F)(F)F.NCC=C[NH3+]");
  }
  SECTION("make sure we don't go too far") {
    auto m = "F[B-](F)(F)F.[NH4+2]CCC"_smiles;  // totally bogus structure
    REQUIRE(m);
    bool canonicalOrdering = true;
    MolStandardize::Uncharger uncharger(canonicalOrdering);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(1)->getFormalCharge() == -1);
    CHECK(outm->getAtomWithIdx(5)->getFormalCharge() == 1);
    CHECK(MolToSmiles(*outm) == "CCC[NH3+].F[B-](F)(F)F");
  }
}

TEST_CASE("github #2610: Uncharger incorrectly modifying a zwitterion.",
          "uncharger") {
  SECTION("demo") {
    auto m = "C1=CC=CC[NH+]1-[O-]"_smiles;
    REQUIRE(m);
    bool canonicalOrdering = true;
    MolStandardize::Uncharger uncharger(canonicalOrdering);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(5)->getFormalCharge() == 1);
    CHECK(outm->getAtomWithIdx(6)->getFormalCharge() == -1);
    CHECK(MolToSmiles(*outm) == "[O-][NH+]1C=CC=CC1");
  }
}

TEST_CASE("problems with ringInfo initialization", "normalizer") {
  std::string tfs =
      R"TXT(Bad amide tautomer1	[C:1]([OH1;D1:2])=;!@[NH1:3]>>[C:1](=[OH0:2])-[NH2:3]
Bad amide tautomer2	[C:1]([OH1;D1:2])=;!@[NH0:3]>>[C:1](=[OH0:2])-[NH1:3])TXT";
  std::stringstream iss(tfs);
  MolStandardize::Normalizer nrml(iss, 20);
  SECTION("example1") {
    auto m = "Cl.Cl.OC(=N)NCCCCCCCCCCCCNC(O)=N"_smiles;
    REQUIRE(m);
    std::unique_ptr<ROMol> res(nrml.normalize(*m));
    REQUIRE(res);
    CHECK(MolToSmiles(*res) == "Cl.Cl.NC(=O)NCCCCCCCCCCCCNC(N)=O");
  }
}