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
          {3, 4, 5}, {3}, {4, 5}, {3, 5}};
      for (const auto &protectedAtoms : protectedAtomsVec) {
        ROMol mcopy(*m);
        for (auto i : protectedAtoms) {
          mcopy.getAtomWithIdx(i)->setProp("_protected", 1);
        }
        auto tauts = te.enumerate(mcopy);
        CHECK(tauts.size() == 1);
      }
    }
    {  // blocking non-participating atoms
      ROMol mcopy(*m);
      for (auto i : {0, 1, 2, 4}) {
        mcopy.getAtomWithIdx(i)->setProp("_protected", 1);
      }
      auto tauts = te.enumerate(mcopy);
      CHECK(tauts.size() == 2);
    }
  }
}