//
// Copyright (C) David Cosgrove and other RDKit contributors 2026.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <GraphMol/SmilesParse/SmilesParse.h>

#include "FMCS.h"

#include <catch2/catch_all.hpp>

using namespace RDKit;

TEST_CASE("Github9034") {
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      "C1CCCCC1",
      "C1CCCC1",
  };

  for (auto &i : smi) {
    mols.emplace_back(v2::SmilesParse::MolFromSmiles(i));
  }
  MCSParameters p;
  p.Verbose = true;
  p.StoreAll = true;
  auto res = findMCS(mols, &p);
  // The bug was a crash caused by p.Verbose = true.  So the real test
  // is whether this completes correctly, but check a few things.
  CHECK(res.DegenerateSmartsQueryMolDict.size() == 4);
  for (const auto &r : res.DegenerateSmartsQueryMolDict) {
    CHECK(r.second->getNumAtoms() == 5);
  }
}