//
//  Copyright (C) 2026 Niels Maeder and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>
#include <string>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include "Embedder.h"

namespace {
void runTest(const bool legacy) {
  auto m = RDKit::v2::SmilesParse::MolFromSmiles(std::string(2000, 'C'));
  REQUIRE(m);
  CHECK(m->getNumAtoms() == 2000);
  RDKit::DGeomHelpers::EmbedParameters params{
      .randomSeed = 0xf00d, .useLegacyImplementation = legacy};
  int cid = RDKit::DGeomHelpers::EmbedMolecule(*m, params);
  CHECK(cid >= 0);
}
}  // namespace

const std::string rdbase = getenv("RDBASE");

TEST_CASE("testGithub3019Legacy") { runTest(true); }
TEST_CASE("testGithub3019AIO") { runTest(false); }
