//
//  Copyright (C) 2019 Greg Landrum and T5 Informatics GmbH
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "RDGeneral/test.h"
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/ScaffoldNetwork/ScaffoldNetwork.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

// declarations of stuff we'll be testing that isn't in the public API
namespace RDKit {
namespace ScaffoldNetwork {
namespace detail {
ROMol *makeScaffoldGeneric(const ROMol &mol, bool doAtoms, bool doBonds);
ROMol *removeAttachmentPoints(const ROMol &mol,
                              const ScaffoldNetworkParams &params);
ROMol *pruneMol(const ROMol &mol, const ScaffoldNetworkParams &params);
ROMol *flattenMol(const ROMol &mol, const ScaffoldNetworkParams &params);
}  // namespace detail
}  // namespace ScaffoldNetwork
}  // namespace RDKit

TEST_CASE("flattenMol", "[unittest, scaffolds]") {
  auto m = "Cl.[13CH3][C@H](F)/C=C/C"_smiles;
  REQUIRE(m);
  SECTION("defaults") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    std::unique_ptr<ROMol> pm(ScaffoldNetwork::detail::flattenMol(*m, ps));
    REQUIRE(pm);
    auto smiles = MolToSmiles(*pm);
    CHECK(smiles == "CC=CC(C)F");
  }
  SECTION("isotopes") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.flattenIsotopes = false;
    std::unique_ptr<ROMol> pm(ScaffoldNetwork::detail::flattenMol(*m, ps));
    REQUIRE(pm);
    auto smiles = MolToSmiles(*pm);
    CHECK(smiles == "CC=CC([13CH3])F");
  }
  SECTION("chirality") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.flattenChirality = false;
    std::unique_ptr<ROMol> pm(ScaffoldNetwork::detail::flattenMol(*m, ps));
    REQUIRE(pm);
    auto smiles = MolToSmiles(*pm);
    CHECK(smiles == "C/C=C/[C@H](C)F");
  }
  SECTION("chirality and isotopes") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.flattenChirality = false;
    ps.flattenIsotopes = false;
    std::unique_ptr<ROMol> pm(ScaffoldNetwork::detail::flattenMol(*m, ps));
    REQUIRE(pm);
    auto smiles = MolToSmiles(*pm);
    CHECK(smiles == "C/C=C/[C@H]([13CH3])F");
  }
  SECTION("keep largest") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.flattenKeepLargest = false;
    std::unique_ptr<ROMol> pm(ScaffoldNetwork::detail::flattenMol(*m, ps));
    REQUIRE(pm);
    auto smiles = MolToSmiles(*pm);
    CHECK(smiles == "CC=CC(C)F.Cl");
  }
  SECTION("turn everything off") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.flattenChirality = false;
    ps.flattenIsotopes = false;
    ps.flattenKeepLargest = false;
    std::unique_ptr<ROMol> pm(ScaffoldNetwork::detail::flattenMol(*m, ps));
    REQUIRE(pm);
    auto smiles = MolToSmiles(*pm);
    CHECK(smiles == "C/C=C/[C@H]([13CH3])F.Cl");
  }
}

TEST_CASE("pruneMol", "[unittest, scaffolds]") {
  auto m = "O=C(O)C1C(=O)CC1"_smiles;
  REQUIRE(m);
  SECTION("defaults") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    std::unique_ptr<ROMol> pm(ScaffoldNetwork::detail::pruneMol(*m, ps));
    REQUIRE(pm);
    auto smiles = MolToSmiles(*pm);
    CHECK(smiles == "O=C1CCC1");
  }
}

TEST_CASE("removeAttachmentPoints", "[unittest, scaffolds]") {
  auto m = "*c1ccc(*)c*1"_smiles;
  REQUIRE(m);
  SECTION("defaults") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    std::unique_ptr<ROMol> pm(
        ScaffoldNetwork::detail::removeAttachmentPoints(*m, ps));
    REQUIRE(pm);
    CHECK(m->getNumAtoms() == 8);
    CHECK(pm->getNumAtoms() == 6);
  }
}

TEST_CASE("makeScaffoldGeneric", "[unittest, scaffolds]") {
  auto m = "c1[nH]ccc1"_smiles;
  REQUIRE(m);
  SECTION("atoms") {
    std::unique_ptr<ROMol> pm(
        ScaffoldNetwork::detail::makeScaffoldGeneric(*m, true, false));
    REQUIRE(pm);
    auto smiles = MolToSmiles(*pm);
    CHECK(smiles == "*1:*:*:*:*:1");
  }
  SECTION("bonds") {
    std::unique_ptr<ROMol> pm(
        ScaffoldNetwork::detail::makeScaffoldGeneric(*m, false, true));
    REQUIRE(pm);
    auto smiles = MolToSmiles(*pm);
    CHECK(smiles == "C1CCNC1");
  }
  SECTION("both") {
    std::unique_ptr<ROMol> pm(
        ScaffoldNetwork::detail::makeScaffoldGeneric(*m, true, true));
    REQUIRE(pm);
    auto smiles = MolToSmiles(*pm);
    CHECK(smiles == "*1****1");
  }
}

TEST_CASE("Network defaults", "[scaffolds]") {
  SECTION("basics") {}
}
