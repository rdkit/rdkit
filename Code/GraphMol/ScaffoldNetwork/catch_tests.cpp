//
//  Copyright (C) 2019 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "catch.hpp"
#include "GraphMol/ScaffoldNetwork/detail.h"
#include "RDGeneral/test.h"
#include <sstream>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ScaffoldNetwork/ScaffoldNetwork.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>

using namespace RDKit;

#if 1
TEST_CASE("flattenMol", "[unittest][scaffolds]") {
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

TEST_CASE("pruneMol", "[unittest][scaffolds]") {
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

TEST_CASE("removeAttachmentPoints", "[unittest][scaffolds]") {
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

TEST_CASE("makeScaffoldGeneric", "[unittest][scaffolds]") {
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

TEST_CASE("getMolFragments", "[unittest][scaffolds]") {
  auto m = "c1ccccc1CC1NC(=O)CCC1"_smiles;
  REQUIRE(m);
  SECTION("defaults") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    auto frags = ScaffoldNetwork::detail::getMolFragments(*m, ps);
    REQUIRE(frags.size() == 2);
    CHECK(frags[0].first == "O=C1CCCC(Cc2ccccc2)N1");
    CHECK(frags[1].first == "O=C1CCCC(Cc2ccccc2)N1");

    auto smi1 = MolToSmiles(*frags[0].second);
    auto smi2 = MolToSmiles(*frags[1].second);
    // don't want to make any assumptions about the order in which the
    // fragments come back:
    CHECK((std::min(smi1, smi2) == "*C1CCCC(=O)N1"));
    CHECK((std::max(smi1, smi2) == "*c1ccccc1"));
  }
  SECTION("keep-linkers") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.keepOnlyFirstFragment = false;
    auto frags = ScaffoldNetwork::detail::getMolFragments(*m, ps);
    REQUIRE(frags.size() == 8);

    std::vector<std::pair<std::string, std::string>> res;
    res.reserve(frags.size());
    for (const auto &frag : frags) {
      res.push_back(std::make_pair(frag.first, MolToSmiles(*frag.second)));
    }
    std::sort(res.begin(), res.end());
    CHECK(res[0].first == "*CC1CCCC(=O)N1");
    CHECK(res[0].second == "*C*");
    CHECK(res[1].first == "*CC1CCCC(=O)N1");
    CHECK(res[1].second == "*C1CCCC(=O)N1");
    CHECK(res[5].first == "O=C1CCCC(Cc2ccccc2)N1");
    CHECK(res[5].second == "*CC1CCCC(=O)N1");
  }
  SECTION("includeScaffoldsWithAttachments=false") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.includeScaffoldsWithAttachments = false;
    auto frags = ScaffoldNetwork::detail::getMolFragments(*m, ps);
    REQUIRE(frags.size() == 2);
    CHECK(frags[0].first == "O=C1CCCC(Cc2ccccc2)N1");
    CHECK(frags[1].first == "O=C1CCCC(Cc2ccccc2)N1");

    auto smi1 = MolToSmiles(*frags[0].second);
    auto smi2 = MolToSmiles(*frags[1].second);
    // don't want to make any assumptions about the order in which the
    // fragments come back:
    CHECK((std::min(smi1, smi2) == "O=C1CCCCN1"));
    CHECK((std::max(smi1, smi2) == "c1ccccc1"));
  }
}

TEST_CASE("addMolToNetwork", "[unittest][scaffolds]") {
  SECTION("defaults") {
    auto m = "c1ccccc1CC1NC(=O)CCC1"_smiles;
    REQUIRE(m);
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ScaffoldNetwork::ScaffoldNetwork net;
    ScaffoldNetwork::detail::addMolToNetwork(*m, net, ps);
    CHECK(net.nodes.size() == 9);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 8);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Fragment;
                        }) == 2);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Generic;
                        }) == 4);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type ==
                                 ScaffoldNetwork::EdgeType::RemoveAttachment;
                        }) == 2);
    CHECK(rdcast<size_t>(std::count(net.counts.begin(), net.counts.end(), 1)) ==
          net.counts.size());

    // make sure adding the same molecule again doesn't do anything except
    // change the counts:
    ScaffoldNetwork::detail::addMolToNetwork(*m, net, ps);
    CHECK(net.nodes.size() == 9);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 8);
    CHECK(rdcast<size_t>(std::count(net.counts.begin(), net.counts.end(), 2)) ==
          net.counts.size());
  }
  SECTION("flucloxacillin") {
    auto m =
        "Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12"_smiles;
    REQUIRE(m);
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.includeGenericScaffolds = false;
    ps.includeScaffoldsWithoutAttachments = false;
    ScaffoldNetwork::ScaffoldNetwork net;
    ScaffoldNetwork::detail::addMolToNetwork(*m, net, ps);
    CHECK(net.nodes.size() == 7);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 9);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Fragment;
                        }) == 8);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type ==
                                 ScaffoldNetwork::EdgeType::Initialize;
                        }) == 1);
  }

  SECTION("generic flattened structures") {
    auto m = "Cc1ccccc1OC1C(C)C1"_smiles;
    REQUIRE(m);
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.includeGenericScaffolds = true;
    ps.includeScaffoldsWithoutAttachments = false;
    ps.keepOnlyFirstFragment = true;
    ScaffoldNetwork::ScaffoldNetwork net;
    ScaffoldNetwork::detail::addMolToNetwork(*m, net, ps);
    CHECK(net.nodes.size() == 7);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 6);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Fragment;
                        }) == 2);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type ==
                                 ScaffoldNetwork::EdgeType::Initialize;
                        }) == 1);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Generic;
                        }) == 3);
    // std::copy(net.nodes.begin(), net.nodes.end(),
    //           std::ostream_iterator<std::string>(std::cerr, " "));
    // std::cerr << std::endl;
    CHECK(std::find(net.nodes.begin(), net.nodes.end(),
                    "*1**1**1:*:*:*:*:*:1") != net.nodes.end());
  }
}
TEST_CASE("Network defaults", "[scaffolds]") {
  auto smis = {"c1ccccc1CC1NC(=O)CCC1", "c1cccnc1CC1NC(=O)CCC1"};
  std::vector<ROMOL_SPTR> ms;
  for (const auto smi : smis) {
    auto m = SmilesToMol(smi);
    REQUIRE(m);
    ms.push_back(ROMOL_SPTR(m));
  }
  SECTION("basics") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ScaffoldNetwork::ScaffoldNetwork net;
    ScaffoldNetwork::updateScaffoldNetwork(ms, net, ps);
    CHECK(net.nodes.size() == 12);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 13);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Fragment;
                        }) == 4);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Generic;
                        }) == 6);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type ==
                                 ScaffoldNetwork::EdgeType::RemoveAttachment;
                        }) == 3);
  }
  SECTION("don't remove attachments (makes sure parameters actually work)") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.includeScaffoldsWithoutAttachments = false;
    ScaffoldNetwork::ScaffoldNetwork net;
    ScaffoldNetwork::updateScaffoldNetwork(ms, net, ps);
    CHECK(net.nodes.size() == 7);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 7);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Fragment;
                        }) == 4);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Generic;
                        }) == 3);
  }
  SECTION("create network basics") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ScaffoldNetwork::ScaffoldNetwork net =
        ScaffoldNetwork::createScaffoldNetwork(ms, ps);
    CHECK(net.nodes.size() == 12);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 13);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Fragment;
                        }) == 4);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Generic;
                        }) == 6);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type ==
                                 ScaffoldNetwork::EdgeType::RemoveAttachment;
                        }) == 3);
  }
}
TEST_CASE("ostream integration", "[scaffolds]") {
  auto smis = {"c1ccccc1CC1NC(=O)CCC1"};
  std::vector<ROMOL_SPTR> ms;
  for (const auto smi : smis) {
    auto m = SmilesToMol(smi);
    REQUIRE(m);
    ms.push_back(ROMOL_SPTR(m));
  }
  SECTION("edges") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ScaffoldNetwork::ScaffoldNetwork net =
        ScaffoldNetwork::createScaffoldNetwork(ms, ps);
    CHECK(net.edges.size() == 8);
    CHECK(net.edges[0].beginIdx == 0);
    CHECK(net.edges[0].endIdx == 1);
    CHECK(net.edges[0].type == ScaffoldNetwork::EdgeType::Fragment);

    std::ostringstream oss;
    oss << net.edges[0];
    auto txt = oss.str();
    CHECK(txt == "NetworkEdge( 0->1, type:Fragment )");
  }
}

TEST_CASE("no attachment points", "[unittest][scaffolds]") {
  auto m = "c1ccccc1CC1NC(=O)CCC1"_smiles;
  REQUIRE(m);
  SECTION("others default") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.includeScaffoldsWithAttachments = false;
    ScaffoldNetwork::ScaffoldNetwork net;
    ScaffoldNetwork::detail::addMolToNetwork(*m, net, ps);
    CHECK(net.nodes.size() == 5);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 4);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Fragment;
                        }) == 2);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Generic;
                        }) == 2);
    CHECK(rdcast<size_t>(std::count(net.counts.begin(), net.counts.end(), 1)) ==
          net.counts.size());
  }
  SECTION("no generic") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.includeScaffoldsWithAttachments = false;
    ps.includeGenericScaffolds = false;
    ScaffoldNetwork::ScaffoldNetwork net;
    ScaffoldNetwork::detail::addMolToNetwork(*m, net, ps);
    CHECK(net.nodes.size() == 3);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 2);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Fragment;
                        }) == 2);
    CHECK(rdcast<size_t>(std::count(net.counts.begin(), net.counts.end(), 1)) ==
          net.counts.size());
  }
  SECTION("generic bonds") {
    auto m = "Cc1ccccc1OC1C(C)C1"_smiles;
    REQUIRE(m);
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.includeGenericScaffolds = true;
    ps.includeGenericBondScaffolds = true;
    ps.includeScaffoldsWithoutAttachments = false;
    ps.keepOnlyFirstFragment = true;
    ScaffoldNetwork::ScaffoldNetwork net;
    ScaffoldNetwork::detail::addMolToNetwork(*m, net, ps);
    CHECK(net.nodes.size() == 9);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 8);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Fragment;
                        }) == 2);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type ==
                                 ScaffoldNetwork::EdgeType::Initialize;
                        }) == 1);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Generic;
                        }) == 3);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type ==
                                 ScaffoldNetwork::EdgeType::GenericBond;
                        }) == 2);
    // std::copy(net.nodes.begin(), net.nodes.end(),
    //           std::ostream_iterator<std::string>(std::cerr, " "));
    // std::cerr << std::endl;
    CHECK(std::find(net.nodes.begin(), net.nodes.end(), "*1**1**1*****1") !=
          net.nodes.end());
    CHECK(std::find(net.nodes.begin(), net.nodes.end(), "**1*****1") !=
          net.nodes.end());
  }
  SECTION("generic bonds 2") {
    // this tests a very particular case where the generic bond scaffold is the
    // same as the generic scaffold that leads to it. Make sure we do not end up
    // with this kind of self edge
    auto m = "Cc1ccccc1OC1C(C)C1"_smiles;
    REQUIRE(m);
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.includeGenericBondScaffolds = true;
    ScaffoldNetwork::ScaffoldNetwork net;
    ScaffoldNetwork::detail::addMolToNetwork(*m, net, ps);
    CHECK(net.nodes.size() == 14);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 13);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Fragment;
                        }) == 2);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type ==
                                 ScaffoldNetwork::EdgeType::Initialize;
                        }) == 1);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Generic;
                        }) == 5);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type ==
                                 ScaffoldNetwork::EdgeType::GenericBond;
                        }) == 3);
  }
  SECTION("generic + no attach") {
    auto m = "Cc1ccccc1OC1C(C)C1"_smiles;
    REQUIRE(m);
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.includeGenericScaffolds = true;
    ps.includeGenericBondScaffolds = false;
    ps.includeScaffoldsWithoutAttachments = true;
    ps.keepOnlyFirstFragment = true;
    ScaffoldNetwork::ScaffoldNetwork net;
    ScaffoldNetwork::detail::addMolToNetwork(*m, net, ps);
    CHECK(net.nodes.size() == 11);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 10);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Fragment;
                        }) == 2);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type ==
                                 ScaffoldNetwork::EdgeType::Initialize;
                        }) == 1);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Generic;
                        }) == 5);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type ==
                                 ScaffoldNetwork::EdgeType::RemoveAttachment;
                        }) == 2);
    // std::copy(net.nodes.begin(), net.nodes.end(),
    //           std::ostream_iterator<std::string>(std::cerr, " "));
    // std::cerr << std::endl;
    CHECK(std::find(net.nodes.begin(), net.nodes.end(), "c1ccccc1") !=
          net.nodes.end());
    CHECK(std::find(net.nodes.begin(), net.nodes.end(), "*1:*:*:*:*:*:1") !=
          net.nodes.end());
  }
}

TEST_CASE("BRICS Fragmenter", "[unittest][scaffolds]") {
  auto m = "c1ccccc1C(=O)NC1NC(=O)CCC1"_smiles;
  REQUIRE(m);
  SECTION("original behavior default") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.includeScaffoldsWithoutAttachments = false;
    ps.includeGenericScaffolds = false;
    ps.keepOnlyFirstFragment = false;
    ScaffoldNetwork::ScaffoldNetwork net;
    ScaffoldNetwork::detail::addMolToNetwork(*m, net, ps);
    // std::copy(net.nodes.begin(), net.nodes.end(),
    //           std::ostream_iterator<std::string>(std::cerr, " "));
    // std::cerr << std::endl;
    CHECK(net.nodes.size() == 6);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 8);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Fragment;
                        }) == 8);
    CHECK(std::count(net.counts.begin(), net.counts.end(), 1) == 3);
    CHECK(std::count(net.counts.begin(), net.counts.end(), 2) == 3);
  }

  SECTION("BRICS fragmenter") {
    ScaffoldNetwork::ScaffoldNetworkParams ps =
        ScaffoldNetwork::getBRICSNetworkParams();
    ps.includeScaffoldsWithoutAttachments = false;
    ps.includeGenericScaffolds = false;
    ScaffoldNetwork::ScaffoldNetwork net;
    ScaffoldNetwork::detail::addMolToNetwork(*m, net, ps);
    CHECK(net.nodes.size() == 10);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 20);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Fragment;
                        }) == 20);
  }
}

TEST_CASE("Implicit Hs on aromatic atoms with attachments",
          "[bug][scaffolds]") {
  auto m = "c1cn(C3CCC3)nc1"_smiles;
  REQUIRE(m);
  SECTION("original behavior default") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.includeScaffoldsWithoutAttachments = true;
    ps.includeGenericScaffolds = false;
    ScaffoldNetwork::ScaffoldNetwork net;
    ScaffoldNetwork::detail::addMolToNetwork(*m, net, ps);
    CHECK(net.nodes.size() == 5);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 4);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Fragment;
                        }) == 2);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type ==
                                 ScaffoldNetwork::EdgeType::RemoveAttachment;
                        }) == 2);
    CHECK(std::count(net.counts.begin(), net.counts.end(), 1) == 5);
    for (auto nd : net.nodes) {
      std::unique_ptr<ROMol> m(SmilesToMol(nd));
      CHECK(m);
    }
  }
}

TEST_CASE("scaffold with attachment when attachments are disabled",
          "[bug][scaffolds]") {
  auto m = "C1CCC1C1CCCC1C1CCCCC1"_smiles;
  REQUIRE(m);
  SECTION("bug report") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.includeScaffoldsWithoutAttachments = true;
    ps.includeScaffoldsWithAttachments = false;
    ps.includeGenericScaffolds = false;
    ScaffoldNetwork::ScaffoldNetwork net;
    ScaffoldNetwork::detail::addMolToNetwork(*m, net, ps);
    CHECK(net.nodes.size() == 6);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 8);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Fragment;
                        }) == 8);

    CHECK(std::count(net.counts.begin(), net.counts.end(), 1) == 3);
    CHECK(std::count(net.counts.begin(), net.counts.end(), 2) == 3);
    for (auto nd : net.nodes) {
      CHECK(nd.find("*") == std::string::npos);
      std::unique_ptr<ROMol> m(SmilesToMol(nd));
      CHECK(m);
    }
  }
}

TEST_CASE("larger multi-mol test", "[regression][scaffold]") {
  std::vector<std::string> smiles{
      "Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]"
      "12",
      "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccccc3)C(=O)N2[C@H]1C(=O)O",
      "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O.[Na]",
      "Cc1onc(-c2ccccc2)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12"};
  std::vector<std::shared_ptr<ROMol>> ms;
  ms.reserve(smiles.size());
  std::transform(smiles.cbegin(), smiles.cend(), std::back_inserter(ms),
                 [](const std::string &smi) {
                   return std::shared_ptr<ROMol>(SmilesToMol(smi));
                 });

  SECTION("basics") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.includeGenericScaffolds = false;
    ps.includeScaffoldsWithoutAttachments = false;
    ScaffoldNetwork::ScaffoldNetwork net =
        ScaffoldNetwork::createScaffoldNetwork(ms, ps);
    CHECK(net.nodes.size() == 11);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 14);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Fragment;
                        }) == 10);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type ==
                                 ScaffoldNetwork::EdgeType::Initialize;
                        }) == 4);
    auto snodes = net.nodes;
    std::sort(snodes.begin(), snodes.end());
    std::vector<std::string> tgt{
        "*C1C(=O)N2CCSC12",
        "*c1ccccc1",
        "*c1conc1*",
        "*c1conc1-c1ccccc1",
        "*c1nocc1C(=O)NC1C(=O)N2CCSC12",
        "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O.[Na]",
        "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccccc3)C(=O)N2[C@H]1C(=O)O",
        "Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@@H]1C(=O)N2[C@@H]1SC(C)(C)[C@@H]2C(=O)O",
        "Cc1onc(-c2ccccc2)c1C(=O)N[C@@H]1C(=O)N2[C@@H]1SC(C)(C)[C@@H]2C(=O)O",
        "O=C(Cc1ccccc1)NC1C(=O)N2CCSC12",
        "O=C(NC1C(=O)N2CCSC12)c1conc1-c1ccccc1"};
    CHECK(snodes == tgt);
  }
  SECTION("generics") {
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.includeGenericScaffolds = true;
    ps.includeScaffoldsWithoutAttachments = false;
    ScaffoldNetwork::ScaffoldNetwork net =
        ScaffoldNetwork::createScaffoldNetwork(ms, ps);
    CHECK(net.nodes.size() == 18);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 21);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Fragment;
                        }) == 10);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type ==
                                 ScaffoldNetwork::EdgeType::Initialize;
                        }) == 4);
    CHECK(std::count_if(net.edges.begin(), net.edges.end(),
                        [](ScaffoldNetwork::NetworkEdge e) {
                          return e.type == ScaffoldNetwork::EdgeType::Generic;
                        }) == 7);
    auto snodes = net.nodes;
    std::sort(snodes.begin(), snodes.end());
    // std::copy(snodes.begin(), snodes.end(),
    //           std::ostream_iterator<std::string>(std::cerr, " "));
    // std::cerr << std::endl;
    std::vector<std::string> tgt{
        "**1*2****2*1=*",
        "**1:*:*:*:*:*:1",
        "**1:*:*:*:*:1*",
        "**1:*:*:*:*:1*(=*)**1*2****2*1=*",
        "**1:*:*:*:*:1*1:*:*:*:*:*:1",
        "*=*(**1*2****2*1=*)**1:*:*:*:*:*:1",
        "*=*1*2****2*1**(=*)*1:*:*:*:*:1*1:*:*:*:*:*:1",
        "*C1C(=O)N2CCSC12",
        "*c1ccccc1",
        "*c1conc1*",
        "*c1conc1-c1ccccc1",
        "*c1nocc1C(=O)NC1C(=O)N2CCSC12",
        "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O.[Na]",
        "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccccc3)C(=O)N2[C@H]1C(=O)O",
        "Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@@H]1C(=O)N2[C@@H]1SC(C)(C)[C@@H]2C(=O)O",
        "Cc1onc(-c2ccccc2)c1C(=O)N[C@@H]1C(=O)N2[C@@H]1SC(C)(C)[C@@H]2C(=O)O",
        "O=C(Cc1ccccc1)NC1C(=O)N2CCSC12",
        "O=C(NC1C(=O)N2CCSC12)c1conc1-c1ccccc1"};
    CHECK(snodes == tgt);
  }
}
#endif

TEST_CASE("BRICS performance problems", "[bug]") {
  SECTION("mol 1") {
    auto m =
        "COc1ccc(cc1)-c1cc(ccc1-c1ccc(OC)cc1)C(=O)OCCOC(=O)CCOC1CCCCC1"_smiles;
    REQUIRE(m);
    ScaffoldNetwork::ScaffoldNetworkParams ps =
        ScaffoldNetwork::getBRICSNetworkParams();
    ps.includeGenericScaffolds = false;
    ps.includeScaffoldsWithoutAttachments = false;
    ScaffoldNetwork::ScaffoldNetwork net;
    ScaffoldNetwork::detail::addMolToNetwork(*m, net, ps);
    CHECK(net.nodes.size() == 59);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 387);
  }
  SECTION("mol 2") {
    auto m =
        "COc1ccc(cc1)-c1cc(ccc1-c1ccc(OC)cc1)C(=O)OCCOC(=O)CCOCCOCCOCCOCCOC1CCCCC1"_smiles;
    REQUIRE(m);
    ScaffoldNetwork::ScaffoldNetworkParams ps =
        ScaffoldNetwork::getBRICSNetworkParams();
    ps.includeGenericScaffolds = false;
    ps.includeScaffoldsWithoutAttachments = false;
    ScaffoldNetwork::ScaffoldNetwork net;
    ScaffoldNetwork::detail::addMolToNetwork(*m, net, ps);
    CHECK(net.nodes.size() == 148);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 1969);
  }
}

TEST_CASE("Github #3177: seg fault with null molecules", "[bug]") {
  SECTION("basics") {
    std::vector<ROMOL_SPTR> mols;
    mols.emplace_back(nullptr);
    ScaffoldNetwork::ScaffoldNetworkParams ps =
        ScaffoldNetwork::getBRICSNetworkParams();
    REQUIRE_THROWS_AS(ScaffoldNetwork::createScaffoldNetwork(mols, ps),
                      ValueErrorException);
  }
  SECTION("including one valid molecule") {
    std::vector<ROMOL_SPTR> mols;
    mols.emplace_back(SmilesToMol("Cc1ccccc1"));
    mols.emplace_back(nullptr);
    ScaffoldNetwork::ScaffoldNetworkParams ps =
        ScaffoldNetwork::getBRICSNetworkParams();
    REQUIRE_THROWS_AS(ScaffoldNetwork::createScaffoldNetwork(mols, ps),
                      ValueErrorException);
  }
}

TEST_CASE("GitHub #3153: Kekulization error in molecules with aromatic C+",
          "[bug]") {
  SECTION("Standard Representation") {
    auto smis = {"O=C1C=CC(CC2=CC=CC2)=CC=C1"};
    std::vector<ROMOL_SPTR> ms;
    for (const auto smi : smis) {
      auto m = SmilesToMol(smi);
      REQUIRE(m);
      ms.push_back(ROMOL_SPTR(m));
    }
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ScaffoldNetwork::ScaffoldNetwork net =
        ScaffoldNetwork::createScaffoldNetwork(ms, ps);
    CHECK(net.nodes.size() == 9);
    CHECK(std::find(net.nodes.begin(), net.nodes.end(), "O=c1cccccc1") !=
          net.nodes.end());
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 8);
  }
  SECTION("Heteroatom inside ring") {
    auto smis = {"c1cccn1CC"};
    std::vector<ROMOL_SPTR> ms;
    for (const auto smi : smis) {
      auto m = SmilesToMol(smi);
      REQUIRE(m);
      ms.push_back(ROMOL_SPTR(m));
    }
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ScaffoldNetwork::ScaffoldNetwork net =
        ScaffoldNetwork::createScaffoldNetwork(ms, ps);

    CHECK(net.nodes.size() == 3);
    CHECK(std::find(net.nodes.begin(), net.nodes.end(), "c1cc[nH]c1") !=
          net.nodes.end());
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 2);
  }
  SECTION("Aromatic Carbocation") {
    auto smis = {"[O-][C+]1C=CC(CC2=CC=CC2)=CC=C1"};
    std::vector<ROMOL_SPTR> ms;
    for (const auto smi : smis) {
      auto m = SmilesToMol(smi);
      REQUIRE(m);
      ms.push_back(ROMOL_SPTR(m));
    }
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ScaffoldNetwork::ScaffoldNetwork net =
        ScaffoldNetwork::createScaffoldNetwork(ms, ps);
    CHECK(std::find(net.nodes.begin(), net.nodes.end(),
                    "C1=CCC(Cc2ccc[cH+]cc2)=C1") != net.nodes.end());
    CHECK(net.nodes.size() == 11);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 10);
  }
}

#ifdef RDK_USE_BOOST_SERIALIZATION

TEST_CASE("Serialization", "[serialization]") {
  auto smis = {"c1ccccc1CC1NC(=O)CCC1", "c1cccnc1CC1NC(=O)CCC1"};
  std::vector<ROMOL_SPTR> ms;
  for (const auto smi : smis) {
    auto m = SmilesToMol(smi);
    REQUIRE(m);
    ms.push_back(ROMOL_SPTR(m));
  }
  SECTION("basics") {
    // start by building the network
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ScaffoldNetwork::ScaffoldNetwork net;
    ScaffoldNetwork::updateScaffoldNetwork(ms, net, ps);
    CHECK(net.nodes.size() == 12);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 13);
    std::stringstream ss;
    boost::archive::text_oarchive oa(ss);
    oa << net;
    std::string pkl = ss.str();
    {
      std::stringstream iss(pkl);
      boost::archive::text_iarchive ia(iss);
      ScaffoldNetwork::ScaffoldNetwork net2;
      ia >> net2;
      CHECK(net2.nodes.size() == 12);
      CHECK(net2.counts.size() == net2.nodes.size());
      CHECK(net2.edges.size() == 13);
      CHECK(net2.nodes == net.nodes);
      CHECK(net2.counts == net.counts);
      CHECK(net2.edges == net.edges);
    }
    {
      std::stringstream iss(pkl);
      boost::archive::text_iarchive ia(iss);
      ScaffoldNetwork::ScaffoldNetwork net2(pkl);
      CHECK(net2.nodes.size() == 12);
      CHECK(net2.counts.size() == net2.nodes.size());
      CHECK(net2.edges.size() == 13);
      CHECK(net2.nodes == net.nodes);
      CHECK(net2.counts == net.counts);
      CHECK(net2.edges == net.edges);
    }
  }
}
#endif

TEST_CASE("molCounts", "[scaffolds]") {
  SECTION("basics") {
    auto smis = {"C1CC(C1)C1C(C1C1NCCCC1)C1OCCC1",
                 "C1CC(C1)C1C(C1C1NCCCC1)C1CCCC1"};
    std::vector<ROMOL_SPTR> ms;
    for (const auto smi : smis) {
      auto m = SmilesToMol(smi);
      REQUIRE(m);
      ms.push_back(ROMOL_SPTR(m));
    }
    ScaffoldNetwork::ScaffoldNetworkParams ps;
    ps.includeGenericScaffolds = false;
    ps.includeScaffoldsWithoutAttachments = false;
    ScaffoldNetwork::ScaffoldNetwork net =
        ScaffoldNetwork::createScaffoldNetwork(ms, ps);
    CHECK(net.nodes.size() == 16);
    CHECK(net.counts.size() == net.nodes.size());
    CHECK(net.edges.size() == 40);
    std::vector<std::pair<std::string, unsigned int>> endps = {
        {"*C1CCC1", 2u},
        {"*C1CCCCN1", 2u},
        {"*C1CCCO1", 1u},
        {"*C1CCCC1", 1u},
        {"*C1C(*)C1*", 2u}};
    for (const auto &endp : endps) {
      auto loc = std::find(net.nodes.begin(), net.nodes.end(), endp.first);
      CHECK(loc != net.nodes.end());
      auto idx = loc - net.nodes.begin();
      CHECK(net.molCounts[idx] == endp.second);
    }
  }
}