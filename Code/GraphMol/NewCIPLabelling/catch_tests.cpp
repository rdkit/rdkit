//
//  Copyright (C) 2020 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do
                          // this in one cpp file

#include <bitset>
#include <string>
#include <vector>

#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include "catch.hpp"

#include "Digraph.hpp"
#include "rules/Pairlist.hpp"
#include "rules/Rule1a.hpp"
#include "rules/Rule2.hpp"

#include "RDKitCipMol.hpp"

using namespace RDKit;
using namespace RDKit::NewCIPLabelling;

std::string toBinaryString(int value) {
  return std::bitset<32>(value).to_string();
}

TEST_CASE("Descriptor lists", "[accurateCip]") {
  auto descriptors = PairList();

  SECTION("IgnoreConstructionNulls") {
    CHECK(!descriptors.add(Descriptor::NONE));
    CHECK(!descriptors.add(Descriptor::UNKNOWN));
    CHECK(!descriptors.add(Descriptor::ns));
  }
  SECTION("IgnoreConstructionPseudo") {
    CHECK(!descriptors.add(Descriptor::r));
    CHECK(!descriptors.add(Descriptor::s));

    CHECK(descriptors.add(Descriptor::R));
    CHECK(descriptors.add(Descriptor::S));
  }
  SECTION("Pairing") {
    REQUIRE(descriptors.getPairing() == 0);

    CHECK("00000000000000000000000000000000" ==
          toBinaryString(descriptors.getPairing()));

    descriptors.add(Descriptor::R);
    CHECK("00000000000000000000000000000000" ==
          toBinaryString(descriptors.getPairing()));

    // like
    descriptors.add(Descriptor::R);
    CHECK("01000000000000000000000000000000" ==
          toBinaryString(descriptors.getPairing()));

    // like
    descriptors.add(Descriptor::R);
    CHECK("01100000000000000000000000000000" ==
          toBinaryString(descriptors.getPairing()));

    // unlike
    descriptors.add(Descriptor::S);
    CHECK("01100000000000000000000000000000" ==
          toBinaryString(descriptors.getPairing()));

    // like
    descriptors.add(Descriptor::R);
    CHECK("01101000000000000000000000000000" ==
          toBinaryString(descriptors.getPairing()));

    // like
    descriptors.add(Descriptor::R);
    CHECK("01101100000000000000000000000000" ==
          toBinaryString(descriptors.getPairing()));

    // like
    descriptors.add(Descriptor::R);
    CHECK("01101110000000000000000000000000" ==
          toBinaryString(descriptors.getPairing()));

    // unlike
    descriptors.add(Descriptor::S);
    CHECK("01101110000000000000000000000000" ==
          toBinaryString(descriptors.getPairing()));

    // like
    descriptors.add(Descriptor::R);
    CHECK("01101110100000000000000000000000" ==
          toBinaryString(descriptors.getPairing()));
  }
  SECTION("Append Empty") {
    descriptors.add(Descriptor::R);
    descriptors.add(Descriptor::R);
    descriptors.add(Descriptor::S);
    descriptors.add(Descriptor::R);

    const auto tmp = std::vector<PairList>(1);
    const auto lists = descriptors.append(tmp);
    REQUIRE(1u == lists.size());
    CHECK(descriptors.getPairing() == lists[0].getPairing());
  }
  SECTION("Append") {
    descriptors.add(Descriptor::R);
    descriptors.add(Descriptor::R);
    descriptors.add(Descriptor::S);
    descriptors.add(Descriptor::R);

    PairList tail1 = PairList();
    tail1.add(Descriptor::R);
    tail1.add(Descriptor::S);
    tail1.add(Descriptor::R);

    PairList tail2 = PairList();
    tail2.add(Descriptor::S);
    tail2.add(Descriptor::S);
    tail2.add(Descriptor::R);

    const auto tmp = std::vector<PairList>({tail1, tail2});
    const auto created = descriptors.append(tmp);

    REQUIRE(2 == created.size());

    CHECK("01011010000000000000000000000000" ==
          toBinaryString(created[0].getPairing()));
    CHECK("01010010000000000000000000000000" ==
          toBinaryString(created[1].getPairing()));
  }
  SECTION("pairRM") {
    PairList list1 = PairList();
    PairList list2 = PairList();
    list1.add(Descriptor::R);
    list1.add(Descriptor::M);
    list1.add(Descriptor::R);
    list1.add(Descriptor::S);
    list2.add(Descriptor::R);
    list2.add(Descriptor::P);
    list2.add(Descriptor::S);
    list2.add(Descriptor::M);
    CHECK(list1.toString() == "R:llu");
    CHECK(list2.toString() == "R:uul");
  }
  SECTION("Clear") {
    PairList descriptors = PairList();
    descriptors.add(Descriptor::R);
    descriptors.add(Descriptor::R);
    descriptors.add(Descriptor::S);
    descriptors.add(Descriptor::R);
    descriptors.add(Descriptor::S);

    CHECK(descriptors.getPairing() > 0);
    descriptors.clear();
    CHECK(descriptors.getPairing() == 0);
  }
}

void check_incoming_edge_count(Node<RdkA, RdkB> *root) {
  auto queue = std::vector<Node<RdkA, RdkB> *>({root});

  for (auto pos = 0u; pos < queue.size(); ++pos) {
    const auto node = queue[pos];

    int incoming_edges = 0;
    for (const auto &e : node->getEdges()) {
      if (!e->isBeg(node)) {
        REQUIRE(e->isEnd(node));
        ++incoming_edges;
      } else if (!e->getEnd()->isTerminal()) {
        queue.push_back(e->getEnd());
      }
    }

    // Only the current root should have no incoming nodes.
    // All other nodes should have exactly one incoming edge.
    if (node == root) {
      CHECK(incoming_edges == 0);
    } else {
      CHECK(incoming_edges == 1);
    }
  }
}

TEST_CASE("Digraph", "[accurateCip]") {
  auto mol = R"(CC\C(\C(\C)=N\O)=N\O)"_smiles; // VS013
  auto cipmol = NewCIPLabelling::RDKitCipMol(mol.get());

  auto first_root_idx = 2;
  auto first_root_atom = cipmol.getAtom(first_root_idx);

  auto g = Digraph<RdkA, RdkB>(&cipmol, first_root_atom);
  g.expandAll();
  REQUIRE(g.getNumNodes() == 21);

  auto current_root = g.getCurrRoot();
  REQUIRE(current_root->getAtom()->getIdx() == first_root_idx);

  check_incoming_edge_count(current_root);

  auto second_root_idx = 3;
  auto second_root_atom = cipmol.getAtom(second_root_idx);
  auto second_root_nodes = g.getNodes(second_root_atom);
  CHECK(second_root_nodes.size() == 2);

  g.changeRoot(second_root_nodes[0]);
  REQUIRE(g.getNumNodes() == 21);

  current_root = g.getCurrRoot();
  REQUIRE(current_root->getAtom()->getIdx() == second_root_idx);

  check_incoming_edge_count(current_root);
}

TEST_CASE("Rule1a", "[accurateCip]") {
  SECTION("Compare equal") {
    auto mol = "COC"_smiles;
    auto cipmol = NewCIPLabelling::RDKitCipMol(mol.get());
    auto g = Digraph<RdkA, RdkB>(&cipmol, cipmol.getAtom(1));
    auto origin = g.getRoot();
    REQUIRE(origin->getAtomicNumNumerator() == 8);

    Rule1a<RdkA, RdkB> rule(&cipmol);

    auto edges = origin->getEdges();
    REQUIRE(edges.size() == 2);
    REQUIRE(edges[0] != edges[1]);

    CHECK(0 == rule.compare(edges[0], edges[1]));

    CHECK(!rule.getSorter()->prioritise(origin, edges).isUnique());
  }

  SECTION("Compare different") {
    auto mol = "CON"_smiles;
    auto cipmol = NewCIPLabelling::RDKitCipMol(mol.get());
    auto g = Digraph<RdkA, RdkB>(&cipmol, cipmol.getAtom(1));
    auto origin = g.getRoot();
    REQUIRE(origin->getAtomicNumNumerator() == 8);

    Rule1a<RdkA, RdkB> rule(&cipmol);

    auto edges = origin->getEdges();
    REQUIRE(edges.size() == 2);
    REQUIRE(edges[0] != edges[1]);

    REQUIRE(edges[0]->getEnd()->getAtomicNumNumerator() == 6);
    REQUIRE(edges[1]->getEnd()->getAtomicNumNumerator() == 7);

    CHECK(rule.compare(edges[0], edges[1]) < 0);
    CHECK(rule.compare(edges[1], edges[0]) > 0);

    CHECK(rule.getSorter()->prioritise(origin, edges).isUnique());
  }
}

TEST_CASE("Rule2", "[accurateCip]") {
  SECTION("Compare equal") {
    auto mol = "COC"_smiles;
    auto cipmol = NewCIPLabelling::RDKitCipMol(mol.get());
    auto g = Digraph<RdkA, RdkB>(&cipmol, cipmol.getAtom(1));
    auto origin = g.getRoot();
    REQUIRE(origin->getAtomicNumNumerator() == 8);

    Rule2<RdkA, RdkB> rule(&cipmol);

    auto edges = origin->getEdges();
    REQUIRE(edges.size() == 2);
    REQUIRE(edges[0] != edges[1]);

    CHECK(0 == rule.compare(edges[0], edges[1]));

    CHECK(!rule.getSorter()->prioritise(origin, edges).isUnique());
  }

  SECTION("Compare different: Zero Mass Num") {
    auto mol = "CO[13C]"_smiles;
    auto cipmol = NewCIPLabelling::RDKitCipMol(mol.get());
    auto g = Digraph<RdkA, RdkB>(&cipmol, cipmol.getAtom(1));
    auto origin = g.getRoot();
    REQUIRE(origin->getAtomicNumNumerator() == 8);

    Rule2<RdkA, RdkB> rule(&cipmol);

    auto edges = origin->getEdges();
    REQUIRE(edges.size() == 2);
    REQUIRE(edges[0] != edges[1]);

    REQUIRE(cipmol.getMassNum(edges[0]->getEnd()->getAtom()) == 0);
    REQUIRE(cipmol.getMassNum(edges[1]->getEnd()->getAtom()) == 13);

    CHECK(rule.compare(edges[0], edges[1]) < 0);
    CHECK(rule.compare(edges[1], edges[0]) > 0);

    CHECK(rule.getSorter()->prioritise(origin, edges).isUnique());
  }

  SECTION("Compare different: Non-Zero Mass Num") {
    auto mol = "[13C]O[14C]"_smiles;
    auto cipmol = NewCIPLabelling::RDKitCipMol(mol.get());
    auto g = Digraph<RdkA, RdkB>(&cipmol, cipmol.getAtom(1));
    auto origin = g.getRoot();
    REQUIRE(origin->getAtomicNumNumerator() == 8);

    Rule2<RdkA, RdkB> rule(&cipmol);

    auto edges = origin->getEdges();
    REQUIRE(edges.size() == 2);
    REQUIRE(edges[0] != edges[1]);

    REQUIRE(cipmol.getMassNum(edges[0]->getEnd()->getAtom()) == 13);
    REQUIRE(cipmol.getMassNum(edges[1]->getEnd()->getAtom()) == 14);

    CHECK(rule.compare(edges[0], edges[1]) < 0);
    CHECK(rule.compare(edges[1], edges[0]) > 0);

    CHECK(rule.getSorter()->prioritise(origin, edges).isUnique());
  }
}

TEST_CASE("Tetrahedral assignment", "[accurateCip]") {
  auto mol = "Br[C@H](Cl)F"_smiles;
  REQUIRE(mol->getNumAtoms() == 4);

  auto chiral_atom = mol->getAtomWithIdx(1);
  chiral_atom->clearProp(common_properties::_CIPCode);
  REQUIRE(chiral_atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);

  MolOps::assignRealCIPStereo(*mol);

  std::string chirality;
  CHECK(chiral_atom->getPropIfPresent(common_properties::_CIPCode, chirality) ==
        true);
  CHECK(chirality == "S");
}

TEST_CASE("Double bond stereo assignment", "[accurateCip]") {
  auto mol = R"(CC\C(\C(\C)=N\O)=N\O)"_smiles; // VS013
  REQUIRE(mol->getNumAtoms() == 9);

  auto bond_1 = mol->getBondWithIdx(4);
  auto bond_2 = mol->getBondWithIdx(6);
  REQUIRE(bond_1->getBondType() == Bond::DOUBLE);
  REQUIRE(bond_2->getBondType() == Bond::DOUBLE);
  REQUIRE(bond_1->getStereo() == Bond::STEREOE);
  REQUIRE(bond_2->getStereo() == Bond::STEREOZ);

  MolOps::assignRealCIPStereo(*mol);

  std::string chirality;
  CHECK(bond_1->getPropIfPresent(common_properties::_CIPCode, chirality) ==
        true);
  CHECK(chirality == "E");

  CHECK(bond_2->getPropIfPresent(common_properties::_CIPCode, chirality) ==
        true);
  CHECK(chirality == "Z");
}

TEST_CASE("RDKit Issues", "[accurateCip]") {
  SECTION("Issue #2984") {
    // RDKit doesn't see this stereo labels, so calculating stereo
    // will wipe this chirality and bond directions
    auto debugParse = false;
    auto sanitize = false;
    auto mol = SmilesToMol(R"(CC/C=C1\C[C@H](O)C1)", debugParse, sanitize);
    REQUIRE(mol->getNumAtoms() == 8);
    MolOps::sanitizeMol(*mol);

    auto bond = mol->getBondWithIdx(2);
    REQUIRE(bond->getBondType() == Bond::DOUBLE);
    REQUIRE(bond->getStereo() == Bond::STEREONONE);
    REQUIRE(mol->getBondWithIdx(1)->getBondDir() == Bond::ENDUPRIGHT);
    REQUIRE(mol->getBondWithIdx(3)->getBondDir() == Bond::ENDDOWNRIGHT);

    auto atom = mol->getAtomWithIdx(5);
    REQUIRE(atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);

    MolOps::assignRealCIPStereo(*mol);

    std::string chirality;
    CHECK(bond->getPropIfPresent(common_properties::_CIPCode, chirality) ==
          true);
    CHECK(chirality == "z");

    CHECK(atom->getPropIfPresent(common_properties::_CIPCode, chirality) ==
          true);
    CHECK(chirality == "R");
  }
}