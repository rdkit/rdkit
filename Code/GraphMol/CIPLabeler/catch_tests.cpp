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
#include <list>
#include <string>
#include <vector>

#include "catch.hpp"

#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include "CIPLabeler.h"
#include "Digraph.h"
#include "rules/Pairlist.hpp"
#include "rules/Rule1a.h"
#include "rules/Rule2.h"

#include "CIPMol.h"

using namespace RDKit;
using namespace RDKit::CIPLabeler;

std::string toBinaryString(int value) {
  return std::bitset<32>(value).to_string();
}

TEST_CASE("Descriptor lists", "[accurateCIP]") {
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
}

void check_incoming_edge_count(Node *root) {
  auto queue = std::list<Node *>({root});

  for (auto itr = queue.begin(); itr != queue.end(); ++itr) {
    const auto &node = *itr;

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

void expandAll(Digraph &g) {
  auto queue = std::list<Node *>({g.getOriginRoot()});

  for (auto itr = queue.begin(); itr != queue.end(); ++itr) {
    const auto &node = *itr;

    for (const auto &e : node->getEdges()) {
      if (!e->isBeg(node)) {
        continue;
      }
      if (!e->getEnd()->isTerminal()) {
        queue.push_back(e->getEnd());
      }
    }
  }
}

TEST_CASE("Digraph", "[accurateCIP]") {
  auto mol =
      R"(CC1(OC2=C(C=3NC[C@@]4(C3C=C2)C([C@@H]5C[C@@]67C(N([C@]5(C4)CN6CC[C@@]7(C)O)C)=O)(C)C)OC=C1)C)"_smiles; // VS040
  auto cipmol = CIPLabeler::CIPMol(mol.get());

  auto initial_root_idx = 1;
  auto initial_root_atom = cipmol.getAtom(initial_root_idx);

  Digraph g(cipmol, initial_root_atom);
  expandAll(g);
  REQUIRE(g.getNumNodes() == 3819);

  auto current_root = g.getCurrentRoot();
  REQUIRE(current_root->getAtom()->getIdx() == initial_root_idx);

  check_incoming_edge_count(current_root);

  auto new_root_idx = 24;
  auto new_root_atom = cipmol.getAtom(new_root_idx);
  auto new_root_nodes = g.getNodes(new_root_atom);
  CHECK(new_root_nodes.size() == 104);

  Node *new_root_node = *new_root_nodes.begin();
  for (const auto &node : new_root_nodes) {
    if (!node->isDuplicate() &&
        node->getDistance() > new_root_node->getDistance()) {
      new_root_node = node;
    }
  }
  REQUIRE(new_root_node->getDistance() == 25);

  g.changeRoot(new_root_node);
  REQUIRE(g.getNumNodes() == 3819);

  current_root = g.getCurrentRoot();
  REQUIRE(current_root->getAtom()->getIdx() == new_root_idx);

  check_incoming_edge_count(current_root);
}

TEST_CASE("Rule1a", "[accurateCIP]") {
  SECTION("Compare equal") {
    auto mol = "COC"_smiles;
    auto cipmol = CIPLabeler::CIPMol(mol.get());
    Digraph g(cipmol, cipmol.getAtom(1));
    auto origin = g.getOriginRoot();

    auto frac = origin->getAtomicNumFraction();
    REQUIRE(frac.numerator() == 8);
    REQUIRE(frac.denominator() == 1);

    Rule1a rule;

    auto edges = origin->getEdges();
    REQUIRE(edges.size() == 2);
    REQUIRE(edges[0] != edges[1]);

    CHECK(0 == rule.compare(edges[0], edges[1]));

    CHECK(!rule.getSorter()->prioritise(origin, edges).isUnique());
  }

  SECTION("Compare different") {
    auto mol = "CON"_smiles;
    auto cipmol = CIPLabeler::CIPMol(mol.get());
    Digraph g(cipmol, cipmol.getAtom(1));
    auto origin = g.getOriginRoot();

    auto frac = origin->getAtomicNumFraction();
    REQUIRE(frac.numerator() == 8);
    REQUIRE(frac.denominator() == 1);

    Rule1a rule;

    auto edges = origin->getEdges();
    REQUIRE(edges.size() == 2);
    REQUIRE(edges[0] != edges[1]);

    auto frac0 = edges[0]->getEnd()->getAtomicNumFraction();
    REQUIRE(frac0.numerator() == 6);
    REQUIRE(frac0.denominator() == 1);

    auto frac1 = edges[1]->getEnd()->getAtomicNumFraction();
    REQUIRE(frac1.numerator() == 7);
    REQUIRE(frac1.denominator() == 1);

    CHECK(rule.compare(edges[0], edges[1]) < 0);
    CHECK(rule.compare(edges[1], edges[0]) > 0);

    CHECK(rule.getSorter()->prioritise(origin, edges).isUnique());
  }
}

TEST_CASE("Rule2", "[accurateCIP]") {
  SECTION("Compare equal") {
    auto mol = "COC"_smiles;
    auto cipmol = CIPLabeler::CIPMol(mol.get());
    Digraph g(cipmol, cipmol.getAtom(1));
    auto origin = g.getOriginRoot();

    auto frac = origin->getAtomicNumFraction();
    REQUIRE(frac.numerator() == 8);
    REQUIRE(frac.denominator() == 1);

    Rule2 rule;

    auto edges = origin->getEdges();
    REQUIRE(edges.size() == 2);
    REQUIRE(edges[0] != edges[1]);

    CHECK(0 == rule.compare(edges[0], edges[1]));

    CHECK(!rule.getSorter()->prioritise(origin, edges).isUnique());
  }

  SECTION("Compare different: Zero Mass Num") {
    auto mol = "CO[13C]"_smiles;
    auto cipmol = CIPLabeler::CIPMol(mol.get());
    Digraph g(cipmol, cipmol.getAtom(1));
    auto origin = g.getOriginRoot();

    auto frac = origin->getAtomicNumFraction();
    REQUIRE(frac.numerator() == 8);
    REQUIRE(frac.denominator() == 1);

    Rule2 rule;

    auto edges = origin->getEdges();
    REQUIRE(edges.size() == 2);
    REQUIRE(edges[0] != edges[1]);

    REQUIRE(edges[0]->getEnd()->getMassNum() == 0);
    REQUIRE(edges[1]->getEnd()->getMassNum() == 13);

    CHECK(rule.compare(edges[0], edges[1]) < 0);
    CHECK(rule.compare(edges[1], edges[0]) > 0);

    CHECK(rule.getSorter()->prioritise(origin, edges).isUnique());
  }

  SECTION("Compare different: Non-Zero Mass Num") {
    auto mol = "[13C]O[14C]"_smiles;
    auto cipmol = CIPLabeler::CIPMol(mol.get());
    Digraph g(cipmol, cipmol.getAtom(1));
    auto origin = g.getOriginRoot();

    auto frac = origin->getAtomicNumFraction();
    REQUIRE(frac.numerator() == 8);
    REQUIRE(frac.denominator() == 1);

    Rule2 rule;

    auto edges = origin->getEdges();
    REQUIRE(edges.size() == 2);
    REQUIRE(edges[0] != edges[1]);

    REQUIRE(edges[0]->getEnd()->getMassNum() == 13);
    REQUIRE(edges[1]->getEnd()->getMassNum() == 14);

    CHECK(rule.compare(edges[0], edges[1]) < 0);
    CHECK(rule.compare(edges[1], edges[0]) > 0);

    CHECK(rule.getSorter()->prioritise(origin, edges).isUnique());
  }
}

TEST_CASE("Tetrahedral assignment", "[accurateCIP]") {
  auto mol = "Br[C@H](Cl)F"_smiles;
  REQUIRE(mol->getNumAtoms() == 4);

  auto chiral_atom = mol->getAtomWithIdx(1);
  chiral_atom->clearProp(common_properties::_CIPCode);
  REQUIRE(chiral_atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);

  CIPLabeler::assignCIPLabels(*mol);

  std::string chirality;
  CHECK(chiral_atom->getPropIfPresent(common_properties::_CIPCode, chirality) ==
        true);
  CHECK(chirality == "S");
}

TEST_CASE("Double bond stereo assignment", "[accurateCIP]") {
  auto mol = R"(CC\C(\C(\C)=N\O)=N\O)"_smiles; // VS013
  REQUIRE(mol->getNumAtoms() == 9);

  auto bond_1 = mol->getBondWithIdx(4);
  auto bond_2 = mol->getBondWithIdx(6);
  REQUIRE(bond_1->getBondType() == Bond::DOUBLE);
  REQUIRE(bond_2->getBondType() == Bond::DOUBLE);
  REQUIRE(bond_1->getStereo() == Bond::STEREOE);
  REQUIRE(bond_2->getStereo() == Bond::STEREOZ);

  CIPLabeler::assignCIPLabels(*mol);

  std::string chirality;
  CHECK(bond_1->getPropIfPresent(common_properties::_CIPCode, chirality) ==
        true);
  CHECK(chirality == "E");

  CHECK(bond_2->getPropIfPresent(common_properties::_CIPCode, chirality) ==
        true);
  CHECK(chirality == "Z");
}
