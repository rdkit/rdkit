//
//  Copyright (C) 2020 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <bitset>
#include <list>
#include <string>
#include <vector>

#include <strstream>

#include <catch2/catch_all.hpp>

#include <GraphMol/MolOps.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MarvinParse/MarvinParser.h>
#include <GraphMol/test_fixtures.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "CIPLabeler.h"
#include "Digraph.h"
#include "rules/Pairlist.h"
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

  for (const auto &node : queue) {
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
  auto queue = std::list<Node *>({g.getOriginalRoot()});

  for (const auto &node : queue) {
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
      R"(CC1(OC2=C(C=3NC[C@@]4(C3C=C2)C([C@@H]5C[C@@]67C(N([C@]5(C4)CN6CC[C@@]7(C)O)C)=O)(C)C)OC=C1)C)"_smiles;  // VS040
  CIPLabeler::CIPMol cipmol(*mol);

  auto initial_root_idx = 1u;
  auto initial_root_atom = cipmol.getAtom(initial_root_idx);

  Digraph g(cipmol, initial_root_atom);
  expandAll(g);
  REQUIRE(g.getNumNodes() == 3819);

  auto current_root = g.getCurrentRoot();
  REQUIRE(current_root->getAtom()->getIdx() == initial_root_idx);

  check_incoming_edge_count(current_root);

  auto new_root_idx = 24u;
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
    CIPLabeler::CIPMol cipmol(*mol);
    Digraph g(cipmol, cipmol.getAtom(1));
    auto origin = g.getOriginalRoot();

    auto frac = origin->getAtomicNumFraction();
    REQUIRE(frac.numerator() == 8);
    REQUIRE(frac.denominator() == 1);

    Rule1a rule;

    auto edges = origin->getEdges();
    REQUIRE(edges.size() == 2);
    REQUIRE(edges[0] != edges[1]);

    CHECK(0 == rule.compare(edges[0], edges[1]));

    CHECK(!rule.getSorter()->prioritize(origin, edges).isUnique());
  }

  SECTION("Compare different") {
    auto mol = "CON"_smiles;
    CIPLabeler::CIPMol cipmol(*mol);
    Digraph g(cipmol, cipmol.getAtom(1));
    auto origin = g.getOriginalRoot();

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

    CHECK(rule.getSorter()->prioritize(origin, edges).isUnique());
  }
}

TEST_CASE("Rule2", "[accurateCIP]") {
  SECTION("Compare equal") {
    auto mol = "COC"_smiles;
    CIPLabeler::CIPMol cipmol(*mol);
    Digraph g(cipmol, cipmol.getAtom(1));
    auto origin = g.getOriginalRoot();

    auto frac = origin->getAtomicNumFraction();
    REQUIRE(frac.numerator() == 8);
    REQUIRE(frac.denominator() == 1);

    Rule2 rule;

    auto edges = origin->getEdges();
    REQUIRE(edges.size() == 2);
    REQUIRE(edges[0] != edges[1]);

    CHECK(0 == rule.compare(edges[0], edges[1]));

    CHECK(!rule.getSorter()->prioritize(origin, edges).isUnique());
  }

  SECTION("Compare different: Zero Mass Num") {
    auto mol = "CO[13C]"_smiles;
    CIPLabeler::CIPMol cipmol(*mol);
    Digraph g(cipmol, cipmol.getAtom(1));
    auto origin = g.getOriginalRoot();

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

    CHECK(rule.getSorter()->prioritize(origin, edges).isUnique());
  }

  SECTION("Compare different: Non-Zero Mass Num") {
    auto mol = "[13C]O[14C]"_smiles;
    CIPLabeler::CIPMol cipmol(*mol);
    Digraph g(cipmol, cipmol.getAtom(1));
    auto origin = g.getOriginalRoot();

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

    CHECK(rule.getSorter()->prioritize(origin, edges).isUnique());
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
  CHECK(chiral_atom->getPropIfPresent(common_properties::_CIPCode, chirality));
  CHECK(chirality == "S");
}

TEST_CASE("Double bond stereo assignment", "[accurateCIP]") {
  auto mol = R"(CC\C(\C(\C)=N\O)=N\O)"_smiles;  // VS013
  REQUIRE(mol->getNumAtoms() == 9);

  auto bond_1 = mol->getBondWithIdx(4);
  auto bond_2 = mol->getBondWithIdx(6);
  REQUIRE(bond_1->getBondType() == Bond::DOUBLE);
  REQUIRE(bond_2->getBondType() == Bond::DOUBLE);
  REQUIRE(bond_1->getStereo() == Bond::STEREOE);
  REQUIRE(bond_2->getStereo() == Bond::STEREOZ);

  CIPLabeler::assignCIPLabels(*mol);

  std::string chirality;
  CHECK(bond_1->getPropIfPresent(common_properties::_CIPCode, chirality));
  CHECK(chirality == "E");

  CHECK(bond_2->getPropIfPresent(common_properties::_CIPCode, chirality));
  CHECK(chirality == "Z");
}

TEST_CASE("phosphine and arsine chirality", "[accurateCIP]") {
  const std::vector<std::pair<std::string, std::string>> mols{
      {"C[P@](C1CCCC1)C1=CC=CC=C1", "R"},
      {"C[As@@](C1CCCC1)C1=CC=CC=C1", "S"},
      {"C[P@H]C1CCCCC1", "R"},
      {"C[P@@H]C1CCCCC1", "S"}};

  for (const auto &ref : mols) {
    std::unique_ptr<RWMol> mol{SmilesToMol(ref.first)};
    REQUIRE(mol);
    CIPLabeler::assignCIPLabels(*mol);

    std::string chirality;
    CHECK(mol->getAtomWithIdx(1)->getPropIfPresent(common_properties::_CIPCode,
                                                   chirality));
    CHECK(chirality == ref.second);
  }
}

TEST_CASE("assign specific atoms and bonds", "[accurateCIP]") {
  SECTION("Assign atoms") {
    auto mol = "C[C@H](Cl)CC[C@H](Cl)C"_smiles;
    REQUIRE(mol);

    auto atom1 = mol->getAtomWithIdx(1);
    auto atom5 = mol->getAtomWithIdx(5);

    REQUIRE(atom1->hasProp(common_properties::_CIPCode));
    REQUIRE(atom5->hasProp(common_properties::_CIPCode));

    atom1->clearProp(common_properties::_CIPCode);
    atom5->clearProp(common_properties::_CIPCode);

    boost::dynamic_bitset<> atoms(mol->getNumAtoms());
    boost::dynamic_bitset<> bonds;
    atoms.set(1);
    CIPLabeler::assignCIPLabels(*mol, atoms, bonds);

    std::string chirality;
    CHECK(atom1->getPropIfPresent(common_properties::_CIPCode, chirality));
    CHECK(chirality == "S");
    CHECK(!atom5->hasProp(common_properties::_CIPCode));
  }
  SECTION("Assign bonds") {
    auto mol = R"(C\C=C\C=C/C)"_smiles;
    REQUIRE(mol);

    auto bond1 = mol->getBondWithIdx(1);
    auto bond3 = mol->getBondWithIdx(3);

    REQUIRE(bond1->getBondType() == Bond::DOUBLE);
    REQUIRE(bond3->getBondType() == Bond::DOUBLE);

    REQUIRE(!bond1->hasProp(common_properties::_CIPCode));
    REQUIRE(!bond3->hasProp(common_properties::_CIPCode));

    boost::dynamic_bitset<> atoms;
    boost::dynamic_bitset<> bonds(mol->getNumBonds());
    bonds.set(3);
    CIPLabeler::assignCIPLabels(*mol, atoms, bonds);

    std::string stereo;
    CHECK(!bond1->hasProp(common_properties::_CIPCode));
    CHECK(bond3->getPropIfPresent(common_properties::_CIPCode, stereo));
    CHECK(stereo == "Z");
  }
}

TEST_CASE("para-stereochemistry", "[accurateCIP]") {
  SECTION("example 1") {
    // slightly simplified complex example from Salome Rieder
    auto mol = "C\\C=C/[C@@H](\\C=C\\O)[C@H](C)[C@H](\\C=C/C)\\C=C\\O"_smiles;
    REQUIRE(mol);
    CIPLabeler::assignCIPLabels(*mol);

    std::string chirality;
    CHECK(mol->getAtomWithIdx(3)->getPropIfPresent(common_properties::_CIPCode,
                                                   chirality));
    CHECK(chirality == "R");
    CHECK(mol->getAtomWithIdx(7)->getPropIfPresent(common_properties::_CIPCode,
                                                   chirality));
    CHECK(chirality == "r");
    CHECK(mol->getAtomWithIdx(9)->getPropIfPresent(common_properties::_CIPCode,
                                                   chirality));
    CHECK(chirality == "S");
  }
  SECTION("example 2") {
    // lovely complex example from Salome Rieder
    auto mol = "C\\C=C/[C@@H](\\C=C\\C)[C@H](C)[C@H](\\C=C/C)\\C=C\\C"_smiles;
    REQUIRE(mol);
    CIPLabeler::assignCIPLabels(*mol);

    std::string chirality;
    CHECK(mol->getAtomWithIdx(3)->getPropIfPresent(common_properties::_CIPCode,
                                                   chirality));
    CHECK(chirality == "R");
    CHECK(mol->getAtomWithIdx(7)->getPropIfPresent(common_properties::_CIPCode,
                                                   chirality));
    CHECK(chirality == "r");
    CHECK(mol->getAtomWithIdx(9)->getPropIfPresent(common_properties::_CIPCode,
                                                   chirality));
    CHECK(chirality == "S");
  }
}
TEST_CASE("para-stereochemistry2", "[accurateCIP]") {
  SECTION("example 1") {
    // example with one pseudo symmetry
    UseLegacyStereoPerceptionFixture useLegacy(false);

    auto mol = "OC(=O)[C@H]1CC[C@@H](CC1)O[C@@H](F)Cl"_smiles;
    REQUIRE(mol);
    CIPLabeler::assignCIPLabels(*mol);

    std::string chirality;
    CHECK(mol->getAtomWithIdx(3)->getPropIfPresent(common_properties::_CIPCode,
                                                   chirality));
    CHECK(chirality == "S");
    CHECK(mol->getAtomWithIdx(6)->getPropIfPresent(common_properties::_CIPCode,
                                                   chirality));
    CHECK(chirality == "r");
    CHECK(mol->getAtomWithIdx(10)->getPropIfPresent(common_properties::_CIPCode,
                                                    chirality));
    CHECK(chirality == "S");
  }

  SECTION("example 2") {
    // example with one pseudo symmetry
    UseLegacyStereoPerceptionFixture useLegacy(false);

    auto mol = "OC(=O)[C@H]1CC[C@@H](CC1)O[CH](Cl)Cl"_smiles;
    REQUIRE(mol);
    CIPLabeler::assignCIPLabels(*mol);

    std::string chirality;
    CHECK(mol->getAtomWithIdx(3)->getPropIfPresent(common_properties::_CIPCode,
                                                   chirality));
    CHECK(chirality == "r");
    CHECK(mol->getAtomWithIdx(6)->getPropIfPresent(common_properties::_CIPCode,
                                                   chirality));
    CHECK(chirality == "r");
    CHECK(!mol->getAtomWithIdx(10)->getPropIfPresent(
        common_properties::_CIPCode, chirality));
    // CHECK(chirality == "S");
  }
}

TEST_CASE(
    "Github #4996: Bad handling of dummy atoms in the CIP assignment code",
    "[accurateCIP]") {
  SECTION("case 1") {
    auto m = "*[C@](F)(Cl)Br"_smiles;
    REQUIRE(m);
    bool cleanit = true;
    bool force = true;
    // original assignment:
    MolOps::assignStereochemistry(*m, cleanit, force);
    std::string cip;
    CHECK(m->getAtomWithIdx(1)->getPropIfPresent(common_properties::_CIPCode,
                                                 cip));
    CHECK(cip == "S");

    m->getAtomWithIdx(1)->clearProp(common_properties::_CIPCode);
    CIPLabeler::assignCIPLabels(*m);
    CHECK(m->getAtomWithIdx(1)->getPropIfPresent(common_properties::_CIPCode,
                                                 cip));
    CHECK(cip == "S");
  }
  SECTION("dummies can match dummies") {
    auto m = "*[C@](*)(Cl)Br"_smiles;
    REQUIRE(m);
    bool cleanit = true;
    bool force = true;
    // original assignment:
    MolOps::assignStereochemistry(*m, cleanit, force);
    CHECK(!m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));

    CIPLabeler::assignCIPLabels(*m);
    CHECK(!m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  }
  SECTION("case 2") {
    auto m = "C1CC[C@](*)2CCCC[C@H]2C1"_smiles;
    REQUIRE(m);

    bool cleanit = true;
    bool force = true;
    // original assignment doesn't work for these:
    MolOps::assignStereochemistry(*m, cleanit, force);
    CHECK(!m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
    CHECK(!m->getAtomWithIdx(9)->hasProp(common_properties::_CIPCode));

    CIPLabeler::assignCIPLabels(*m);
    std::string cip;
    CHECK(m->getAtomWithIdx(3)->getPropIfPresent(common_properties::_CIPCode,
                                                 cip));
    CHECK(cip == "s");
    cip = "";
    CHECK(m->getAtomWithIdx(9)->getPropIfPresent(common_properties::_CIPCode,
                                                 cip));
    CHECK(cip == "s");
  }
}

TEST_CASE("CIP code errors on fragments which cannot be kekulized",
          "[accurateCIP]") {
  SECTION("fragment not affecting the stereochem") {
    auto m = "F[C@H](CNC)CCc(:c):c"_smarts;
    m->getAtomWithIdx(1)->clearProp(common_properties::_CIPCode);
    m->updatePropertyCache();
    CIPLabeler::assignCIPLabels(*m);
    std::string cip;
    CHECK(m->getAtomWithIdx(1)->getPropIfPresent(common_properties::_CIPCode,
                                                 cip));
    CHECK(cip == "S");
  }
  SECTION("fragment, unique bits") {
    auto m = "F[C@H](C(N)C)c(:c):c"_smarts;
    m->getAtomWithIdx(1)->clearProp(common_properties::_CIPCode);
    m->updatePropertyCache();
    CIPLabeler::assignCIPLabels(*m);
    std::string cip;
    CHECK(m->getAtomWithIdx(1)->getPropIfPresent(common_properties::_CIPCode,
                                                 cip));
  }
  SECTION("fragment, non-unique bits") {
    auto m = "F[C@H]([C]([NH])[CH2])c(:n):c"_smarts;
    m->getAtomWithIdx(1)->clearProp(common_properties::_CIPCode);
    m->updatePropertyCache();
    CIPLabeler::assignCIPLabels(*m);
    std::string cip;
    CHECK(!m->getAtomWithIdx(1)->getPropIfPresent(common_properties::_CIPCode,
                                                  cip));
  }
}

TEST_CASE("GitHub Issue #5142", "[bug][accurateCIP]") {
  auto mol = "*C1C[C@H](CCC)[C@@H](C)[C@H](C)C1"_smiles;
  REQUIRE(mol);
  CIPLabeler::assignCIPLabels(*mol);
}

TEST_CASE("CIP max iterations test", "[accurateCIP]") {
  SECTION("Exceeds max iterations") {
    std::string molBlock = R"(
  Mrv2117 11112217353D          

 40 50  0  0  0  0            999 V2000
    7.5483   -7.7451   -3.3419 H   0  0  0  0  0  0  0  0  0  0  0  0
    9.3984   -5.8792   -3.4318 H   0  0  0  0  0  0  0  0  0  0  0  0
    9.3399   -6.5425   -1.9898 C   0  0  2  0  0  0  0  0  0  0  0  0
    8.6820   -3.9425   -1.6474 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.8163   -7.1034   -1.6635 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.5722   -4.7343   -0.6304 H   0  0  0  0  0  0  0  0  0  0  0  0
    7.1701   -7.2380   -0.8890 C   0  0  1  0  0  0  0  0  0  0  0  0
   11.8300   -6.9765   -2.5544 H   0  0  0  0  0  0  0  0  0  0  0  0
   10.6631   -7.1690   -1.5375 C   0  0  1  0  0  0  0  0  0  0  0  0
   11.4077   -9.5239   -2.1472 H   0  0  0  0  0  0  0  0  0  0  0  0
   11.9006  -10.0139    0.4722 H   0  0  0  0  0  0  0  0  0  0  0  0
   10.6600   -8.8749    0.2914 C   0  0  2  0  0  0  0  0  0  0  0  0
   12.5136   -5.9115   -0.4190 H   0  0  0  0  0  0  0  0  0  0  0  0
   12.5338   -7.8406    1.4379 H   0  0  0  0  0  0  0  0  0  0  0  0
    8.6981   -8.8654    3.4177 H   0  0  0  0  0  0  0  0  0  0  0  0
    9.4046  -10.6907    1.7143 H   0  0  0  0  0  0  0  0  0  0  0  0
    9.3185   -9.3106    0.9700 C   0  0  1  0  0  0  0  0  0  0  0  0
    7.4958  -10.8075   -0.1783 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.5552   -7.8780    2.6983 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.8522   -9.1147    0.4514 H   0  0  0  0  0  0  0  0  0  0  0  0
    7.1601   -8.2997    0.2353 C   0  0  2  0  0  0  0  0  0  0  0  0
    7.5774   -7.5989    1.5700 C   0  0  1  0  0  0  0  0  0  0  0  0
   10.6577   -4.0099    0.1683 H   0  0  0  0  0  0  0  0  0  0  0  0
   11.0808   -6.4894   -0.1982 C   0  0  2  0  0  0  0  0  0  0  0  0
    8.9123   -5.4699   -0.9424 C   0  0  1  0  0  0  0  0  0  0  0  0
    8.9155   -8.2329    2.0184 C   0  0  2  0  0  0  0  0  0  0  0  0
   11.0766   -7.5513    0.9147 C   0  0  1  0  0  0  0  0  0  0  0  0
   10.6780   -6.9857    3.3021 H   0  0  0  0  0  0  0  0  0  0  0  0
    7.5861   -5.8939   -0.2667 C   0  0  2  0  0  0  0  0  0  0  0  0
    6.9973   -5.0597    1.9875 H   0  0  0  0  0  0  0  0  0  0  0  0
    9.3193   -5.8429    1.5048 C   0  0  2  0  0  0  0  0  0  0  0  0
    9.4740   -4.7173    2.5743 H   0  0  0  0  0  0  0  0  0  0  0  0
    9.9809   -7.1457    1.9811 C   0  0  1  0  0  0  0  0  0  0  0  0
    7.8270   -6.1245    1.2483 C   0  0  2  0  0  0  0  0  0  0  0  0
   10.0023   -5.4308    0.1366 C   0  0  1  0  0  0  0  0  0  0  0  0
    8.9125   -8.9385   -1.4737 C   0  0  2  0  0  0  0  0  0  0  0  0
    8.2374   -9.3522   -0.1173 C   0  0  1  0  0  0  0  0  0  0  0  0
   10.4116   -8.6461   -1.2148 C   0  0  1  0  0  0  0  0  0  0  0  0
    8.7726   -9.9758   -2.5047 H   0  0  0  0  0  0  0  0  0  0  0  0
    8.2577   -7.6250   -1.9478 C   0  0  2  0  0  0  0  0  0  0  0  0
  1 40  1  0  0  0  0
  2  3  1  0  0  0  0
  3  9  1  0  0  0  0
  3 25  1  0  0  0  0
  3 40  1  0  0  0  0
  4 25  1  0  0  0  0
  5  7  1  0  0  0  0
  6 29  1  0  0  0  0
  7 21  1  0  0  0  0
  7 29  1  0  0  0  0
  7 40  1  0  0  0  0
  8  9  1  0  0  0  0
  9 24  1  0  0  0  0
  9 38  1  0  0  0  0
 10 38  1  0  0  0  0
 11 12  1  0  0  0  0
 12 17  1  0  0  0  0
 12 27  1  0  0  0  0
 12 38  1  0  0  0  0
 13 24  1  0  0  0  0
 14 27  1  0  0  0  0
 15 26  1  0  0  0  0
 16 17  1  0  0  0  0
 17 26  1  0  0  0  0
 17 37  1  0  0  0  0
 18 37  1  0  0  0  0
 19 22  1  0  0  0  0
 20 21  1  0  0  0  0
 21 22  1  0  0  0  0
 21 37  1  0  0  0  0
 22 26  1  0  0  0  0
 22 34  1  0  0  0  0
 23 35  1  0  0  0  0
 24 27  1  0  0  0  0
 24 35  1  0  0  0  0
 25 29  1  0  0  0  0
 25 35  1  0  0  0  0
 26 33  1  0  0  0  0
 27 33  1  0  0  0  0
 28 33  1  0  0  0  0
 29 34  1  0  0  0  0
 30 34  1  0  0  0  0
 31 32  1  0  0  0  0
 31 33  1  0  0  0  0
 31 35  1  0  0  0  0
 31 34  1  0  0  0  0
 36 37  1  0  0  0  0
 36 38  1  0  0  0  0
 36 39  1  0  0  0  0
 36 40  1  0  0  0  0
M  END

  )";

    auto mol =
        std::unique_ptr<RDKit::RWMol>(MolBlockToMol(molBlock, true, false));
    REQUIRE(mol);
    CHECK_THROWS_AS(CIPLabeler::assignCIPLabels(*mol, 100000),
                    CIPLabeler::MaxIterationsExceeded);

    // try a second call to test that the timer is reset

    CHECK_THROWS_AS(CIPLabeler::assignCIPLabels(*mol, 100000),
                    CIPLabeler::MaxIterationsExceeded);
  }
}

void testOneAtropIomerMandP(std::string inputText, const std::string &expected,
                            bool isSmiles = true) {
  std::unique_ptr<RWMol> mol;

  if (isSmiles) {
    SmilesParserParams smilesParserParams;
    smilesParserParams.sanitize = true;
    smilesParserParams.removeHs = false;

    mol.reset(SmilesToMol(inputText, smilesParserParams));
  } else {
    mol.reset(MrvBlockToMol(inputText));
  }
  REQUIRE(mol);
  CIPLabeler::assignCIPLabels(*mol, 100000);

  std::ostrstream out;
  bool foundOne = false;
  for (auto bond : mol->bonds()) {
    if (bond->hasProp(common_properties::_CIPCode)) {
      out << bond->getBeginAtomIdx() << "-" << bond->getEndAtomIdx() << "="
          << bond->getProp<std::string>(common_properties::_CIPCode) << ":";
      foundOne = true;
    }
  }

  if (!foundOne) {
    out << "none ";
  }
  out << std::ends;

  CHECK(out.str() == expected);
}

TEST_CASE("AssignMandP", "[accurateCIP]") {
  SECTION("Suzzana") {
    testOneAtropIomerMandP(
        "<cml xmlns=\"http://www.chemaxon.com\" version=\"ChemAxon file format v20.20.0, generated by vunknown\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.chemaxon.com http://www.chemaxon.com/marvin/schema/mrvSchema_20_20_0.xsd\"><MDocument><MChemicalStruct><molecule molID=\"m1\" absStereo=\"true\"><atomArray><atom id=\"a1\" elementType=\"C\" x2=\"4.070453300769707\" y2=\"-1.47671998818624\"/><atom id=\"a2\" elementType=\"C\" x2=\"4.372853298350507\" y2=\"0.03322666640085333\"/><atom id=\"a3\" elementType=\"C\" x2=\"5.8316532866801065\" y2=\"0.5265866624539733\"/><atom id=\"a4\" elementType=\"C\" x2=\"6.134053284260906\" y2=\"2.0365333170410667\"/><atom id=\"a5\" elementType=\"C\" x2=\"4.977466626846933\" y2=\"3.053306642240213\"/><atom id=\"a6\" elementType=\"O\" x2=\"5.27967995776256\" y2=\"4.56343996349248\" lonePair=\"2\"/><atom id=\"a7\" elementType=\"C\" x2=\"3.51847997185216\" y2=\"2.5601333128522668\"/><atom id=\"a8\" elementType=\"C\" x2=\"2.361893314438187\" y2=\"3.576906638051413\"/><atom id=\"a9\" elementType=\"C\" x2=\"3.2162666409365337\" y2=\"1.0501866582651733\"/><atom id=\"a10\" elementType=\"N\" x2=\"-1.0962933245629867\" y2=\"-0.40786666340373334\" lonePair=\"1\"/><atom id=\"a11\" elementType=\"C\" x2=\"-2.3527466478446932\" y2=\"0.48253332947306665\"/><atom id=\"a12\" elementType=\"N\" x2=\"-2.37103998103168\" y2=\"2.022346650487893\" lonePair=\"1\"/><atom id=\"a13\" elementType=\"C\" x2=\"-3.58791997129664\" y2=\"-0.43735999650112\"/><atom id=\"a14\" elementType=\"C\" x2=\"-3.09455997524352\" y2=\"-1.89615998483072\"/><atom id=\"a15\" elementType=\"C\" x2=\"-1.5547466542286932\" y2=\"-1.8780533183089065\"/></atomArray><bondArray><bond id=\"b1\" atomRefs2=\"a1 a2\" order=\"1\"/><bond id=\"b2\" atomRefs2=\"a2 a3\" order=\"2\"/><bond id=\"b3\" atomRefs2=\"a3 a4\" order=\"1\"/><bond id=\"b4\" atomRefs2=\"a4 a5\" order=\"2\"/><bond id=\"b5\" atomRefs2=\"a5 a6\" order=\"1\"/><bond id=\"b6\" atomRefs2=\"a5 a7\" order=\"1\"/><bond id=\"b7\" atomRefs2=\"a7 a8\" order=\"1\"/><bond id=\"b8\" atomRefs2=\"a7 a9\" order=\"2\"/><bond id=\"b9\" atomRefs2=\"a9 a2\" order=\"1\"><bondStereo>W</bondStereo></bond><bond id=\"b10\" atomRefs2=\"a9 a10\" order=\"1\"/><bond id=\"b11\" atomRefs2=\"a10 a11\" order=\"1\"/><bond id=\"b12\" atomRefs2=\"a11 a12\" order=\"1\"/><bond id=\"b13\" atomRefs2=\"a11 a13\" order=\"2\"/><bond id=\"b14\" atomRefs2=\"a13 a14\" order=\"1\"/><bond id=\"b15\" atomRefs2=\"a14 a15\" order=\"2\"/><bond id=\"b16\" atomRefs2=\"a10 a15\" order=\"1\"/></bondArray></molecule></MChemicalStruct></MDocument></cml>",
        "8-9=M:", false);
  }

  SECTION("Suzzana2") {
    testOneAtropIomerMandP(
        "<cml xmlns=\"http://www.chemaxon.com\" version=\"ChemAxon file format v20.20.0, generated by vunknown\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.chemaxon.com http://www.chemaxon.com/marvin/schema/mrvSchema_20_20_0.xsd\"><MDocument><MChemicalStruct><molecule molID=\"m1\" absStereo=\"true\"><atomArray><atom id=\"a1\" elementType=\"C\" x2=\"4.070453300769707\" y2=\"-1.47671998818624\"/><atom id=\"a2\" elementType=\"C\" x2=\"4.372853298350507\" y2=\"0.03322666640085333\"/><atom id=\"a3\" elementType=\"C\" x2=\"5.8316532866801065\" y2=\"0.5265866624539733\"/><atom id=\"a4\" elementType=\"C\" x2=\"6.134053284260906\" y2=\"2.0365333170410667\"/><atom id=\"a5\" elementType=\"C\" x2=\"4.977466626846933\" y2=\"3.053306642240213\"/><atom id=\"a6\" elementType=\"O\" x2=\"5.27967995776256\" y2=\"4.56343996349248\" lonePair=\"2\"/><atom id=\"a7\" elementType=\"C\" x2=\"3.51847997185216\" y2=\"2.5601333128522668\"/><atom id=\"a8\" elementType=\"C\" x2=\"2.361893314438187\" y2=\"3.576906638051413\"/><atom id=\"a9\" elementType=\"C\" x2=\"3.2162666409365337\" y2=\"1.0501866582651733\"/><atom id=\"a10\" elementType=\"N\" x2=\"-1.0962933245629867\" y2=\"-0.40786666340373334\" lonePair=\"1\"/><atom id=\"a11\" elementType=\"C\" x2=\"-2.3527466478446932\" y2=\"0.48253332947306665\"/><atom id=\"a12\" elementType=\"N\" x2=\"-2.37103998103168\" y2=\"2.022346650487893\" lonePair=\"1\"/><atom id=\"a13\" elementType=\"C\" x2=\"-3.58791997129664\" y2=\"-0.43735999650112\"/><atom id=\"a14\" elementType=\"C\" x2=\"-3.09455997524352\" y2=\"-1.89615998483072\"/><atom id=\"a15\" elementType=\"C\" x2=\"-1.5547466542286932\" y2=\"-1.8780533183089065\"/></atomArray><bondArray><bond id=\"b1\" atomRefs2=\"a1 a2\" order=\"1\"/><bond id=\"b2\" atomRefs2=\"a2 a3\" order=\"2\"/><bond id=\"b3\" atomRefs2=\"a3 a4\" order=\"1\"/><bond id=\"b4\" atomRefs2=\"a4 a5\" order=\"2\"/><bond id=\"b5\" atomRefs2=\"a5 a6\" order=\"1\"/><bond id=\"b6\" atomRefs2=\"a5 a7\" order=\"1\"/><bond id=\"b7\" atomRefs2=\"a7 a8\" order=\"1\"/><bond id=\"b8\" atomRefs2=\"a7 a9\" order=\"2\"/><bond id=\"b9\" atomRefs2=\"a9 a2\" order=\"1\"><bondStereo>H</bondStereo></bond><bond id=\"b10\" atomRefs2=\"a9 a10\" order=\"1\"/><bond id=\"b11\" atomRefs2=\"a10 a11\" order=\"1\"/><bond id=\"b12\" atomRefs2=\"a11 a12\" order=\"1\"/><bond id=\"b13\" atomRefs2=\"a11 a13\" order=\"2\"/><bond id=\"b14\" atomRefs2=\"a13 a14\" order=\"1\"/><bond id=\"b15\" atomRefs2=\"a14 a15\" order=\"2\"/><bond id=\"b16\" atomRefs2=\"a10 a15\" order=\"1\"/></bondArray></molecule></MChemicalStruct></MDocument></cml>",
        "8-9=P:", false);
  }

  SECTION("BMS-986142") {
    testOneAtropIomerMandP(
        "FC1=C(C2=C(C)C(N3C(=O)C4=C(N(C)C3=O)C(F)=CC=C4)=CC=C2)C2=C(NC3=C2CC[C@H](C(O)(C)C)C3)C(C(=O)N)=C1 |(2.1158,0.5489,;1.4029,0.9642,;0.6554,0.5402,;0.6459,-0.2846,;1.3556,-0.7053,;2.0747,-0.3011,;1.346,-1.5302,;2.3682,-1.8618,;2.9748,-1.3027,;2.794,-0.4978,;3.7623,-1.5486,;3.9431,-2.3535,;3.3364,-2.9126,;3.5173,-3.7175,;2.549,-2.6667,;1.9423,-3.2258,;4.7594,-2.6222,;4.9309,-3.4292,;5.3958,-2.0448,;5.2074,-1.2063,;4.3851,-0.9565,;0.6268,-1.9345,;-0.0827,-1.5137,;-0.0732,-0.6889,;-0.0819,0.9813,;-0.0819,1.8063,;-0.8665,2.0612,;-1.3515,1.3938,;-0.8665,0.7264,;-1.2039,-0.0639,;-2.0578,-0.1602,;-2.5629,0.5349,;-3.3837,0.4518,;-4.2045,0.3687,;-3.4668,1.2725,;-3.3006,-0.369,;-2.2074,1.3172,;0.6555,2.2474,;0.6459,3.0723,;-0.0731,3.4765,;1.3556,3.4931,;1.403,1.8235,),wU:7.14wD:31.35|",
        "6-7=P:");
  }
  SECTION("AtropManyChirals") {
    testOneAtropIomerMandP(
        "O=C1C=CNC(=O)N1C(=C)[C@]([C@H](C)Br)([C@@H](F)C)[C@H](Cl)C |(1.6038,0.5917,;1.2842,-0.2899,;1.9435,-0.5613,;2.0181,-1.3829,;1.344,-1.8584,;0.5951,-1.5123,;0.0595,-2.2018,;0.5204,-0.6907,;-0.3117,-0.3431,;-0.861,-1.0085,;-0.5013,0.5405,;-0.4918,1.4387,;-0.0783,2.2018,;-1.0961,2.0398,;0.3347,0.736,;0.6313,-0.0436,;1.118,1.0786,;-1.3746,0.4021,;-2.0181,0.9732,;-1.9327,-0.2608,),wU:7.7,11.13,10.14,17.18wD:14.15|",
        "7-8=M:");
  }
  SECTION("BMS-986142_3d") {
    testOneAtropIomerMandP(
        "FC1=C(C2=C(C([H])([H])[H])C(N3C(=O)C4=C(N(C([H])([H])[H])C3=O)C(F)=C([H])C([H])=C4[H])=C([H])C([H])=C2[H])C2=C(N([H])C3=C2C([H])([H])C([H])([H])C(C(O[H])(C([H])([H])[H])C([H])([H])[H])([H])C3([H])[H])C(C(=O)N([H])[H])=C1[H] |(-1.4248,-1.9619,-2.2208;-1.8741,-1.6494,-1.6034;-1.6332,-0.9186,-1.2186;-0.9031,-0.4893,-1.4844;-0.138,-0.6553,-1.1342;-0.0582,-1.2887,-0.4699;-0.2862,-1.8851,-0.6843;0.5723,-1.396,-0.2643;-0.4068,-1.1014,0.0747;0.5573,-0.2316,-1.4028;1.3588,-0.3932,-1.0482;1.5905,0.0773,-0.3585;1.2002,0.6554,-0.0665;2.3933,-0.1194,0.0095;2.9068,-0.6967,-0.3594;2.6425,-1.1179,-1.0794;3.1982,-1.7169,-1.457;3.3106,-2.2227,-1.0278;3.7704,-1.4129,-1.6385;2.9318,-1.9829,-2.0144;1.8575,-1.0047,-1.4246;1.6168,-1.4364,-2.0026;3.6768,-0.8489,0.0041;4.2097,-1.3952,-0.3066;3.9217,-0.4363,0.7167;4.5184,-0.5645,0.9877;3.4017,0.1354,1.0788;3.5896,0.457,1.6328;2.6386,0.2956,0.7262;2.2428,0.7475,1.0192;0.4878,0.3584,-2.0214;1.0254,0.6927,-2.2354;-0.2773,0.5242,-2.3716;-0.3315,0.9837,-2.8531;-0.9727,0.1003,-2.1031;-1.5642,0.2383,-2.3831;-2.1142,-0.6018,-0.5711;-2.8218,-1.0304,-0.3345;-3.1799,-0.5981,0.3004;-3.7005,-0.7644,0.5964;-2.7239,0.0857,0.4696;-2.0583,0.1027,-0.0552;-1.4396,0.7778,-0.0395;-0.8137,0.538,-0.0869;-1.5466,1.1855,-0.5657;-1.4992,1.2797,0.7621;-1.2099,0.9291,1.2614;-1.1302,1.8345,0.67;-2.4007,1.4836,0.9934;-2.4335,2.0277,1.7578;-2.112,1.5747,2.4382;-2.128,1.9287,2.9199;-1.9128,2.8014,1.6418;-2.0588,3.1185,1.0662;-1.2527,2.6718,1.666;-2.0211,3.2354,2.1466;-3.3209,2.2556,1.9643;-3.6328,2.5429,1.4415;-3.3403,2.6896,2.48;-3.6795,1.7235,2.1686;-2.6689,1.821,0.4719;-2.9136,0.6925,1.1267;-2.7776,0.4065,1.7216;-3.5722,0.834,1.1115;-3.069,-1.7619,-0.7143;-3.8054,-2.1982,-0.4602;-4.2532,-1.9338,0.1012;-4.0043,-2.9212,-0.8738;-3.6777,-3.1773,-1.3373;-4.5184,-3.2354,-0.7104;-2.581,-2.0639,-1.3539;-2.7202,-2.6244,-1.6865),wD:10.20|", "9-10=P:");
  }
  SECTION("JDQ443_atrop2") {
    testOneAtropIomerMandP(
        "ClC1=C(C2=C(C)N(C3CC4(C3)CN(C(=O)C=C)C4)N=C2C2=CC3=C(N(C)N=C3)C=C2)C2=C(NN=C2)C=C1C |(-2.3128,3.1804,;-2.4843,2.3736,;-1.8712,1.8214,;-1.0865,2.0764,;-0.8316,2.861,;-1.3166,3.5284,;-0.0066,2.861,;0.478,3.5284,;0.3491,4.3433,;1.164,4.4723,;1.2929,3.6576,;1.0348,5.2872,;1.8497,5.4162,;2.3346,6.0837,;1.9992,6.8374,;3.1551,5.9974,;3.6401,6.6649,;1.9788,4.6014,;0.2482,2.0765,;-0.4191,1.5915,;-0.0067,0.8771,;0.8182,0.877,;1.2308,0.1626,;0.8182,-0.5518,;1.3686,-1.1598,;1.1985,-1.9672,;2.1174,-0.8258,;2.0325,-0.01,;-0.0067,-0.5518,;-0.4193,0.1625,;-2.0427,1.0145,;-2.8273,0.7596,;-2.8254,-0.0605,;-2.0451,-0.3132,;-1.5624,0.3498,;-3.4405,1.3115,;-3.2689,2.1186,;-3.882,2.6705,),wD:3.21|",
        "2-3=M:");
  }
  SECTION("Mrtx1719") {
    testOneAtropIomerMandP(
        "ClC1=C(F)C(C2=C(C3=CC4=C(C=C3)C(=O)NN=C4CN)C=NN2C)=C(C#N)C(OC2CC2)=C1 |(0.9486,-2.0481,;0.2357,-1.6329,;0.2388,-0.808,;0.955,-0.3981,;-0.474,-0.3926,;-0.4707,0.4321,;0.1985,0.9146,;0.9821,0.6566,;1.254,1.4719,;2.0972,1.6376,;2.6493,1.0245,;2.3964,0.2031,;1.5571,0.0181,;3.4562,1.196,;4.0082,0.5829,;3.7112,1.9806,;3.1591,2.5937,;2.3522,2.4222,;1.8001,3.0353,;2.0551,3.8199,;-0.0533,1.7003,;-0.8784,1.7034,;-1.1363,0.9197,;-1.9219,0.6678,;-1.1899,-0.8024,;-1.9029,-0.3871,;-2.6157,0.0279,;-1.1932,-1.6273,;-1.9093,-2.0371,;-1.9124,-2.8622,;-2.3276,-3.5751,;-1.5026,-3.5781,;-0.4803,-2.0427,),wU:4.3|",
        "4-5=M:");
  }

  SECTION("RP-6306") {
    testOneAtropIomerMandP(
        "OC1=C(C)C(N2C3=C(C(C(=O)N)=C2N)C=C(C)C(C)=N3)=C(C)C=C1 |(2.0992,3.5589,;2.2694,2.7517,;1.6554,2.2008,;0.8712,2.4571,;1.8255,1.3935,;1.0637,0.4452,;0.8122,-0.3353,;-0.0126,-0.3353,;-0.264,0.4452,;-1.0482,0.7015,;-1.6621,0.1505,;-1.2182,1.5087,;0.3998,0.9266,;0.3998,1.7516,;-0.425,-1.0498,;-0.0126,-1.7642,;-0.425,-2.4786,;0.8122,-1.7642,;1.2247,-2.4786,;1.2247,-1.0498,;2.6096,1.1372,;2.7797,0.33,;3.2236,1.6882,;3.0535,2.4954,),wD:4.3|",
        "4-5=P:");
  }
  SECTION("Sotorasib") {
    testOneAtropIomerMandP(
        "FC1=C(C2=C(O)C=CC=C2F)N=C2N(C3=C(C)C=CN=C3C(C)C)C(=O)N=C(N3[C@H](C)CN(C(=O)C=C)CC3)C2=C1 |(1.4145,2.306,;1.4178,1.4809,;2.1619,1.0513,;2.8781,1.461,;2.8811,2.2859,;2.1683,2.7012,;3.5973,2.6956,;4.3101,2.2805,;4.3069,1.4555,;3.5908,1.0457,;3.5877,0.2207,;2.1553,0.1919,;1.4045,-0.226,;1.4045,-1.0512,;2.1191,-1.4637,;2.1191,-2.2887,;1.4045,-2.7012,;2.8335,-2.7012,;3.5479,-2.2887,;3.5479,-1.4637,;2.8336,-1.0512,;3.3638,-0.4191,;4.1763,-0.5623,;3.0816,0.3561,;0.6902,-1.4637,;0.6902,-2.2887,;-0.0243,-1.0512,;-0.0243,-0.2262,;-0.7388,0.1862,;-1.4531,-0.2262,;-1.4531,-1.0512,;-2.1677,0.1862,;-2.1677,1.0112,;-2.8821,1.4236,;-2.8821,2.2487,;-3.5966,1.0112,;-4.3112,1.4236,;-1.4531,1.4237,;-0.7386,1.0114,;0.6902,0.1864,;0.6768,1.0455,),wU:14.21wD:29.31|",
        "13-14=M:");
  }
  SECTION("ZM374979") {
    testOneAtropIomerMandP(
        "C(N(C[C@H](CCN1CCC(C2=C(S(C)=O)C=CC=C2)CC1)C1=CC(Cl)=C(Cl)C=C1)C)(=O)C1=C(CC)C(C#N)=CC2=C1C=CC=C2 |(2.5006,0.2061,;1.7861,-0.2061,;1.0716,0.2061,;0.3572,-0.2061,;-0.3572,0.2061,;-1.0716,-0.2061,;-1.7861,0.2061,;-1.7861,1.0313,;-2.5006,1.4438,;-3.215,1.0313,;-3.9296,1.4438,;-3.9296,2.2686,;-3.215,2.6811,;-2.5006,2.2686,;-3.215,3.5063,;-4.6441,2.6811,;-5.3584,2.2686,;-5.3584,1.4438,;-4.6441,1.0313,;-3.215,0.2061,;-2.5006,-0.2061,;0.3572,-1.0313,;-0.3572,-1.4438,;-0.3572,-2.2686,;-1.0716,-2.6811,;0.3572,-2.6811,;0.3572,-3.5063,;1.0716,-2.2686,;1.0716,-1.4438,;1.7861,-1.0313,;2.5006,1.0313,;3.215,-0.2061,;3.215,-1.0313,;2.5006,-1.4438,;2.5006,-2.2686,;3.9296,-1.4438,;3.9296,-2.2686,;3.9296,-3.0938,;4.6441,-1.0313,;4.6441,-0.2061,;3.9296,0.2061,;3.9296,1.0313,;4.6441,1.4438,;5.3584,1.0313,;5.3584,0.2061,),wU:3.22wD:0.0|",
        "0-31=P:");
  }
  SECTION("two atropisomers") {
    // auto m = MolFileToMol(
    testOneAtropIomerMandP(
        "C1(N2C(C)=CC=C2Br)=C(C)C(C)=C(N2C(C)=CC=C2Br)C(C)=C1C |(-0.0002,1.5403,;-0.0002,3.0805,;-1.334,3.8508,;-2.6678,3.0807,;-1.334,5.391,;1.3338,5.391,;1.3338,3.8508,;2.6676,3.0807,;-1.3338,0.7702,;-2.6678,1.5403,;-1.3338,-0.7702,;-2.6678,-1.5401,;-0.0002,-1.5403,;-0.0002,-3.0805,;1.3338,-3.8508,;2.6676,-3.0805,;1.3338,-5.391,;-1.334,-5.391,;-1.334,-3.8508,;-2.6678,-3.0805,;1.3338,-0.7702,;2.6678,-1.5403,;1.3338,0.7702,;2.6678,1.5404,),wU:1.6,13.14|",
        "0-1=m:12-13=m:");
  }
  SECTION("psuedo") {
    testOneAtropIomerMandP(
        "N1(n2c(C)ccc2Br)C(=O)[C@H](C)[C@H](C)C1=O |(-11.1517,1.8306,;-11.1517,3.3708,;-12.4855,4.1411,;-13.8193,3.371,;-12.4855,5.6813,;-9.8177,5.6813,;-9.8177,4.1411,;-8.4839,3.371,;-12.3975,0.9252,;-13.8622,1.4011,;-11.9217,-0.5394,;-12.8269,-1.7852,;-10.3817,-0.5394,;-9.4765,-1.7852,;-9.9059,0.9252,;-8.4413,1.4011,),wU:0.8,10.11,12.13|",
        "0-1=p:");
  }
}

TEST_CASE("upsertTest", "[basic]") {
  SECTION("upsertTest1") {
    auto ps = SmilesParserParams();
    ps.allowCXSMILES = true;
    ps.sanitize = true;
    ps.removeHs = false;

    auto m = SmilesToMol(
        "CC1=C(C2=C(F)C=C(C(N)=O)C3=C2C2=C(C[C@@H](C(C)(C)O)CC2)N3)C=CC=C1N1C(=O)C2=C(C(F)=CC=C2)N(C)C1=O |(-0.0582,-1.2887,-0.4699;-0.138,-0.6553,-1.1342;-0.9031,-0.4893,-1.4844;-1.6332,-0.9186,-1.2186;-1.8741,-1.6494,-1.6034;-1.4248,-1.9619,-2.2208;-2.581,-2.0639,-1.3539;-3.069,-1.7619,-0.7143;-3.8054,-2.1982,-0.4602;-4.0043,-2.9212,-0.8738;-4.2532,-1.9338,0.1012;-2.8218,-1.0304,-0.3345;-2.1142,-0.6018,-0.5711;-2.0583,0.1027,-0.0552;-2.7239,0.0857,0.4696;-2.9136,0.6925,1.1267;-2.4007,1.4836,0.9934;-2.4335,2.0277,1.7578;-1.9128,2.8014,1.6418;-3.3209,2.2556,1.9643;-2.112,1.5747,2.4382;-1.4992,1.2797,0.7621;-1.4396,0.7778,-0.0395;-3.1799,-0.5981,0.3004;-0.9727,0.1003,-2.1031;-0.2773,0.5242,-2.3716;0.4878,0.3584,-2.0214;0.5573,-0.2316,-1.4028;1.3588,-0.3932,-1.0482;1.5905,0.0773,-0.3585;1.2002,0.6554,-0.0665;2.3933,-0.1194,0.0095;2.9068,-0.6967,-0.3594;3.6768,-0.8489,0.0041;4.2097,-1.3952,-0.3066;3.9217,-0.4363,0.7167;3.4017,0.1354,1.0788;2.6386,0.2956,0.7262;2.6425,-1.1179,-1.0794;3.1982,-1.7169,-1.457;1.8575,-1.0047,-1.4246;1.6168,-1.4364,-2.0026),wU:27.30|",
        ps);

    unsigned int flags = RDKit::SmilesWrite::CXSmilesFields::CX_COORDS |
                         RDKit::SmilesWrite::CXSmilesFields::CX_MOLFILE_VALUES |
                         RDKit::SmilesWrite::CXSmilesFields::CX_ATOM_PROPS |
                         RDKit::SmilesWrite::CXSmilesFields::CX_BOND_CFG |
                         RDKit::SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO;

    auto smilesWriteParams = SmilesWriteParams();
    smilesWriteParams.canonical = true;
    auto smilesOut =
        MolToCXSmiles(*m, smilesWriteParams, flags, RestoreBondDirOptionClear);
    CHECK(smilesOut != "");
  }
}

TEST_CASE("atropisomers", "[basic]") {
  SECTION("atropisomers1") {
    UseLegacyStereoPerceptionFixture useLegacy(false);

    std::string rdbase = getenv("RDBASE");

    std::vector<std::string> files = {
        "BMS-986142_atrop5.sdf", "BMS-986142_atrop1.sdf",
        "BMS-986142_atrop2.sdf", "JDQ443_atrop1.sdf",
        "JDQ443_atrop2.sdf",     "Mrtx1719_atrop1.sdf",
        "Mrtx1719_atrop2.sdf",   "RP-6306_atrop1.sdf",
        "RP-6306_atrop2.sdf",    "Sotorasib_atrop1.sdf",
        "Sotorasib_atrop2.sdf",  "ZM374979_atrop1.sdf",
        "ZM374979_atrop2.sdf"};

    for (auto file : files) {
      auto fName =
          rdbase + "/Code/GraphMol/FileParsers/test_data/atropisomers/" + file;

      auto molsdf = MolFileToMol(fName, true, true, true);
      REQUIRE(molsdf);
      CIPLabeler::assignCIPLabels(*molsdf, 100000);

      std::map<std::pair<unsigned int, unsigned int>, std::string> CIPVals;
      for (auto bond : molsdf->bonds()) {
        auto a1 = bond->getBeginAtomIdx();
        auto a2 = bond->getEndAtomIdx();
        if (a1 > a2) {
          std::swap(a1, a2);
        }

        if (bond->hasProp(common_properties::_CIPCode)) {
          CIPVals[std::make_pair(a1, a2)] =
              bond->getProp<std::string>(common_properties::_CIPCode);
        } else {
          CIPVals[std::make_pair(a1, a2)] = "N/A";
        }
      }

      molsdf->clearConformers();

      SmilesWriteParams wp;
      wp.canonical = false;
      wp.doKekule = false;
      unsigned int flags =
          RDKit::SmilesWrite::CXSmilesFields::CX_MOLFILE_VALUES |
          RDKit::SmilesWrite::CXSmilesFields::CX_ATOM_PROPS |
          RDKit::SmilesWrite::CXSmilesFields::CX_BOND_CFG |
          RDKit::SmilesWrite::CXSmilesFields::CX_BOND_ATROPISOMER;
      SmilesParserParams pp;
      pp.allowCXSMILES = true;
      pp.sanitize = true;

      auto smi = MolToCXSmiles(*molsdf, wp, flags);
      auto newMol = SmilesToMol(smi, pp);
      CIPLabeler::assignCIPLabels(*newMol, 100000);

      std::map<std::pair<unsigned int, unsigned int>, std::string> newCIPVals;
      for (auto bond : newMol->bonds()) {
        auto a1 = bond->getBeginAtomIdx();
        auto a2 = bond->getEndAtomIdx();
        if (a1 > a2) {
          std::swap(a1, a2);
        }
        if (bond->hasProp(common_properties::_CIPCode)) {
          newCIPVals[std::make_pair(a1, a2)] =
              bond->getProp<std::string>(common_properties::_CIPCode);
        } else {
          newCIPVals[std::make_pair(a1, a2)] = "N/A";
        }
      }

      RDKit::SubstructMatchParameters params;

      auto match = RDKit::SubstructMatch(*molsdf, *newMol, params);

      for (auto thisBond : newMol->bonds()) {
        unsigned int a1 = thisBond->getBeginAtomIdx();
        unsigned int a2 = thisBond->getEndAtomIdx();
        if (a1 > a2) {
          std::swap(a1, a2);
        }

        unsigned int a1match = match[0][a1].second;
        unsigned int a2match = match[0][a2].second;

        if (a1match > a2match) {
          std::swap(a1match, a2match);
        }
        CHECK(CIPVals[std::make_pair(a1match, a2match)] ==
              newCIPVals[std::make_pair(a1, a2)]);
      }
    }
  }
}