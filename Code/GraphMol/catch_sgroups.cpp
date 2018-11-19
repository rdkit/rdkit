//
//
//  Copyright (C) 2018 Greg Landrum and T5 Informatics GmbH
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
#include <GraphMol/Sgroup.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;


/* Auxiliary functions */
template <class T>
void testIdxVector(const std::vector<T> &groupVector,
                   const std::vector<unsigned int> &reference) {
  size_t vecSize = reference.size();
  CHECK(groupVector.size() == vecSize);

  auto sgItr = groupVector.begin();
  for (auto refItr = reference.begin(); refItr != reference.end();
       ++sgItr, ++refItr) {
    CHECK(1 + (*sgItr)->getIdx() == *refItr);
  }
}

void testBrackets(
    const std::vector<SGroup::Bracket> &brackets,
    const std::vector<std::array<std::array<double, 3>, 3>> &reference) {
  CHECK(brackets.size() == 2);
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        CHECK(std::abs(brackets[i][j][k] - reference[i][j][k]) < 1.e-6);
      }
    }
  }
}

std::shared_ptr<RWMol> buildSampleMolecule() {
  // This builds a RDKit::RWMol with all implemented SGroup features in order to
  // test them. SGroups and features probably do not make any sense.

  //// Initialize Molecule ////
  auto m = std::make_shared<RWMol>();

  // Add some atoms and bonds
  for (unsigned i = 0; i < 6; ++i) {
    m->addAtom(new Atom(6), false, true);

    if (i > 0) {
      m->addBond(i - 1, i, Bond::SINGLE);
    }
  }

  //// First SGroup ////
  {
    SGroup *sg = new SGroup(m.get(), "MUL");
    m->addSGroup(sg);

    sg->setProp("SUBTYPE", "BLO");
    sg->setProp("MULT", "n");
    sg->setProp("CONNECT", "HH");

    // Add some atoms and bonds
    for (unsigned i = 0; i < 3; ++i) {
      sg->addAtomWithIdx(i);
      sg->addPAtomWithIdx(i);
      sg->addBondWithIdx(i);  // add 2 CBONDs + 1 XBOND
    }

    sg->setCompNo(7);
    sg->setProp("ESTATE", "E");

    SGroup::Bracket bracket1 = {{RDGeom::Point3D(1., 3., 0.),
                                 RDGeom::Point3D(5., 7., 0.),
                                 RDGeom::Point3D(0., 0., 0.)}};
    sg->addBracket(bracket1);

    SGroup::Bracket bracket2 = {{RDGeom::Point3D(2., 4., 0.),
                                 RDGeom::Point3D(6., 8., 0.),
                                 RDGeom::Point3D(0., 0., 0.)}};
    sg->addBracket(bracket2);
    sg->addCState(m->getBondWithIdx(2),
                  nullptr);  // Vector should not be parsed (not a SUP group)

    sg->setProp("CLASS", "TEST CLASS");

    sg->addAttachPoint(m->getAtomWithIdx(0), m->getAtomWithIdx(0), "XX");

    sg->setProp("BRKTYP", "PAREN");
  }
  //// Second SGroup ////
  {
    SGroup *sg = new SGroup(m.get(), "SUP");
    m->addSGroup(sg);

    // Add some atoms and bonds
    for (unsigned i = 3; i < 6; ++i) {
      sg->addAtomWithIdx(i);
      sg->addPAtomWithIdx(i);
      sg->addBondWithIdx(i - 1);  // add 1 XBOND + 2 CBONDs
    }

    sg->setProp("LABEL", "TEST LABEL");

    // V2000 has only x and y coords; z value restricted to 0.
    auto *vector = new RDGeom::Point3D(3., 4., 0.);
    sg->addCState(m->getBondWithIdx(2),
                  vector);  // Vector should be parsed now!

    sg->addAttachPoint(m->getAtomWithIdx(3), nullptr, "YY");
  }
  //// Third SGroup ////
  {
    SGroup *sg = new SGroup(m.get(), "DAT");
    m->addSGroup(sg);

    sg->setProp("FIELDNAME", "SAMPLE FIELD NAME");  // 30 char max
    // Field Type is ignored in V3000
    sg->setProp("FIELDINFO", "SAMPLE FIELD INFO");  // 20 char max
    sg->setProp("QUERYTYPE", "PQ");                 // 2 char max
    sg->setProp("QUERYOP", "SAMPLE QUERY OP");  // 15 char max (rest of line)

    // This should be properly formatted, but format is not checked
    sg->setProp("FIELDDISP", "SAMPLE FIELD DISP");

    sg->addDataField("SAMPLE DATA FIELD 1");
    sg->addDataField("SAMPLE DATA FIELD 2");
    sg->addDataField("SAMPLE DATA FIELD 3");
  }

  // Set a parent with higher index
  auto sg0 = getMolSGroup(*m,0);
  auto sg2 = getMolSGroup(*m,2);
  sg0->setParent(sg2);

  return m;
}

/* End Auxiliary functions */

TEST_CASE("Basic Sgroup creation", "[Sgroups]") {
// Create two SGroups and add them to a molecule
  RWMol m;

  SGroup *sg = new SGroup(&m, "DAT");
  m.addSGroup(sg);

  sg = new SGroup(&m, "SUP");
  m.addSGroup(sg);

  sg = nullptr;

  REQUIRE(getMolNumSGroups(m) == 2);
  CHECK(getMolSGroup(m,0)->getType() == "DAT");
  CHECK(getMolSGroup(m,1)->getType() == "SUP");
}


TEST_CASE("Build and test sample molecule", "[Sgroups]") {
  auto mol = buildSampleMolecule();
  REQUIRE(mol);
  CHECK(getMolNumSGroups(*mol) == 3);
  
  SECTION("first sgroup"){
    auto sg = getMolSGroup(*mol,0);
    CHECK(sg->getType() == "MUL");

    CHECK(sg->getProp("SUBTYPE") == "BLO");
    CHECK(sg->getProp("MULT") == "n");
    CHECK(sg->getProp("CONNECT") == "HH");

    std::vector<unsigned int> atoms_reference = {1, 2, 3};
    auto atoms = sg->getAtoms();
    testIdxVector(atoms, atoms_reference);

    std::vector<unsigned int> patoms_reference = {1, 2, 3};
    testIdxVector(sg->getPAtoms(), patoms_reference);

    std::vector<unsigned int> bonds_reference = {1, 2, 3};
    auto bonds = sg->getBonds();

    // bonds are not sorted in V3000; sort them here
    auto cmpOutputIdx = [](Bond *a, Bond *b) {
      return a->getIdx() < b->getIdx();
    };
    std::sort(bonds.begin(), bonds.end(), cmpOutputIdx);

    testIdxVector(bonds, bonds_reference);

    CHECK(sg->getBondType(bonds[0]) == SGroup::BondType::CBOND);
    CHECK(sg->getBondType(bonds[1]) == SGroup::BondType::CBOND);
    CHECK(sg->getBondType(bonds[2]) == SGroup::BondType::XBOND);

    CHECK(sg->getCompNo() == 7);
    CHECK(sg->getProp("ESTATE") == "E");

    std::vector<std::array<std::array<double, 3>, 3>> brackets_reference = {{
        {{{{1., 3., 0.}}, {{5., 7., 0.}}, {{0., 0., 0.}}}},
        {{{{2., 4., 0.}}, {{6., 8., 0.}}, {{0., 0., 0.}}}},
    }};
    testBrackets(sg->getBrackets(), brackets_reference);

    auto cstates = sg->getCStates();
    CHECK(cstates.size() == 1);
    CHECK(cstates[0].bond == bonds[2]);
    CHECK(cstates[0].vector == nullptr);

    CHECK(sg->getProp("CLASS") == "TEST CLASS");

    auto ap = sg->getAttachPoints();
    CHECK(ap.size() == 1);
    CHECK(ap[0].aAtom == atoms[0]);
    CHECK(ap[0].lvAtom == atoms[0]);
    CHECK(ap[0].id == "XX");

    CHECK(sg->getProp("BRKTYP") == "PAREN");

    auto parent = sg->getParent();
    CHECK(parent == getMolSGroup(*mol,2));
  }

  SECTION("second sgroup"){
    auto sg = getMolSGroup(*mol,1);
    CHECK(sg->getType() == "SUP");

    std::vector<unsigned int> atoms_reference = {4, 5, 6};
    auto atoms = sg->getAtoms();
    testIdxVector(atoms, atoms_reference);

    std::vector<unsigned int> patoms_reference = {4, 5, 6};
    testIdxVector(sg->getPAtoms(), patoms_reference);

    std::vector<unsigned int> bonds_reference = {3, 4, 5};
    auto bonds = sg->getBonds();

    // bonds are not sorted in V3000; sort them here
    auto cmpOutputIdx = [](Bond *a, Bond *b) {
      return a->getIdx() < b->getIdx();
    };
    std::sort(bonds.begin(), bonds.end(), cmpOutputIdx);

    testIdxVector(bonds, bonds_reference);
    CHECK(sg->getBondType(bonds[0]) == SGroup::BondType::XBOND);
    CHECK(sg->getBondType(bonds[1]) == SGroup::BondType::CBOND);
    CHECK(sg->getBondType(bonds[2]) == SGroup::BondType::CBOND);

    CHECK(sg->getProp("LABEL") == "TEST LABEL");

    auto cstates = sg->getCStates();
    CHECK(cstates.size() == 1);
    CHECK(cstates[0].bond == bonds[0]);
    CHECK(cstates[0].vector != nullptr);
    CHECK(cstates[0].vector->x == 3.);
    CHECK(cstates[0].vector->y == 4.);
    CHECK(cstates[0].vector->z == 0.);

    auto ap = sg->getAttachPoints();
    CHECK(ap.size() == 1);
    CHECK(ap[0].aAtom == atoms[0]);
    CHECK(ap[0].lvAtom == nullptr);
    CHECK(ap[0].id == "YY");
  }

  SECTION("third sgroup"){
    auto sg = getMolSGroup(*mol,2);
    CHECK(sg->getType() == "DAT");

    CHECK(sg->getProp("FIELDNAME") == "SAMPLE FIELD NAME");
    CHECK(sg->getProp("FIELDINFO") == "SAMPLE FIELD INFO");
    CHECK(sg->getProp("QUERYTYPE") == "PQ");
    CHECK(sg->getProp("QUERYOP") == "SAMPLE QUERY OP");

    CHECK(sg->getProp("FIELDDISP") == "SAMPLE FIELD DISP");

    auto dataFields = sg->getDataFields();
    CHECK(dataFields.size() == 3);
    CHECK(dataFields[0] == "SAMPLE DATA FIELD 1");
    CHECK(dataFields[1] == "SAMPLE DATA FIELD 2");
    CHECK(dataFields[2] == "SAMPLE DATA FIELD 3");
  }

}

