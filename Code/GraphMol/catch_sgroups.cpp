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
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include "GraphMol/FileParsers/FileParsers.h"

using namespace RDKit;

/* Auxiliary functions */
void testIdxVector(const std::vector<unsigned int> &groupVector,
                   const std::vector<unsigned int> &reference) {
  size_t vecSize = reference.size();
  CHECK(groupVector.size() == vecSize);

  auto sgItr = groupVector.begin();
  for (auto refItr = reference.begin(); refItr != reference.end();
       ++sgItr, ++refItr) {
    CHECK(1 + *sgItr == *refItr);
  }
}

void testBrackets(
    const std::vector<SubstanceGroup::Bracket> &brackets,
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

RWMol buildSampleMolecule() {
  // This builds a RDKit::RWMol with all implemented SubstanceGroup features in
  // order to test them. SubstanceGroups and features probably do not make any
  // sense.

  //// Initialize Molecule ////
  RWMol mol;

  // Add some atoms and bonds
  for (unsigned i = 0; i < 6; ++i) {
    mol.addAtom(new Atom(6), false, true);

    if (i > 0) {
      mol.addBond(i - 1, i, Bond::SINGLE);
    }
  }

  //// First SubstanceGroup ////
  {
    SubstanceGroup sg(&mol, "MUL");

    sg.setProp("SUBTYPE", "BLO");
    sg.setProp("MULT", "n");
    sg.setProp("CONNECT", "HH");

    // Add some atoms and bonds
    for (unsigned i = 0; i < 3; ++i) {
      sg.addAtomWithIdx(i);
      sg.addParentAtomWithIdx(i);
      sg.addBondWithIdx(i);  // add 2 CBONDs + 1 XBOND
    }

    sg.setProp("COMPNO", 7u);
    sg.setProp("ESTATE", "E");

    SubstanceGroup::Bracket bracket1 = {{RDGeom::Point3D(1., 3., 0.),
                                         RDGeom::Point3D(5., 7., 0.),
                                         RDGeom::Point3D(0., 0., 0.)}};
    sg.addBracket(bracket1);

    SubstanceGroup::Bracket bracket2 = {{RDGeom::Point3D(2., 4., 0.),
                                         RDGeom::Point3D(6., 8., 0.),
                                         RDGeom::Point3D(0., 0., 0.)}};
    sg.addBracket(bracket2);

    // Vector should not be parsed (not a SUP group)
    sg.addCState(2, RDGeom::Point3D());

    sg.setProp("CLASS", "TEST CLASS");

    sg.addAttachPoint(0, 0, "XX");

    sg.setProp("BRKTYP", "PAREN");

    addSubstanceGroup(mol, sg);
  }
  //// Second SubstanceGroup ////
  {
    SubstanceGroup sg(&mol, "SUP");

    // Add some atoms and bonds
    for (unsigned i = 3; i < 6; ++i) {
      sg.addAtomWithIdx(i);
      sg.addParentAtomWithIdx(i);
      sg.addBondWithIdx(i - 1);  // add 1 XBOND + 2 CBONDs
    }

    sg.setProp("LABEL", "TEST LABEL");

    // V2000 has only x and y coords; z value restricted to 0.
    RDGeom::Point3D vector(3., 4., 0.);
    sg.addCState(2, vector);  // Vector should be parsed now!

    sg.addAttachPoint(3, -1, "YY");

    addSubstanceGroup(mol, sg);
  }
  //// Third SubstanceGroup ////
  {
    SubstanceGroup sg(&mol, "DAT");

    sg.setProp("FIELDNAME", "SAMPLE FIELD NAME");  // 30 char max
    // Field Type is ignored in V3000
    sg.setProp("FIELDINFO", "SAMPLE FIELD INFO");  // 20 char max
    sg.setProp("QUERYTYPE", "PQ");                 // 2 char max
    sg.setProp("QUERYOP", "SAMPLE QUERY OP");      // 15 char max (rest of line)

    // This should be properly formatted, but format is not checked
    sg.setProp("FIELDDISP", "SAMPLE FIELD DISP");

    STR_VECT dataFields = {"SAMPLE DATA FIELD 1", "SAMPLE DATA FIELD 2",
                           "SAMPLE DATA FIELD 3"};
    sg.setProp("DATAFIELDS", dataFields);

    addSubstanceGroup(mol, sg);
  }

  // Set a parent with higher index
  const auto &sgroups = getSubstanceGroups(mol);
  sgroups.at(0).setProp<unsigned int>("PARENT", 2);

  return mol;
}

/* End Auxiliary functions */

TEST_CASE("Basic Sgroup creation", "[Sgroups]") {
  // Create two SubstanceGroups and add them to a molecule
  RWMol mol;

  {
    SubstanceGroup sg0(&mol, "DAT");
    SubstanceGroup sg1(&mol, "SUP");
    addSubstanceGroup(mol, sg0);
    addSubstanceGroup(mol, sg1);
  }

  const auto &sgroups = getSubstanceGroups(mol);
  REQUIRE(sgroups.size() == 2);
  CHECK(sgroups.at(0).getProp<std::string>("TYPE") == "DAT");
  CHECK(sgroups.at(1).getProp<std::string>("TYPE") == "SUP");
}

TEST_CASE("Build and test sample molecule", "[Sgroups]") {
  auto mol = buildSampleMolecule();

  const auto &sgroups = getSubstanceGroups(mol);
  CHECK(sgroups.size() == 3);

  SECTION("first sgroup") {
    const auto &sg = sgroups.at(0);
    CHECK(sg.getProp<std::string>("TYPE") == "MUL");

    CHECK(sg.getProp<std::string>("SUBTYPE") == "BLO");
    CHECK(sg.getProp<std::string>("MULT") == "n");
    CHECK(sg.getProp<std::string>("CONNECT") == "HH");

    std::vector<unsigned int> atoms_reference = {1, 2, 3};
    auto atoms = sg.getAtoms();
    testIdxVector(atoms, atoms_reference);

    std::vector<unsigned int> patoms_reference = {1, 2, 3};
    testIdxVector(sg.getParentAtoms(), patoms_reference);

    std::vector<unsigned int> bonds_reference = {1, 2, 3};
    auto bonds = sg.getBonds();

    // bonds are not sorted in V3000; sort them here
    std::sort(bonds.begin(), bonds.end());

    testIdxVector(bonds, bonds_reference);

    CHECK(sg.getBondType(bonds[0]) == SubstanceGroup::BondType::CBOND);
    CHECK(sg.getBondType(bonds[1]) == SubstanceGroup::BondType::CBOND);
    CHECK(sg.getBondType(bonds[2]) == SubstanceGroup::BondType::XBOND);

    CHECK(sg.getProp<unsigned int>("COMPNO") == 7u);
    CHECK(sg.getProp<std::string>("ESTATE") == "E");

    std::vector<std::array<std::array<double, 3>, 3>> brackets_reference = {{
        {{{{1., 3., 0.}}, {{5., 7., 0.}}, {{0., 0., 0.}}}},
        {{{{2., 4., 0.}}, {{6., 8., 0.}}, {{0., 0., 0.}}}},
    }};
    testBrackets(sg.getBrackets(), brackets_reference);

    auto cstates = sg.getCStates();
    CHECK(cstates.size() == 1);
    CHECK(cstates[0].bondIdx == bonds[2]);
    CHECK(cstates[0].vector.x == 0.);
    CHECK(cstates[0].vector.y == 0.);
    CHECK(cstates[0].vector.z == 0.);

    CHECK(sg.getProp<std::string>("CLASS") == "TEST CLASS");

    auto ap = sg.getAttachPoints();
    CHECK(ap.size() == 1);
    CHECK(ap[0].aIdx == atoms[0]);
    CHECK(ap[0].lvIdx == static_cast<int>(atoms[0]));
    CHECK(ap[0].id == "XX");

    CHECK(sg.getProp<std::string>("BRKTYP") == "PAREN");

    CHECK(sg.getProp<unsigned int>("PARENT") == 2u);
  }

  SECTION("second sgroup") {
    const auto &sg = sgroups.at(1);
    CHECK(sg.getProp<std::string>("TYPE") == "SUP");

    std::vector<unsigned int> atoms_reference = {4, 5, 6};
    auto atoms = sg.getAtoms();
    testIdxVector(atoms, atoms_reference);

    std::vector<unsigned int> patoms_reference = {4, 5, 6};
    testIdxVector(sg.getParentAtoms(), patoms_reference);

    std::vector<unsigned int> bonds_reference = {3, 4, 5};
    auto bonds = sg.getBonds();

    // bonds are not sorted in V3000; sort them here
    std::sort(bonds.begin(), bonds.end());

    testIdxVector(bonds, bonds_reference);
    CHECK(sg.getBondType(bonds[0]) == SubstanceGroup::BondType::XBOND);
    CHECK(sg.getBondType(bonds[1]) == SubstanceGroup::BondType::CBOND);
    CHECK(sg.getBondType(bonds[2]) == SubstanceGroup::BondType::CBOND);

    CHECK(sg.getProp<std::string>("LABEL") == "TEST LABEL");

    auto cstates = sg.getCStates();
    CHECK(cstates.size() == 1);
    CHECK(cstates[0].bondIdx == bonds[0]);
    CHECK(cstates[0].vector.x == 3.);
    CHECK(cstates[0].vector.y == 4.);
    CHECK(cstates[0].vector.z == 0.);

    auto ap = sg.getAttachPoints();
    CHECK(ap.size() == 1);
    CHECK(ap[0].aIdx == atoms[0]);
    CHECK(ap[0].lvIdx == -1);
    CHECK(ap[0].id == "YY");
  }

  SECTION("third sgroup") {
    const auto &sg = sgroups.at(2);
    CHECK(sg.getProp<std::string>("TYPE") == "DAT");

    CHECK(sg.getProp<std::string>("FIELDNAME") == "SAMPLE FIELD NAME");
    CHECK(sg.getProp<std::string>("FIELDINFO") == "SAMPLE FIELD INFO");
    CHECK(sg.getProp<std::string>("QUERYTYPE") == "PQ");
    CHECK(sg.getProp<std::string>("QUERYOP") == "SAMPLE QUERY OP");

    CHECK(sg.getProp<std::string>("FIELDDISP") == "SAMPLE FIELD DISP");

    auto dataFields = sg.getProp<STR_VECT>("DATAFIELDS");
    CHECK(dataFields.size() == 3);
    CHECK(dataFields[0] == "SAMPLE DATA FIELD 1");
    CHECK(dataFields[1] == "SAMPLE DATA FIELD 2");
    CHECK(dataFields[2] == "SAMPLE DATA FIELD 3");
  }
}
