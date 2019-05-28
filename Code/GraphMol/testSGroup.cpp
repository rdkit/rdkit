//
//
//  Copyright (C) 2002-2018 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SubstanceGroup.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include "GraphMol/FileParsers/FileParsers.h"
#include "GraphMol/FileParsers/MolSupplier.h"
#include "GraphMol/FileParsers/MolWriters.h"

#include <memory>
#include <cstdlib>

using namespace RDKit;

/* Auxiliary functions */
void testIdxVector(const std::vector<unsigned int> &groupVector,
                   const std::vector<unsigned int> &reference) {
  size_t vecSize = reference.size();
  TEST_ASSERT(groupVector.size() == vecSize);

  auto sgItr = groupVector.begin();
  for (auto refItr = reference.begin(); refItr != reference.end();
       ++sgItr, ++refItr) {
    TEST_ASSERT(1 + *sgItr == *refItr);
  }
}

void testBrackets(
    const std::vector<SubstanceGroup::Bracket> &brackets,
    const std::vector<std::array<std::array<double, 3>, 3>> &reference) {
  TEST_ASSERT(brackets.size() == 2);
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        TEST_ASSERT(std::abs(brackets[i][j][k] - reference[i][j][k]) < 1.e-6);
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

void checkSampleMolecule(const RWMol &mol) {
  // Test a molecule created by buildSampleMolecule (or a copy)

  const auto &sgroups = getSubstanceGroups(mol);
  TEST_ASSERT(sgroups.size() == 3);

  {
    // First SubstanceGroup
    const auto &sg = sgroups.at(0);
    TEST_ASSERT(sg.getProp<std::string>("TYPE") == "MUL");

    TEST_ASSERT(sg.getProp<std::string>("SUBTYPE") == "BLO");
    TEST_ASSERT(sg.getProp<std::string>("MULT") == "n");
    TEST_ASSERT(sg.getProp<std::string>("CONNECT") == "HH");

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

    TEST_ASSERT(sg.getBondType(bonds[0]) == SubstanceGroup::BondType::CBOND);
    TEST_ASSERT(sg.getBondType(bonds[1]) == SubstanceGroup::BondType::CBOND);
    TEST_ASSERT(sg.getBondType(bonds[2]) == SubstanceGroup::BondType::XBOND);

    TEST_ASSERT(sg.getProp<unsigned int>("COMPNO") == 7);
    TEST_ASSERT(sg.getProp<std::string>("ESTATE") == "E");

    std::vector<std::array<std::array<double, 3>, 3>> brackets_reference = {{
        {{{{1., 3., 0.}}, {{5., 7., 0.}}, {{0., 0., 0.}}}},
        {{{{2., 4., 0.}}, {{6., 8., 0.}}, {{0., 0., 0.}}}},
    }};
    testBrackets(sg.getBrackets(), brackets_reference);

    auto cstates = sg.getCStates();
    TEST_ASSERT(cstates.size() == 1);
    TEST_ASSERT(cstates[0].bondIdx == bonds[2]);
    TEST_ASSERT(cstates[0].vector.x == 0.);
    TEST_ASSERT(cstates[0].vector.y == 0.);
    TEST_ASSERT(cstates[0].vector.z == 0.);

    TEST_ASSERT(sg.getProp<std::string>("CLASS") == "TEST CLASS");

    auto ap = sg.getAttachPoints();
    TEST_ASSERT(ap.size() == 1);
    TEST_ASSERT(ap[0].aIdx == atoms[0]);
    TEST_ASSERT(ap[0].lvIdx == static_cast<int>(atoms[0]));
    TEST_ASSERT(ap[0].id == "XX");

    TEST_ASSERT(sg.getProp<std::string>("BRKTYP") == "PAREN");

    TEST_ASSERT(sg.getProp<unsigned int>("PARENT") == 2u);
  }

  {
    // Second SubstanceGroup
    const auto &sg = sgroups.at(1);
    TEST_ASSERT(sg.getProp<std::string>("TYPE") == "SUP");

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
    TEST_ASSERT(sg.getBondType(bonds[0]) == SubstanceGroup::BondType::XBOND);
    TEST_ASSERT(sg.getBondType(bonds[1]) == SubstanceGroup::BondType::CBOND);
    TEST_ASSERT(sg.getBondType(bonds[2]) == SubstanceGroup::BondType::CBOND);

    TEST_ASSERT(sg.getProp<std::string>("LABEL") == "TEST LABEL");

    auto cstates = sg.getCStates();
    TEST_ASSERT(cstates.size() == 1);
    TEST_ASSERT(cstates[0].bondIdx == bonds[0]);
    TEST_ASSERT(cstates[0].vector.x == 3.);
    TEST_ASSERT(cstates[0].vector.y == 4.);
    TEST_ASSERT(cstates[0].vector.z == 0.);

    auto ap = sg.getAttachPoints();
    TEST_ASSERT(ap.size() == 1);
    TEST_ASSERT(ap[0].aIdx == atoms[0]);
    TEST_ASSERT(ap[0].lvIdx == -1);
    TEST_ASSERT(ap[0].id == "YY");
  }

  {
    // Third SubstanceGroup
    const auto &sg = sgroups.at(2);
    TEST_ASSERT(sg.getProp<std::string>("TYPE") == "DAT");

    TEST_ASSERT(sg.getProp<std::string>("FIELDNAME") == "SAMPLE FIELD NAME");
    TEST_ASSERT(sg.getProp<std::string>("FIELDINFO") == "SAMPLE FIELD INFO");
    TEST_ASSERT(sg.getProp<std::string>("QUERYTYPE") == "PQ");
    TEST_ASSERT(sg.getProp<std::string>("QUERYOP") == "SAMPLE QUERY OP");

    TEST_ASSERT(sg.getProp<std::string>("FIELDDISP") == "SAMPLE FIELD DISP");

    auto dataFields = sg.getProp<STR_VECT>("DATAFIELDS");
    TEST_ASSERT(dataFields.size() == 3);
    TEST_ASSERT(dataFields[0] == "SAMPLE DATA FIELD 1");
    TEST_ASSERT(dataFields[1] == "SAMPLE DATA FIELD 2");
    TEST_ASSERT(dataFields[2] == "SAMPLE DATA FIELD 3");
  }
}

/* End Auxiliary functions */

void testCreateSubstanceGroups() {
  BOOST_LOG(rdInfoLog) << " ----------> Testing basic SubstanceGroup creation"
                       << std::endl;

  // Create two SubstanceGroups and add them to a molecule
  RWMol mol;

  {
    SubstanceGroup sg0(&mol, "DAT");
    SubstanceGroup sg1(&mol, "SUP");
    addSubstanceGroup(mol, sg0);
    addSubstanceGroup(mol, sg1);
  }

  const auto &sgroups = getSubstanceGroups(mol);
  TEST_ASSERT(sgroups.size() == 2);
  TEST_ASSERT(sgroups.at(0).getProp<std::string>("TYPE") == "DAT");
  TEST_ASSERT(sgroups.at(1).getProp<std::string>("TYPE") == "SUP");
}

void testParseSubstanceGroups(const std::string &rdbase) {
  BOOST_LOG(rdInfoLog) << " ----------> Parsing Issue3432136_1.mol (V2000)"
                       << std::endl;
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3432136_1.mol";

    std::unique_ptr<RWMol> mol(MolFileToMol(fName));
    TEST_ASSERT(mol);

    const auto &sgroups = getSubstanceGroups(*mol);
    TEST_ASSERT(sgroups.size() == 1);

    const auto &sgroup = sgroups.at(0);

    TEST_ASSERT(sgroup.getProp<std::string>("TYPE") == "MON");

    std::vector<unsigned int> atoms_reference = {2, 3, 4, 1, 5};

    testIdxVector(sgroup.getAtoms(), atoms_reference);

    std::vector<unsigned int> bonds_reference =
        {};  // No bonds defined in this mol
    testIdxVector(sgroup.getBonds(), bonds_reference);

    std::vector<std::array<std::array<double, 3>, 3>> brackets_reference = {{
        {{{{-3.9679, -0.1670, 0.}}, {{-3.9679, 2.1705, 0.}}, {{0., 0., 0.}}}},
        {{{{-0.7244, 2.1705, 0.}}, {{-0.7244, -0.1670, 0.}}, {{0., 0., 0.}}}},
    }};
    testBrackets(sgroup.getBrackets(), brackets_reference);
  }

  BOOST_LOG(rdInfoLog) << " ----------> Parsing Issue3432136_1.v3k.mol (V3000) "
                       << std::endl;
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3432136_1.v3k.mol";
    std::unique_ptr<RWMol> mol(MolFileToMol(fName));
    TEST_ASSERT(mol);

    const auto &sgroups = getSubstanceGroups(*mol);
    TEST_ASSERT(sgroups.size() == 1);

    const auto sgroup = sgroups.at(0);

    TEST_ASSERT(sgroup.getProp<std::string>("TYPE") == "MON");

    std::vector<unsigned int> atoms_reference = {2, 3, 4, 1, 5};
    testIdxVector(sgroup.getAtoms(), atoms_reference);

    std::vector<unsigned int> bonds_reference =
        {};  // No bonds defined in this mol
    testIdxVector(sgroup.getBonds(), bonds_reference);
  }

  BOOST_LOG(rdInfoLog) << " ----------> Parsing Issue3432136_2.v3k.mol (V3000) "
                       << std::endl;
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3432136_2.v3k.mol";
    std::unique_ptr<RWMol> mol(MolFileToMol(fName));
    TEST_ASSERT(mol);

    const auto &sgroups = getSubstanceGroups(*mol);
    TEST_ASSERT(sgroups.size() == 1);

    const auto sgroup = sgroups.at(0);

    TEST_ASSERT(sgroup.getProp<std::string>("TYPE") == "SUP");
    TEST_ASSERT(sgroup.getProp<std::string>("CLASS") == "DEMOCLASS");
    TEST_ASSERT(sgroup.getProp<std::string>("LABEL") == "abbrev");

    std::vector<unsigned int> atoms_reference = {6, 7, 8, 9, 11, 12};
    testIdxVector(sgroup.getAtoms(), atoms_reference);

    std::vector<unsigned int> bonds_reference = {5};
    testIdxVector(sgroup.getBonds(), bonds_reference);

    auto bond = sgroup.getBonds()[0];
    TEST_ASSERT(sgroup.getBondType(bond) == SubstanceGroup::BondType::XBOND);
  }

  BOOST_LOG(rdInfoLog) << " ----------> Parsing Issue3432136_2.mol (V2000) "
                       << std::endl;
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3432136_2.mol";
    std::unique_ptr<RWMol> mol(MolFileToMol(fName));
    TEST_ASSERT(mol);

    const auto &sgroups = getSubstanceGroups(*mol);
    TEST_ASSERT(sgroups.size() == 1);

    const auto sgroup = sgroups.at(0);

    TEST_ASSERT(sgroup.getProp<std::string>("TYPE") == "SUP");

    std::vector<unsigned int> atoms_reference = {6, 7, 8, 9, 11, 12};
    testIdxVector(sgroup.getAtoms(), atoms_reference);

    std::vector<unsigned int> bonds_reference = {5};
    testIdxVector(sgroup.getBonds(), bonds_reference);

    auto bond = sgroup.getBonds()[0];
    TEST_ASSERT(sgroup.getBondType(bond) == SubstanceGroup::BondType::XBOND);
  }
}

void testSubstanceGroupsRoundTrip(const std::string &rdbase, bool forceV3000) {
  BOOST_LOG(rdInfoLog)
      << " ----------> Testing SubstanceGroup writing & parsing Roundtrip ("
      << (forceV3000 ? "V3000" : "V2000") << ')' << std::endl;

  std::string fName =
      rdbase +
      "/Code/GraphMol/FileParsers/test_data/testSubstanceGroupsSample_" +
      (forceV3000 ? "V3000" : "V2000") + ".mol";

  {
    auto sampleMol = buildSampleMolecule();

    const auto &sgroups = getSubstanceGroups(sampleMol);
    TEST_ASSERT(sgroups.size() == 3);

    auto writer = SDWriter(fName);
    writer.setForceV3000(forceV3000);
    writer.write(sampleMol);
    writer.close();
  }

  std::unique_ptr<RWMol> roundtripMol(MolFileToMol(fName));
  checkSampleMolecule(*roundtripMol);
}

void testPickleSubstanceGroups() {
  BOOST_LOG(rdInfoLog)
      << " ----------> Testing SubstanceGroup pickling & unpickling Roundtrip"
      << std::endl;

  std::string pkl;

  {
    auto sampleMol = buildSampleMolecule();

    MolPickler::pickleMol(sampleMol, pkl);
  }

  RWMol roundtripMol(pkl);
  checkSampleMolecule(roundtripMol);
}

void testModifyMol() {
  BOOST_LOG(rdInfoLog)
      << " ----------> Test dropping SubstanceGroups on modification"
      << std::endl;

  auto mol = buildSampleMolecule();
  auto mol_copy = mol;

  const auto &sgroups = getSubstanceGroups(mol);
  TEST_ASSERT(sgroups.size() == 3);

  {  // insertion will drop SubstanceGroups
    const auto &sgroups = getSubstanceGroups(mol_copy);
    TEST_ASSERT(sgroups.size() == 3);

    mol_copy.insertMol(mol);

    TEST_ASSERT(sgroups.size() == 0);
  }
  {
    // adding an atom will drop SubstanceGroups
    mol_copy = mol;

    const auto &sgroups = getSubstanceGroups(mol_copy);
    TEST_ASSERT(sgroups.size() == 3);

    mol_copy.addAtom();

    TEST_ASSERT(sgroups.size() == 0);
  }
  {
    // replacing an atom will drop SubstanceGroups
    mol_copy = mol;

    const auto &sgroups = getSubstanceGroups(mol_copy);
    TEST_ASSERT(sgroups.size() == 3);

    auto new_atom = Atom();
    mol_copy.replaceAtom(1, &new_atom);

    TEST_ASSERT(sgroups.size() == 0);
  }
  {
    // replacing a new bond will drop SubstanceGroups
    mol_copy = mol;

    const auto &sgroups = getSubstanceGroups(mol_copy);
    TEST_ASSERT(sgroups.size() == 3);

    auto new_bond = Bond(Bond::SINGLE);
    mol_copy.replaceBond(1, &new_bond);

    TEST_ASSERT(sgroups.size() == 0);
  }
  {
    // removing an atom will drop SubstanceGroups
    mol_copy = mol;

    const auto &sgroups = getSubstanceGroups(mol_copy);
    TEST_ASSERT(sgroups.size() == 3);

    mol_copy.removeAtom(1);

    TEST_ASSERT(sgroups.size() == 0);
  }
  {
    // creating a new bond between existing atoms will drop SubstanceGroups
    mol_copy = mol;

    const auto &sgroups = getSubstanceGroups(mol_copy);
    TEST_ASSERT(sgroups.size() == 3);

    mol_copy.addBond(1, 3, Bond::SINGLE);

    TEST_ASSERT(sgroups.size() == 0);
  }
  {
    // removing a bond will drop SubstanceGroups
    mol_copy = mol;

    const auto &sgroups = getSubstanceGroups(mol_copy);
    TEST_ASSERT(sgroups.size() == 3);

    mol_copy.removeBond(1, 2);

    TEST_ASSERT(sgroups.size() == 0);
  }
  {
    // creating a partial bond will drop SubstanceGroups
    mol_copy = mol;

    const auto &sgroups = getSubstanceGroups(mol_copy);
    TEST_ASSERT(sgroups.size() == 3);

    auto *b = mol_copy.createPartialBond(1, Bond::SINGLE);

    TEST_ASSERT(sgroups.size() == 0);
    delete b;
  }
}

void testSubstanceGroupChanges(const std::string &rdbase) {
  BOOST_LOG(rdInfoLog) << " ----------> Test SubstanceGroup property changes"
                       << std::endl;

  std::string fName =
      rdbase +
      "/Code/GraphMol/FileParsers/sgroup_test_data/Sgroups_Data_01.mol";
  std::unique_ptr<ROMol> mol(MolFileToMol(fName));
  TEST_ASSERT(mol);
  auto &sgroups1 = getSubstanceGroups(*mol);
  TEST_ASSERT(sgroups1.size() == 2);

  TEST_ASSERT(sgroups1[0].hasProp("FIELDNAME"));
  TEST_ASSERT(sgroups1[0].getProp<std::string>("FIELDNAME") == "pH");
  sgroups1[0].setProp("FIELDNAME", "pKa");

  const auto &sgroups2 = getSubstanceGroups(*mol);
  TEST_ASSERT(sgroups2.size() == 2);
  TEST_ASSERT(sgroups2[0].hasProp("FIELDNAME"));
  TEST_ASSERT(sgroups2[0].getProp<std::string>("FIELDNAME") == "pKa");
}

int main() {
  std::string rdbase = std::string(getenv("RDBASE"));
  if (rdbase.empty()) {
    std::cerr << "\n\n RDBASE has not been set, aborting.\n\n";
    return 1;
  }

  RDLog::InitLogs();

  testCreateSubstanceGroups();
  testParseSubstanceGroups(rdbase);
  testSubstanceGroupsRoundTrip(rdbase, false);  // test V2000
  testSubstanceGroupsRoundTrip(rdbase, true);   // test V3000
  testPickleSubstanceGroups();
  testModifyMol();
  testSubstanceGroupChanges(rdbase);
  return 0;
}
