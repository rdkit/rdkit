//
//
//  Copyright (C) 2018-2021 Greg Landrum and other RDKit contributors
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
#include <GraphMol/MolOps.h>
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
    sg.setProp("index", 1u);

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
    sg.setProp("index", 2u);

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
    sg.setProp("index", 3u);

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

  // We have to set a parent with a lower index in V2000 mol blocks:
  const auto &sgroups = getSubstanceGroups(mol);
  sgroups.at(1).setProp<unsigned int>("PARENT", 1u);

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

    std::vector<std::array<std::array<double, 3>, 3>> brackets_reference = {
        {{{{1., 3., 0.}}, {{5., 7., 0.}}, {{0., 0., 0.}}}},
        {{{{2., 4., 0.}}, {{6., 8., 0.}}, {{0., 0., 0.}}}},
    };
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
    TEST_ASSERT(sg.getProp<unsigned int>("PARENT") == 1u);
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

    std::vector<std::array<std::array<double, 3>, 3>> brackets_reference = {
        {{{{-3.9679, -0.1670, 0.}}, {{-3.9679, 2.1705, 0.}}, {{0., 0., 0.}}}},
        {{{{-0.7244, 2.1705, 0.}}, {{-0.7244, -0.1670, 0.}}, {{0., 0., 0.}}}},
    };
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

  // make sure that calling clear() on the molecule actually clears out
  // SubstanceGroups
  //  This was GitHub #3167
  {
    auto tmol = mol;
    TEST_ASSERT(getSubstanceGroups(tmol).size() == 3);
    tmol.clear();
    TEST_ASSERT(getSubstanceGroups(tmol).size() == 0);
  }

  auto mol_copy = mol;

  const auto &sgroups = getSubstanceGroups(mol);
  TEST_ASSERT(sgroups.size() == 3);

  {  // insertion does not affect SubstanceGroups
    const auto &sgroups = getSubstanceGroups(mol_copy);
    TEST_ASSERT(sgroups.size() == 3);

    mol_copy.insertMol(mol);

    TEST_ASSERT(sgroups.size() == 3);
  }
  {
    // adding an atom does not affect SubstanceGroups
    mol_copy = mol;

    const auto &sgroups = getSubstanceGroups(mol_copy);
    TEST_ASSERT(sgroups.size() == 3);

    mol_copy.addAtom();

    TEST_ASSERT(sgroups.size() == 3);
  }
  {
    // replacing an atom does not drop SubstanceGroups that include that atom
    mol_copy = mol;

    const auto &sgroups = getSubstanceGroups(mol_copy);
    TEST_ASSERT(sgroups.size() == 3);

    auto new_atom = Atom();
    mol_copy.replaceAtom(1, &new_atom);

    TEST_ASSERT(sgroups.size() == 3);
  }
  {
    // replacing a bond does not drop SubstanceGroups that include that bond
    mol_copy = mol;

    const auto &sgroups = getSubstanceGroups(mol_copy);
    TEST_ASSERT(sgroups.size() == 3);

    auto new_bond = Bond(Bond::SINGLE);
    mol_copy.replaceBond(1, &new_bond);

    TEST_ASSERT(sgroups.size() == 3);
  }
  {
    // removing an atom will drop SubstanceGroups that include that atom
    mol_copy = mol;

    const auto &sgroups = getSubstanceGroups(mol_copy);
    TEST_ASSERT(sgroups.size() == 3);

    mol_copy.removeAtom(1);

    TEST_ASSERT(sgroups.size() == 1);
  }
  {
    // creating a new bond between existing atoms does not affect
    // SubstanceGroups
    mol_copy = mol;

    const auto &sgroups = getSubstanceGroups(mol_copy);
    TEST_ASSERT(sgroups.size() == 3);

    mol_copy.addBond(1, 3, Bond::SINGLE);

    TEST_ASSERT(sgroups.size() == 3);
  }
  {
    // removing a bond will drop SubstanceGroups that involve that bond
    mol_copy = mol;

    const auto &sgroups = getSubstanceGroups(mol_copy);
    TEST_ASSERT(sgroups.size() == 3);

    mol_copy.removeBond(1, 2);

    TEST_ASSERT(sgroups.size() == 1);
  }
  {
    // creating a partial bond does not effect SubstanceGroups
    mol_copy = mol;

    const auto &sgroups = getSubstanceGroups(mol_copy);
    TEST_ASSERT(sgroups.size() == 3);

    auto *b = mol_copy.createPartialBond(1, Bond::SINGLE);

    TEST_ASSERT(sgroups.size() == 3);
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

void testSubstanceGroupsAndRemoveAtoms(const std::string &rdbase) {
  BOOST_LOG(rdInfoLog)
      << " ----------> Test impact of removeAtom on SubstanceGroups"
      << std::endl;
  {
    std::string fName =
        rdbase + "/Code/GraphMol/test_data/sgroups_and_remove_atoms_1.mol";
    std::unique_ptr<RWMol> mol(MolFileToMol(fName));
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 13);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.size() == 1);
      TEST_ASSERT(sgroups[0].getAtoms().size() == 3);
      std::vector<unsigned int> tgt{10, 11, 12};
      TEST_ASSERT(sgroups[0].getAtoms() == tgt);
      TEST_ASSERT(sgroups[0].getBonds().size() == 1);
      tgt = {9};
      TEST_ASSERT(sgroups[0].getBonds() == tgt);
      auto aps = sgroups[0].getAttachPoints();
      TEST_ASSERT(aps.size() == 1);
      TEST_ASSERT(aps[0].aIdx == 11);
      TEST_ASSERT(aps[0].lvIdx == 3);
    }
    // remove an atom that's not in an S-group
    mol->removeAtom(9);
    TEST_ASSERT(mol->getNumAtoms() == 12);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.size() == 1);
      TEST_ASSERT(sgroups[0].getAtoms().size() == 3);
      std::vector<unsigned int> tgt{9, 10, 11};
      TEST_ASSERT(sgroups[0].getAtoms() == tgt);
      TEST_ASSERT(sgroups[0].getBonds().size() == 1);
      tgt = {8};
      TEST_ASSERT(sgroups[0].getBonds() == tgt);
      auto aps = sgroups[0].getAttachPoints();
      TEST_ASSERT(aps.size() == 1);
      TEST_ASSERT(aps[0].aIdx == 10);
      TEST_ASSERT(aps[0].lvIdx == 3);
    }
    // remove an atom that is in an S-group
    mol->removeAtom(10);
    TEST_ASSERT(mol->getNumAtoms() == 11);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.empty());
    }
  }
  {
    // example with hs to be removed
    std::string fName =
        rdbase + "/Code/GraphMol/test_data/sgroups_and_remove_atoms_2.mol";
    std::unique_ptr<RWMol> mol(MolFileToMol(fName));
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 14);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.size() == 2);
      TEST_ASSERT(sgroups[0].getAtoms().size() == 3);
      std::vector<unsigned int> tgt{9, 11, 12};
      TEST_ASSERT(sgroups[0].getAtoms() == tgt);
      TEST_ASSERT(sgroups[0].getBonds().size() == 1);
      TEST_ASSERT(sgroups[0].getBonds()[0] == 9);

      TEST_ASSERT(sgroups[1].getAtoms().size() == 2);
      tgt = {10, 13};
      TEST_ASSERT(sgroups[1].getAtoms() == tgt);
      TEST_ASSERT(sgroups[1].getBonds().size() == 1);
      TEST_ASSERT(sgroups[1].getBonds()[0] == 10);
    }
    // remove an atom in the first S group, make sure the second one survives
    mol->removeAtom(11);
    TEST_ASSERT(mol->getNumAtoms() == 13);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.size() == 1);
      TEST_ASSERT(sgroups[0].getAtoms().size() == 2);
      std::vector<unsigned int> tgt{10, 12};
      TEST_ASSERT(sgroups[0].getAtoms() == tgt);
      TEST_ASSERT(sgroups[0].getBonds().size() == 1);
      TEST_ASSERT(sgroups[0].getBonds()[0] == 9);
    }
  }
  {  // example with CSTATE
    std::string fName =
        rdbase + "/Code/GraphMol/test_data/sgroups_and_remove_atoms_3.mol";
    std::unique_ptr<RWMol> mol(MolFileToMol(fName));
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 14);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.size() == 1);
      TEST_ASSERT(sgroups[0].getAtoms().size() == 7);
      std::vector<unsigned int> tgt{7, 8, 9, 10, 11, 12, 13};
      TEST_ASSERT(sgroups[0].getAtoms() == tgt);
      TEST_ASSERT(sgroups[0].getBonds().size() == 1);
      tgt = {8};
      TEST_ASSERT(sgroups[0].getBonds() == tgt);
      auto cstates = sgroups[0].getCStates();
      TEST_ASSERT(cstates.size() == 1);
      TEST_ASSERT(cstates[0].bondIdx == 8);
    }
    // remove an atom that's not in an S-group
    mol->removeAtom(1);
    TEST_ASSERT(mol->getNumAtoms() == 13);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.size() == 1);
      TEST_ASSERT(sgroups[0].getAtoms().size() == 7);
      std::vector<unsigned int> tgt{6, 7, 8, 9, 10, 11, 12};
      TEST_ASSERT(sgroups[0].getAtoms() == tgt);
      TEST_ASSERT(sgroups[0].getBonds().size() == 1);
      tgt = {6};
      TEST_ASSERT(sgroups[0].getBonds() == tgt);
      auto cstates = sgroups[0].getCStates();
      TEST_ASSERT(cstates.size() == 1);
      TEST_ASSERT(cstates[0].bondIdx == 6);
    }
    // remove an atom that is in an S-group
    mol->removeAtom(10);
    TEST_ASSERT(mol->getNumAtoms() == 12);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.empty());
    }
  }
  {  // example with PATOMS
    std::string fName =
        rdbase + "/Code/GraphMol/test_data/sgroups_and_remove_atoms_4.mol";
    std::unique_ptr<RWMol> mol(MolFileToMol(fName));
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 9);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.size() == 1);
      TEST_ASSERT(sgroups[0].getAtoms().size() == 6);
      std::vector<unsigned int> tgt{2, 3, 5, 6, 7, 8};
      TEST_ASSERT(sgroups[0].getAtoms() == tgt);
      TEST_ASSERT(sgroups[0].getBonds().size() == 2);
      tgt = {7, 6};
      TEST_ASSERT(sgroups[0].getBonds() == tgt);
      TEST_ASSERT(sgroups[0].getParentAtoms().size() == 3);
      tgt = {2, 3, 5};
      TEST_ASSERT(sgroups[0].getParentAtoms() == tgt);
    }
    // remove an atom that's not in an S-group
    mol->removeAtom(0u);
    TEST_ASSERT(mol->getNumAtoms() == 8);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.size() == 1);
      TEST_ASSERT(sgroups[0].getAtoms().size() == 6);
      std::vector<unsigned int> tgt{1, 2, 4, 5, 6, 7};
      TEST_ASSERT(sgroups[0].getAtoms() == tgt);
      TEST_ASSERT(sgroups[0].getBonds().size() == 2);
      tgt = {6, 5};
      TEST_ASSERT(sgroups[0].getBonds() == tgt);
      TEST_ASSERT(sgroups[0].getParentAtoms().size() == 3);
      tgt = {1, 2, 4};
      TEST_ASSERT(sgroups[0].getParentAtoms() == tgt);
    }
    // remove an atom that is in an S-group
    mol->removeAtom(7);
    TEST_ASSERT(mol->getNumAtoms() == 7);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.empty());
    }
  }

  {  // example with parent
    std::string fName =
        rdbase + "/Code/GraphMol/test_data/sgroups_and_remove_atoms_5.mol";
    std::unique_ptr<RWMol> mol(MolFileToMol(fName));
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 18);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.size() == 3);
      TEST_ASSERT(sgroups[0].hasProp("index"))
      TEST_ASSERT(sgroups[0].getProp<unsigned int>("index") == 1);
      TEST_ASSERT(sgroups[0].getAtoms().size() == 3);
      std::vector<unsigned int> tgt{3, 2, 7};
      TEST_ASSERT(sgroups[0].getAtoms() == tgt);
      TEST_ASSERT(sgroups[0].getBonds().size() == 2);
      tgt = {1, 8};
      TEST_ASSERT(sgroups[0].getBonds() == tgt);

      TEST_ASSERT(sgroups[1].hasProp("index"))
      TEST_ASSERT(sgroups[1].getProp<unsigned int>("index") == 2);
      TEST_ASSERT(sgroups[1].getAtoms().size() == 6);
      tgt = {5, 4, 10, 15, 16, 17};
      TEST_ASSERT(sgroups[1].getAtoms() == tgt);
      TEST_ASSERT(sgroups[1].getBonds().size() == 2);
      tgt = {8, 16};
      TEST_ASSERT(sgroups[1].getBonds() == tgt);
      TEST_ASSERT(sgroups[1].hasProp("PARENT"))
      TEST_ASSERT(sgroups[1].getProp<unsigned int>("PARENT") == 1);
    }
    // remove an atom that's not in an S-group
    mol->removeAtom(0u);
    TEST_ASSERT(mol->getNumAtoms() == 17);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.size() == 3);
      TEST_ASSERT(sgroups[0].getAtoms().size() == 3);
      std::vector<unsigned int> tgt{2, 1, 7};
    }
    // remove an atom that is in a parent S-group
    mol->removeAtom(1);
    TEST_ASSERT(mol->getNumAtoms() == 16);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.size() == 1);
      TEST_ASSERT(sgroups[0].getAtoms().size() == 4);
      TEST_ASSERT(sgroups[0].hasProp("index"))
      TEST_ASSERT(sgroups[0].getProp<unsigned int>("index") == 3);
    }
  }

  {  // example with sgroup hierarchy
    std::string fName =
        rdbase + "/Code/GraphMol/test_data/sgroups_and_remove_atoms_6.mol";
    std::unique_ptr<RWMol> mol(MolFileToMol(fName));
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 13);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.size() == 3);
      TEST_ASSERT(sgroups[0].hasProp("index"))
      TEST_ASSERT(sgroups[0].getProp<unsigned int>("index") == 1);
      TEST_ASSERT(sgroups[0].getAtoms().size() == 3);
      std::vector<unsigned int> tgt{3, 2, 7};
      TEST_ASSERT(sgroups[0].getAtoms() == tgt);
      TEST_ASSERT(sgroups[0].getBonds().size() == 2);
      tgt = {1, 8};
      TEST_ASSERT(sgroups[0].getBonds() == tgt);

      TEST_ASSERT(sgroups[1].hasProp("index"))
      TEST_ASSERT(sgroups[1].getProp<unsigned int>("index") == 2);
      TEST_ASSERT(sgroups[1].getAtoms().size() == 3);
      tgt = {5, 4, 10};
      TEST_ASSERT(sgroups[1].getAtoms() == tgt);
      TEST_ASSERT(sgroups[1].getBonds().size() == 2);
      tgt = {8, 11};
      TEST_ASSERT(sgroups[1].getBonds() == tgt);
      TEST_ASSERT(sgroups[1].hasProp("PARENT"))
      TEST_ASSERT(sgroups[1].getProp<unsigned int>("PARENT") == 1);

      TEST_ASSERT(sgroups[2].hasProp("index"))
      TEST_ASSERT(sgroups[2].getProp<unsigned int>("index") == 3);
      TEST_ASSERT(sgroups[2].getAtoms().size() == 2);
      tgt = {9, 8};
      TEST_ASSERT(sgroups[2].getAtoms() == tgt);
      TEST_ASSERT(sgroups[2].getBonds().size() == 2);
      tgt = {5, 10};
      TEST_ASSERT(sgroups[2].getBonds() == tgt);
      TEST_ASSERT(sgroups[2].hasProp("PARENT"))
      TEST_ASSERT(sgroups[2].getProp<unsigned int>("PARENT") == 2);
    }
    // remove an atom that's not in an S-group
    mol->removeAtom(0u);
    TEST_ASSERT(mol->getNumAtoms() == 12);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.size() == 3);
      TEST_ASSERT(sgroups[0].getAtoms().size() == 3);
      std::vector<unsigned int> tgt{2, 1, 6};
    }
    // remove an atom that is in an S-group at the top of the hierarchy
    mol->removeAtom(1);
    TEST_ASSERT(mol->getNumAtoms() == 11);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.size() == 0);
    }
  }

  {  // example with things in odd order and large id values
     // NOTE that biovia draw doesn't parse this file properly
    std::string fName =
        rdbase + "/Code/GraphMol/test_data/sgroups_and_remove_atoms_7.mol";
    std::unique_ptr<RWMol> mol(MolFileToMol(fName));
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 18);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.size() == 3);
      TEST_ASSERT(sgroups[2].hasProp("index"))
      TEST_ASSERT(sgroups[2].getProp<unsigned int>("index") == 20);
      TEST_ASSERT(sgroups[2].getAtoms().size() == 6);
      std::vector<unsigned int> tgt{5, 4, 10, 15, 16, 17};
      TEST_ASSERT(sgroups[2].getAtoms() == tgt);
      TEST_ASSERT(sgroups[2].getBonds().size() == 2);
      tgt = {8, 16};
      TEST_ASSERT(sgroups[2].getBonds() == tgt);
      TEST_ASSERT(sgroups[2].getParentAtoms().size() == 3);
      tgt = {5, 4, 10};
      TEST_ASSERT(sgroups[2].getParentAtoms() == tgt);
      TEST_ASSERT(sgroups[2].hasProp("PARENT"))
      TEST_ASSERT(sgroups[2].getProp<unsigned int>("PARENT") == 10);

      TEST_ASSERT(sgroups[1].hasProp("index"))
      TEST_ASSERT(sgroups[1].getProp<unsigned int>("index") == 10);
      TEST_ASSERT(sgroups[1].getAtoms().size() == 3);
      tgt = {3, 2, 7};
      TEST_ASSERT(sgroups[1].getAtoms() == tgt);
      TEST_ASSERT(sgroups[1].getBonds().size() == 2);
      tgt = {1, 8};
      TEST_ASSERT(sgroups[1].getBonds() == tgt);
      TEST_ASSERT(sgroups[1].getParentAtoms().size() == 0);
    }
  }
  {  // copolymer example with PARENT
    std::string fName =
        rdbase + "/Code/GraphMol/test_data/sgroups_copolymer.mol";
    std::unique_ptr<RWMol> mol(MolFileToMol(fName));
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 9);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.size() == 3);
      TEST_ASSERT(sgroups[2].hasProp("index"))
      TEST_ASSERT(sgroups[2].getProp<unsigned int>("index") == 10);
      TEST_ASSERT(sgroups[2].getAtoms().size() == 5);
      std::vector<unsigned int> tgt{3, 2, 4, 5, 7};
      TEST_ASSERT(sgroups[2].getAtoms() == tgt);
      TEST_ASSERT(sgroups[2].getBonds().size() == 2);
      tgt = {1, 5};
      TEST_ASSERT(sgroups[2].getBonds() == tgt);
      TEST_ASSERT(!sgroups[2].hasProp("PARENT"))

      TEST_ASSERT(sgroups[0].hasProp("index"))
      TEST_ASSERT(sgroups[0].getProp<unsigned int>("index") == 2);
      TEST_ASSERT(sgroups[0].getAtoms().size() == 2);
      tgt = {3, 2};
      TEST_ASSERT(sgroups[0].getAtoms() == tgt);
      TEST_ASSERT(sgroups[0].getBonds().size() == 2);
      tgt = {1, 3};
      TEST_ASSERT(sgroups[0].getBonds() == tgt);
      TEST_ASSERT(sgroups[0].hasProp("PARENT"))
      TEST_ASSERT(sgroups[0].getProp<unsigned int>("PARENT") == 10);
    }
    // remove an atom that's not in an S-group
    mol->removeAtom(0u);
    TEST_ASSERT(mol->getNumAtoms() == 8);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.size() == 3);
      TEST_ASSERT(sgroups[2].hasProp("index"))
      TEST_ASSERT(sgroups[2].getProp<unsigned int>("index") == 10);
      TEST_ASSERT(sgroups[2].getAtoms().size() == 5);
      std::vector<unsigned int> tgt{2, 1, 3, 4, 6};
      TEST_ASSERT(sgroups[2].getAtoms() == tgt);
      TEST_ASSERT(sgroups[2].getBonds().size() == 2);
      tgt = {0, 4};
      TEST_ASSERT(sgroups[2].getBonds() == tgt);
      TEST_ASSERT(!sgroups[2].hasProp("PARENT"))

      TEST_ASSERT(sgroups[0].hasProp("index"))
      TEST_ASSERT(sgroups[0].getProp<unsigned int>("index") == 2);
      TEST_ASSERT(sgroups[0].getAtoms().size() == 2);
      tgt = {2, 1};
      TEST_ASSERT(sgroups[0].getAtoms() == tgt);
      TEST_ASSERT(sgroups[0].getBonds().size() == 2);
      tgt = {0, 2};
      TEST_ASSERT(sgroups[0].getBonds() == tgt);
      TEST_ASSERT(sgroups[0].hasProp("PARENT"))
      TEST_ASSERT(sgroups[0].getProp<unsigned int>("PARENT") == 10);
    }
    // remove an atom from parent, make sure children also get deleted
    mol->removeAtom(1u);
    TEST_ASSERT(mol->getNumAtoms() == 7);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.size() == 0);
    }
  }
  {  // copolymer example 2 with PARENT, same as the previous but with a
     // different ordering in the input file
    std::string fName =
        rdbase + "/Code/GraphMol/test_data/sgroups_copolymer2.mol";
    std::unique_ptr<RWMol> mol(MolFileToMol(fName));
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 9);
    {
      auto &sgroups = getSubstanceGroups(*mol);
      TEST_ASSERT(sgroups.size() == 3);
      TEST_ASSERT(sgroups[2].hasProp("index"))
      TEST_ASSERT(sgroups[2].getProp<unsigned int>("index") == 10);
      TEST_ASSERT(sgroups[2].getAtoms().size() == 5);
      std::vector<unsigned int> tgt{3, 2, 4, 5, 7};
      TEST_ASSERT(sgroups[2].getAtoms() == tgt);
      TEST_ASSERT(sgroups[2].getBonds().size() == 2);
      tgt = {1, 5};
      TEST_ASSERT(sgroups[2].getBonds() == tgt);
      TEST_ASSERT(!sgroups[2].hasProp("PARENT"))

      TEST_ASSERT(sgroups[0].hasProp("index"))
      TEST_ASSERT(sgroups[0].getProp<unsigned int>("index") == 2);
      TEST_ASSERT(sgroups[0].getAtoms().size() == 2);
      tgt = {3, 2};
      TEST_ASSERT(sgroups[0].getAtoms() == tgt);
      TEST_ASSERT(sgroups[0].getBonds().size() == 2);
      tgt = {1, 3};
      TEST_ASSERT(sgroups[0].getBonds() == tgt);
      TEST_ASSERT(sgroups[0].hasProp("PARENT"))
      TEST_ASSERT(sgroups[0].getProp<unsigned int>("PARENT") == 10);
    }
  }
}

void testSubstanceGroupsAndRemoveHs(const std::string &rdbase) {
  BOOST_LOG(rdInfoLog) << " ----------> Test impact of SubstanceGroups on "
                          "removeHs (GitHub #3169)"
                       << std::endl;
  {
    std::string fName =
        rdbase + "/Code/GraphMol/test_data/sgroups_and_remove_Hs_1.mol";
    bool sanitize = true;
    bool removeHs = false;
    std::unique_ptr<RWMol> mol(MolFileToMol(fName, sanitize, removeHs));
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 8);
    TEST_ASSERT(getSubstanceGroups(*mol).size() == 1);

    {
      RWMol mol_copy = *mol;
      MolOps::RemoveHsParameters ps;
      ps.removeInSGroups = true;
      MolOps::removeHs(mol_copy, ps);
      TEST_ASSERT(mol_copy.getNumAtoms() == 6);
      TEST_ASSERT(getSubstanceGroups(mol_copy).size() == 1);
    }
    {  // check removeAllHs() too
      RWMol mol_copy = *mol;
      MolOps::removeAllHs(mol_copy);
      TEST_ASSERT(mol_copy.getNumAtoms() == 6);
      TEST_ASSERT(getSubstanceGroups(mol_copy).size() == 1);
    }
  }
}

void testKeepSpecialHsOnRemoval() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing github 5399: "
                       << "Allow removal of non-special H atoms in SGroups"
                       << std::endl;

  // Keep XBOND Hydrogens
  {
    auto mol = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -3.914286 2.714286 0.000000 0
M  V30 2 H -2.914286 2.714286 0.000000 0
M  V30 3 H -4.247615 1.771475 0.000000 0
M  V30 4 H -4.491638 3.530781 0.000000 0
M  V30 5 H -4.904803 2.576897 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(3 1 3 4) XBONDS=(2 1 4) CONNECT=HT
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;

    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 3);

    auto sgs = getSubstanceGroups(*mol);
    TEST_ASSERT(sgs.size() == 1);
    std::vector<unsigned> ref_atoms{0};
    TEST_ASSERT(sgs[0].getAtoms() == ref_atoms);
  }

  // Keep SAP aIdx (but remove lvIdx).
  // (This molecule is completely bogus).
  {
    auto molblock = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 2 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.8551 -2.605 0 0
M  V30 2 C -0.5214 -3.375 0 0
M  V30 3 H -0.5214 -4.9149 0 0
M  V30 4 H 0.8123 -5.6849 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 0 ATOMS=(2 3 4) SAP=(3 3 4 1) XBONDS=(1 2) LABEL="bogus"
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB";

    bool sanitize = false;  // One H has a weird valence
    bool removeHs = false;
    std::unique_ptr<RWMol> mol(MolBlockToMol(molblock, sanitize, removeHs));
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 4);

    MolOps::removeAllHs(*mol);
    TEST_ASSERT(mol->getNumAtoms() == 3);

    auto sgs = getSubstanceGroups(*mol);
    TEST_ASSERT(sgs.size() == 1);
    std::vector<unsigned> ref_atoms{2};
    TEST_ASSERT(sgs[0].getAtoms() == ref_atoms);

    auto saps = sgs[0].getAttachPoints();
    TEST_ASSERT(saps.size() == 1);
    TEST_ASSERT(saps[0].lvIdx == -1);
  }

  // H atom part of a CState bond
  {
    auto mol = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 1 0 0
M  V30 BEGIN ATOM
M  V30 1 H 12.0001 -4.402 0 0
M  V30 2 C 10.4601 -4.402 0 0
M  V30 3 H 8.1485 -3.075 0 0
M  V30 4 O 9.6883 -3.0729 0 0
M  V30 5 O 9.6919 -5.725 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 4 3
M  V30 3 2 2 5
M  V30 4 1 4 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 0 ATOMS=(4 2 3 4 5) SAP=(3 1 2 1) XBONDS=(1 1) CSTATE=(4 1 0 0.82 0) LABEL=Boc
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;

    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 4);

    auto sgs = getSubstanceGroups(*mol);
    TEST_ASSERT(sgs.size() == 1);

    std::vector<unsigned> ref_atoms{{1, 2, 3}};
    TEST_ASSERT(sgs[0].getAtoms() == ref_atoms);

    auto csts = sgs[0].getCStates();
    TEST_ASSERT(csts.size() == 1);
    TEST_ASSERT(csts[0].bondIdx == 0);
  }

  // SGroups with only H atoms
  {
    auto mol = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 2 0 0
M  V30 BEGIN ATOM
M  V30 1 C -3.514286 2.200000 0.000000 0
M  V30 2 H -2.514286 2.200000 0.000000 0
M  V30 3 H -3.847615 1.257190 0.000000 0
M  V30 4 H -4.091638 3.016495 0.000000 0
M  V30 5 H -4.504803 2.062612 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 4) FIELDNAME="dummy1" FIELDDATA="H1"
M  V30 2 DAT 0 ATOMS=(1 5) FIELDNAME="dummy2" FIELDDATA="H2"
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;

    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 3);

    auto sgs = getSubstanceGroups(*mol);
    TEST_ASSERT(sgs.size() == 2);

    for (unsigned i = 0; i < 2; ++i) {
      std::vector<unsigned> ref_atoms{i + 1};
      TEST_ASSERT(sgs[i].getAtoms() == ref_atoms);
    }
  }
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
  testSubstanceGroupsAndRemoveAtoms(rdbase);
  testSubstanceGroupsAndRemoveHs(rdbase);
  testKeepSpecialHsOnRemoval();

  return 0;
}
