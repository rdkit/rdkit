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
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include "GraphMol/FileParsers/FileParsers.h"
#include "GraphMol/FileParsers/MolSupplier.h"
#include "GraphMol/FileParsers/MolWriters.h"

#include <memory>
#include <cstdlib>

using namespace RDKit;

/* Auxiliary functions */
template <class T>
void testIdxVector(const std::vector<T> &groupVector,
                   const std::vector<unsigned int> &reference) {
  size_t vecSize = reference.size();
  TEST_ASSERT(groupVector.size() == vecSize);

  auto sgItr = groupVector.begin();
  for (auto refItr = reference.begin(); refItr != reference.end();
       ++sgItr, ++refItr) {
    TEST_ASSERT(1 + (*sgItr)->getIdx() == *refItr);
  }
}

void testBrackets(
    const std::vector<SGroup::Bracket> &brackets,
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
  auto sg0 = m->getSGroup(0);
  auto sg2 = m->getSGroup(2);
  sg0->setParent(sg2);

  return m;
}

void checkSampleMolecule(RWMol *mol) {
  // Test a molecule created by buildSampleMolecule (or a copy)

  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumSGroups() == 3);

  {
    // First SGroup
    auto sg = mol->getSGroup(0);
    TEST_ASSERT(sg->getType() == "MUL");

    TEST_ASSERT(sg->getProp("SUBTYPE") == "BLO");
    TEST_ASSERT(sg->getProp("MULT") == "n");
    TEST_ASSERT(sg->getProp("CONNECT") == "HH");

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

    TEST_ASSERT(sg->getBondType(bonds[0]) == SGroup::BondType::CBOND);
    TEST_ASSERT(sg->getBondType(bonds[1]) == SGroup::BondType::CBOND);
    TEST_ASSERT(sg->getBondType(bonds[2]) == SGroup::BondType::XBOND);

    TEST_ASSERT(sg->getCompNo() == 7);
    TEST_ASSERT(sg->getProp("ESTATE") == "E");

    std::vector<std::array<std::array<double, 3>, 3>> brackets_reference = {{
        {{{{1., 3., 0.}}, {{5., 7., 0.}}, {{0., 0., 0.}}}},
        {{{{2., 4., 0.}}, {{6., 8., 0.}}, {{0., 0., 0.}}}},
    }};
    testBrackets(sg->getBrackets(), brackets_reference);

    auto cstates = sg->getCStates();
    TEST_ASSERT(cstates.size() == 1);
    TEST_ASSERT(cstates[0].bond == bonds[2]);
    TEST_ASSERT(cstates[0].vector == nullptr);

    TEST_ASSERT(sg->getProp("CLASS") == "TEST CLASS");

    auto ap = sg->getAttachPoints();
    TEST_ASSERT(ap.size() == 1);
    TEST_ASSERT(ap[0].aAtom == atoms[0]);
    TEST_ASSERT(ap[0].lvAtom == atoms[0]);
    TEST_ASSERT(ap[0].id == "XX");

    TEST_ASSERT(sg->getProp("BRKTYP") == "PAREN");

    auto parent = sg->getParent();
    TEST_ASSERT(parent == mol->getSGroup(2));
  }

  {
    // Second SGroup
    auto sg = mol->getSGroup(1);
    TEST_ASSERT(sg->getType() == "SUP");

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
    TEST_ASSERT(sg->getBondType(bonds[0]) == SGroup::BondType::XBOND);
    TEST_ASSERT(sg->getBondType(bonds[1]) == SGroup::BondType::CBOND);
    TEST_ASSERT(sg->getBondType(bonds[2]) == SGroup::BondType::CBOND);

    TEST_ASSERT(sg->getProp("LABEL") == "TEST LABEL");

    auto cstates = sg->getCStates();
    TEST_ASSERT(cstates.size() == 1);
    TEST_ASSERT(cstates[0].bond == bonds[0]);
    TEST_ASSERT(cstates[0].vector != nullptr);
    TEST_ASSERT(cstates[0].vector->x == 3.);
    TEST_ASSERT(cstates[0].vector->y == 4.);
    TEST_ASSERT(cstates[0].vector->z == 0.);

    auto ap = sg->getAttachPoints();
    TEST_ASSERT(ap.size() == 1);
    TEST_ASSERT(ap[0].aAtom == atoms[0]);
    TEST_ASSERT(ap[0].lvAtom == nullptr);
    TEST_ASSERT(ap[0].id == "YY");
  }

  {
    // Third SGroup
    auto sg = mol->getSGroup(2);
    TEST_ASSERT(sg->getType() == "DAT");

    TEST_ASSERT(sg->getProp("FIELDNAME") == "SAMPLE FIELD NAME");
    TEST_ASSERT(sg->getProp("FIELDINFO") == "SAMPLE FIELD INFO");
    TEST_ASSERT(sg->getProp("QUERYTYPE") == "PQ");
    TEST_ASSERT(sg->getProp("QUERYOP") == "SAMPLE QUERY OP");

    TEST_ASSERT(sg->getProp("FIELDDISP") == "SAMPLE FIELD DISP");

    auto dataFields = sg->getDataFields();
    TEST_ASSERT(dataFields.size() == 3);
    TEST_ASSERT(dataFields[0] == "SAMPLE DATA FIELD 1");
    TEST_ASSERT(dataFields[1] == "SAMPLE DATA FIELD 2");
    TEST_ASSERT(dataFields[2] == "SAMPLE DATA FIELD 3");
  }
}

/* End Auxiliary functions */

void testCreateSGroups() {
  BOOST_LOG(rdInfoLog) << " ----------> Testing basic SGroup creation"
                       << std::endl;

  // Create two SGroups and add them to a molecule
  RWMol m;

  SGroup *sg = new SGroup(&m, "DAT");
  m.addSGroup(sg);

  sg = new SGroup(&m, "SUP");
  m.addSGroup(sg);

  sg = nullptr;

  TEST_ASSERT(m.getNumSGroups() == 2);
  TEST_ASSERT(m.getSGroup(0)->getType() == "DAT");
  TEST_ASSERT(m.getSGroup(1)->getType() == "SUP");
}

void testParseSGroups(const std::string &rdbase) {
  BOOST_LOG(rdInfoLog) << " ----------> Parsing Issue3432136_1.mol (V2000)"
                       << std::endl;
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3432136_1.mol";

    auto m = std::unique_ptr<RWMol>(MolFileToMol(fName));

    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumSGroups() == 1);

    auto sgroup = m->getSGroup(0);

    TEST_ASSERT(sgroup->getType() == "MON");

    std::vector<unsigned int> atoms_reference = {2, 3, 4, 1, 5};

    testIdxVector(sgroup->getAtoms(), atoms_reference);

    std::vector<unsigned int> bonds_reference =
        {};  // No bonds defined in this mol
    testIdxVector(sgroup->getBonds(), bonds_reference);

    std::vector<std::array<std::array<double, 3>, 3>> brackets_reference = {{
        {{{{-3.9679, -0.1670, 0.}}, {{-3.9679, 2.1705, 0.}}, {{0., 0., 0.}}}},
        {{{{-0.7244, 2.1705, 0.}}, {{-0.7244, -0.1670, 0.}}, {{0., 0., 0.}}}},
    }};
    testBrackets(sgroup->getBrackets(), brackets_reference);
  }

  BOOST_LOG(rdInfoLog) << " ----------> Parsing Issue3432136_1.v3k.mol (V3000) "
                       << std::endl;
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3432136_1.v3k.mol";
    auto m = std::unique_ptr<RWMol>(MolFileToMol(fName));

    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumSGroups() == 1);

    auto sgroup = m->getSGroup(0);

    TEST_ASSERT(sgroup->getType() == "MON");

    std::vector<unsigned int> atoms_reference = {2, 3, 4, 1, 5};
    testIdxVector(sgroup->getAtoms(), atoms_reference);

    std::vector<unsigned int> bonds_reference =
        {};  // No bonds defined in this mol
    testIdxVector(sgroup->getBonds(), bonds_reference);
  }

  BOOST_LOG(rdInfoLog) << " ----------> Parsing Issue3432136_2.v3k.mol (V3000) "
                       << std::endl;
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3432136_2.v3k.mol";
    auto m = std::unique_ptr<RWMol>(MolFileToMol(fName));

    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumSGroups() == 1);

    auto sgroup = m->getSGroup(0);

    TEST_ASSERT(sgroup->getType() == "SUP");
    TEST_ASSERT(sgroup->getProp("CLASS") == "DEMOCLASS");
    TEST_ASSERT(sgroup->getProp("LABEL") == "abbrev");

    std::vector<unsigned int> atoms_reference = {6, 7, 8, 9, 11, 12};
    testIdxVector(sgroup->getAtoms(), atoms_reference);

    std::vector<unsigned int> bonds_reference = {5};
    testIdxVector(sgroup->getBonds(), bonds_reference);

    auto *bond = sgroup->getBonds()[0];
    TEST_ASSERT(sgroup->getBondType(bond) == SGroup::BondType::XBOND);
  }

  BOOST_LOG(rdInfoLog) << " ----------> Parsing Issue3432136_2.mol (V2000) "
                       << std::endl;
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3432136_2.mol";
    auto m = std::unique_ptr<RWMol>(MolFileToMol(fName));

    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumSGroups() == 1);

    auto sgroup = m->getSGroup(0);

    TEST_ASSERT(sgroup->getType() == "SUP");

    std::vector<unsigned int> atoms_reference = {6, 7, 8, 9, 11, 12};
    testIdxVector(sgroup->getAtoms(), atoms_reference);

    std::vector<unsigned int> bonds_reference = {5};
    testIdxVector(sgroup->getBonds(), bonds_reference);

    auto *bond = sgroup->getBonds()[0];
    TEST_ASSERT(sgroup->getBondType(bond) == SGroup::BondType::XBOND);
  }
}

void testSGroupsRoundTrip(const std::string &rdbase, bool forceV3000) {
  BOOST_LOG(rdInfoLog)
      << " ----------> Testing SGroup writing & parsing Roundtrip ("
      << (forceV3000 ? "V3000" : "V2000") << ')' << std::endl;

  std::string fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/testSGroupsSample_" +
      (forceV3000 ? "V3000" : "V2000") + ".mol";

  {
    auto sampleMol = buildSampleMolecule();

    TEST_ASSERT(sampleMol->getNumSGroups() == 3);

    auto writer = SDWriter(fName);
    writer.setForceV3000(forceV3000);
    writer.write(*sampleMol);
    writer.close();
  }
  auto roundtripMol = std::shared_ptr<RWMol>(MolFileToMol(fName));
  checkSampleMolecule(roundtripMol.get());
}

void testPickleSGroups() {
  BOOST_LOG(rdInfoLog)
      << " ----------> Testing SGroup pickling & unpickling Roundtrip"
      << std::endl;

  std::string pkl;

  {
    auto sampleMol = buildSampleMolecule();

    MolPickler::pickleMol(*sampleMol, pkl);
  }

  auto roundtripMol = std::shared_ptr<RWMol>(new RWMol(pkl));
  checkSampleMolecule(roundtripMol.get());
}

void testModifyMol() {
  BOOST_LOG(rdInfoLog) << " ----------> Test dropping SGroups on modification"
                       << std::endl;

  auto mol = *(buildSampleMolecule());
  auto mol_copy = mol;

  {  // insertion will drop SGroups
    TEST_ASSERT(mol_copy.getNumSGroups() == 3);
    mol_copy.insertMol(mol);
    TEST_ASSERT(mol_copy.getNumSGroups() == 0);
  }
  {
    // adding an atom will drop SGroups
    mol_copy = mol;
    TEST_ASSERT(mol_copy.getNumSGroups() == 3);
    mol_copy.addAtom();
    TEST_ASSERT(mol_copy.getNumSGroups() == 0);
  }
  {
    // replacing an atom will drop SGroups
    mol_copy = mol;
    auto new_atom = Atom();
    TEST_ASSERT(mol_copy.getNumSGroups() == 3);
    mol_copy.replaceAtom(1, &new_atom);
    TEST_ASSERT(mol_copy.getNumSGroups() == 0);
  }
  {
    // replacing a new bond will drop SGroups
    mol_copy = mol;
    auto new_bond = Bond(Bond::SINGLE);
    TEST_ASSERT(mol_copy.getNumSGroups() == 3);
    mol_copy.replaceBond(1, &new_bond);
    TEST_ASSERT(mol_copy.getNumSGroups() == 0);
  }
  {
    // removing an atom will drop SGroups
    mol_copy = mol;
    TEST_ASSERT(mol_copy.getNumSGroups() == 3);
    mol_copy.removeAtom(1);
    TEST_ASSERT(mol_copy.getNumSGroups() == 0);
  }
  {
    // creating a new bond between existing atoms will drop SGroups
    mol_copy = mol;
    TEST_ASSERT(mol_copy.getNumSGroups() == 3);
    mol_copy.addBond(1, 3, Bond::SINGLE);
    TEST_ASSERT(mol_copy.getNumSGroups() == 0);
  }
  {
    // removing a bond will drop SGroups
    mol_copy = mol;
    TEST_ASSERT(mol_copy.getNumSGroups() == 3);
    mol_copy.removeBond(1, 2);
    TEST_ASSERT(mol_copy.getNumSGroups() == 0);
  }
  {
    // creating a partial bond will drop SGroups
    mol_copy = mol;
    TEST_ASSERT(mol_copy.getNumSGroups() == 3);
    mol_copy.createPartialBond(1, Bond::SINGLE);
    TEST_ASSERT(mol_copy.getNumSGroups() == 0);
  }
}

int main() {
  std::string rdbase = std::string(getenv("RDBASE"));
  if (rdbase.empty()) {
    std::cerr << "\n\n RDBASE has not been set, aborting.\n\n";
    return 1;
  }

  RDLog::InitLogs();

  testCreateSGroups();
  testParseSGroups(rdbase);
  testSGroupsRoundTrip(rdbase, false);  // test V2000
  testSGroupsRoundTrip(rdbase, true);   // test V3000
  testPickleSGroups();
  testModifyMol();

  return 0;
}
