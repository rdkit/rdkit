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
  sgroups.at(0).setProp<unsigned int>("PARENT", 3);

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

    std::vector<std::array<std::array<double, 3>, 3>> brackets_reference = {
        {{{{1., 3., 0.}}, {{5., 7., 0.}}, {{0., 0., 0.}}}},
        {{{{2., 4., 0.}}, {{6., 8., 0.}}, {{0., 0., 0.}}}},
    };
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

    CHECK(sg.getProp<unsigned int>("PARENT") == 3u);
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

TEST_CASE("Removing sgroups", "[Sgroups]") {
  SECTION("basics") {
    auto m1 = R"CTAB(
  Mrv2014 07312005252D          
 
  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 6 3 0 0
M  V30 BEGIN ATOM
M  V30 1 * -12.75 11.5 0 0
M  V30 2 O -11.4163 12.27 0 0
M  V30 3 C -10.0826 11.5 0 0
M  V30 4 C -8.749 12.27 0 0
M  V30 5 O -10.0826 9.96 0 0
M  V30 6 N -7.4153 11.5 0 0
M  V30 7 C -6.0816 12.27 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 2 3 5
M  V30 5 1 4 6
M  V30 6 1 6 7
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(3 2 3 5) XBONDS=(2 1 3) BRKXYZ=(9 -9.9955 12.6173 0 -
M  V30 -9.0715 11.0169 0 0 0 0) BRKXYZ=(9 -11.5035 11.1527 0 -12.4275 12.7531 -
M  V30 0 0 0 0) CONNECT=HT LABEL=n
M  V30 2 DAT 0 ATOMS=(1 6) FIELDNAME=foo_data -
M  V30 FIELDDISP="   -7.4153   11.5000    DAU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 FIELDDATA=bar
M  V30 3 DAT 0 ATOMS=(1 7) FIELDNAME=bar_data -
M  V30 FIELDDISP="   -6.0816   12.2700    DAU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 FIELDDATA=baz
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m1);
    auto &sgs = getSubstanceGroups(*m1);
    CHECK(sgs.size() == 3);
    sgs.erase(++sgs.begin());
    CHECK(sgs.size() == 2);
    CHECK(getSubstanceGroups(*m1).size() == 2);
    auto molb = MolToV3KMolBlock(*m1);
    CHECK(molb.find("foo_data") == std::string::npos);
    CHECK(molb.find("M  V30 2 DAT 0 ATOMS=(1 7) FIELDNAME=bar_data") !=
          std::string::npos);
  }
}

TEST_CASE(
    "Github #3315: SubstanceGroups should not be written with quotes around "
    "missing fields",
    "[Sgroups][bug]") {
  SECTION("basics") {
    auto m1 = R"CTAB(
  Mrv2014 07312012022D          

  2  1  0  0  0  0            999 V2000
    1.4295    0.1449    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4266   -0.6801    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
M  STY  1   1 DAT
M  SAL   1  1   1
M  SDT   1 test sgroup                                           
M  SDD   1     0.5348   -0.3403    DRU   ALL  0       0  
M  SED   1 sgroupval
M  END
)CTAB"_ctab;
    REQUIRE(m1);
    auto molb = MolToV3KMolBlock(*m1);
    CHECK(molb.find("FIELDINFO") == std::string::npos);
    CHECK(molb.find("QUERYTYPE") == std::string::npos);
    CHECK(molb.find("QUERYOP") == std::string::npos);
    CHECK(molb.find("FIELDNAME") != std::string::npos);
  }
}

TEST_CASE("Allow brackets, cstates, and attachment points to be removed",
          "[Sgroups]") {
  auto m1 = R"CTAB(example
 -ISIS-  10171405052D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 14 15 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 6.4292 -1.1916 0 0 CFG=3
M  V30 2 C 7.0125 -0.6042 0 0
M  V30 3 N 6.4292 -0.0250999 0 0
M  V30 4 C 5.8416 -0.6042 0 0
M  V30 5 C 5.8416 -1.7708 0 0
M  V30 6 N 6.4292 -2.3584 0 0 CFG=3
M  V30 7 C 7.0125 -1.7708 0 0
M  V30 8 O 5.7166 -3.5875 0 0
M  V30 9 C 5.7166 -4.4125 0 0 CFG=3
M  V30 10 C 4.8875 -4.4125 0 0
M  V30 11 C 6.5376 -4.4166 0 0
M  V30 12 C 5.7166 -5.2376 0 0
M  V30 13 C 6.4292 -3.175 0 0
M  V30 14 O 7.1375 -3.5875 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 1
M  V30 5 1 1 5
M  V30 6 1 5 6
M  V30 7 1 6 7
M  V30 8 1 7 1
M  V30 9 1 6 13
M  V30 10 1 8 9
M  V30 11 1 9 10
M  V30 12 1 9 11
M  V30 13 1 9 12
M  V30 14 2 13 14
M  V30 15 1 8 13
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 0 ATOMS=(7 8 9 10 11 12 13 14) XBONDS=(1 9) BRKXYZ=(9 6.24 -2.9 0 -
M  V30 6.24 -2.9 0 0 0 0) CSTATE=(4 9 0 0.82 0) LABEL=Boc SAP=(3 13 6 1)
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;

  REQUIRE(m1);

  SECTION("brackets") {
    REQUIRE(m1);
    auto &sgs = getSubstanceGroups(*m1);
    REQUIRE(sgs.size() == 1);
    CHECK(sgs[0].getBrackets().size() == 1);
    sgs[0].clearBrackets();
    CHECK(sgs[0].getBrackets().size() == 0);
  }
  SECTION("cstates") {
    auto &sgs = getSubstanceGroups(*m1);
    REQUIRE(sgs.size() == 1);
    CHECK(sgs[0].getCStates().size() == 1);
    sgs[0].clearCStates();
    CHECK(sgs[0].getCStates().size() == 0);
  }
  SECTION("attachpts") {
    auto &sgs = getSubstanceGroups(*m1);
    REQUIRE(sgs.size() == 1);
    CHECK(sgs[0].getAttachPoints().size() == 1);
    sgs[0].clearAttachPoints();
    CHECK(sgs[0].getAttachPoints().size() == 0);
  }
}

TEST_CASE(
    "GitHub Issue #4434: v2000 SGroups do not generate an \"index\" property",
    "[Sgroups][bug]") {
  SECTION("basics") {
    auto v2000_mol = R"CTAB(
  MJ211300                      

  2  1  0  0  0  0  0  0  0  0999 V2000
    1.4295    0.1449    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4266   -0.6801    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
M  STY  2   1 DAT   2 DAT
M  SAL   1  1   1
M  SDT   1 test sgroup 1                                         
M  SDD   1     0.0000    0.0000    DRU   ALL  0       0  
M  SED   1 sgroupval1
M  SAL   2  1   2
M  SDT   2 test sgroup 2                                         
M  SDD   2     0.0000    0.0000    DRU   ALL  0       0  
M  SED   2 sgroupval2
M  END
)CTAB"_ctab;
    auto v3000_mol = R"CTAB(
  Mrv2113 08202115042D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 2 0 0
M  V30 BEGIN ATOM
M  V30 1 C 2.6684 0.2705 0 0
M  V30 2 C 2.663 -1.2695 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 1) FIELDNAME="test sgroup 1" -
M  V30 FIELDDISP="    0.0000    0.0000    DRU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 FIELDDATA=sgroupval1
M  V30 2 DAT 0 ATOMS=(1 2) FIELDNAME="test sgroup 2" -
M  V30 FIELDDISP="    0.0000    0.0000    DRU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 FIELDDATA=sgroupval2
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;

    REQUIRE(v2000_mol);
    REQUIRE(v3000_mol);

    unsigned int index;
    unsigned int count = 0;

    for (const auto &sg : getSubstanceGroups(*v2000_mol)) {
      ++count;
      TEST_ASSERT(sg.getPropIfPresent("index", index));
      TEST_ASSERT(index == count);
    }

    count = 0;
    for (const auto &sg : getSubstanceGroups(*v3000_mol)) {
      ++count;
      TEST_ASSERT(sg.getPropIfPresent("index", index));
      TEST_ASSERT(index == count);
    }
  }
}