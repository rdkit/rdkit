//
//  Copyright (C) 2025 Tad Hurst and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <streambuf>

#include "RDGeneral/test.h"
#include <catch2/catch_all.hpp>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Atropisomers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/lexical_cast.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <filesystem>
using namespace RDKit;
using namespace v2::FileParsers;

using namespace RDKit;

class ScsiMolTest {
 public:
 public:
  bool generateExpectedFiles = false;

  ScsiMolTest() {}

  class ScsiTest {
   public:
    std::string fileName;
    unsigned int totalAtomCount;
    unsigned int totalBondCount;
    unsigned int sgroupCount;
    unsigned int totalQueryAtomCount;
    unsigned int totalQueryBondCount;
    unsigned int querySgroupCount;
    bool scsrExpandResult;
    SCSRBaseHbondOptions scsrBaseHbondOptions;

    ScsiTest(std::string fileNameInit, bool scsrExpandResult,
             SCSRBaseHbondOptions scsrBaseHbondOptions,
             unsigned int totalAtomCountInit, unsigned int totalBondCountInit,
             unsigned int sgroupCountInit, unsigned int totalQueryAtomCountInit,
             unsigned int totalQueryBondCountInit,
             unsigned int querySgroupCountInit = 0)
        : fileName(fileNameInit),

          totalAtomCount(totalAtomCountInit),
          totalBondCount(totalBondCountInit),
          sgroupCount(sgroupCountInit),
          totalQueryAtomCount(totalQueryAtomCountInit),
          totalQueryBondCount(totalQueryBondCountInit),
          querySgroupCount(querySgroupCountInit),
          scsrExpandResult(scsrExpandResult),
          scsrBaseHbondOptions(scsrBaseHbondOptions) {};
  };

  void testScsiFiles(const ScsiTest *scsiTest) {
    BOOST_LOG(rdInfoLog) << "testing scsr  files" << std::endl;

    INFO(scsiTest->fileName);

    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase +
                        "/Code/GraphMol/FileParsers/test_data/macromols/" +
                        scsiTest->fileName;

    RDKit::v2::FileParsers::MolFileParserParams pp;
    pp.sanitize = true;
    pp.removeHs = false;
    pp.strictParsing = true;

    RDKit::v2::FileParsers::MolFromSCSRParams molFromSCSRParams;
    molFromSCSRParams.includeLeavingGroups = true;
    molFromSCSRParams.scsrBaseHbondOptions = scsiTest->scsrBaseHbondOptions;

    std::unique_ptr<RDKit::RWMol> mol;
    if (scsiTest->scsrExpandResult) {
      REQUIRE_NOTHROW(mol = MolFromSCSRFile(fName, pp, molFromSCSRParams));
    } else {
      REQUIRE_THROWS(mol = MolFromSCSRFile(fName, pp, molFromSCSRParams));
      return;
    }

    RDKit::Chirality::removeNonExplicit3DChirality(*(mol.get()));

    CHECK(mol != nullptr);
    CHECK(mol->getNumAtoms() == scsiTest->totalAtomCount);
    CHECK(mol->getNumBonds() == scsiTest->totalBondCount);
    CHECK(getSubstanceGroups(*mol).size() == scsiTest->sgroupCount);

    // now make the expanded mol in "query" mode - not including any leaving
    // groups

    molFromSCSRParams.includeLeavingGroups = false;
    std::unique_ptr<RDKit::RWMol> molNoLeavingGroups;
    if (scsiTest->scsrExpandResult) {
      REQUIRE_NOTHROW(molNoLeavingGroups =
                          MolFromSCSRFile(fName, pp, molFromSCSRParams));
    } else {
      REQUIRE_THROWS(molNoLeavingGroups =
                         MolFromSCSRFile(fName, pp, molFromSCSRParams));
      return;
    }

    CHECK(molNoLeavingGroups != nullptr);
    CHECK(molNoLeavingGroups->getNumAtoms() == scsiTest->totalQueryAtomCount);
    CHECK(molNoLeavingGroups->getNumBonds() == scsiTest->totalQueryBondCount);
    CHECK(getSubstanceGroups(*molNoLeavingGroups).size() ==
          scsiTest->querySgroupCount);
  }

  void threeLetterCodeTest(const ScsiTest *scsiTest) {
    BOOST_LOG(rdInfoLog) << "testing scsr  files with three letter codes"
                         << std::endl;

    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase +
                        "/Code/GraphMol/FileParsers/test_data/macromols/" +
                        scsiTest->fileName;

    RDKit::v2::FileParsers::MolFileParserParams pp;
    pp.sanitize = false;
    pp.removeHs = false;
    pp.strictParsing = false;

    RDKit::v2::FileParsers::MolFromSCSRParams molFromSCSRParams;
    molFromSCSRParams.includeLeavingGroups = true;
    molFromSCSRParams.scsrTemplateNames =
        RDKit::v2::FileParsers::SCSRTemplateNames::AsEntered;
    molFromSCSRParams.scsrBaseHbondOptions = scsiTest->scsrBaseHbondOptions;

    std::unique_ptr<RDKit::RWMol> mol;
    if (scsiTest->scsrExpandResult) {
      REQUIRE_NOTHROW(mol = MolFromSCSRFile(fName, pp, molFromSCSRParams));
    } else {
      REQUIRE_THROWS(mol = MolFromSCSRFile(fName, pp, molFromSCSRParams));
      return;
    }

    RDKit::Chirality::removeNonExplicit3DChirality(*(mol.get()));

    CHECK(mol);
    CHECK(mol->getNumAtoms() == scsiTest->totalAtomCount);
    CHECK(mol->getNumBonds() == scsiTest->totalBondCount);
    CHECK(getSubstanceGroups(*mol).size() == scsiTest->sgroupCount);

    molFromSCSRParams.includeLeavingGroups = true;
    molFromSCSRParams.scsrTemplateNames =
        RDKit::v2::FileParsers::SCSRTemplateNames::UseFirstName;
    std::unique_ptr<RWMol> mol2;
    if (scsiTest->scsrExpandResult) {
      REQUIRE_NOTHROW(mol2 = MolFromSCSRFile(fName, pp, molFromSCSRParams));
    } else {
      REQUIRE_THROWS(mol2 = MolFromSCSRFile(fName, pp, molFromSCSRParams));
    }

    // const std::unique_ptr<RDKit::RWMol> mol;
    CHECK(mol2);
    CHECK(mol2->getNumAtoms() == scsiTest->totalQueryAtomCount);
    CHECK(mol2->getNumBonds() == scsiTest->totalQueryBondCount);
    CHECK(getSubstanceGroups(*mol2).size() == scsiTest->querySgroupCount);
  };
};

TEST_CASE("scsiTests", "scsiTests") {
  SECTION("basics") {
    std::list<ScsiMolTest::ScsiTest> scsiTests{
        ScsiMolTest::ScsiTest("ModifiedPeptide2.mol", true,
                              SCSRBaseHbondOptions::Auto, 438, 444, 81, 407,
                              413, 50),
        ScsiMolTest::ScsiTest("DnaBadPairs_NoCh.mol", true,
                              SCSRBaseHbondOptions::Auto, 84, 94, 14, 80, 90,
                              10),
        ScsiMolTest::ScsiTest("DnaBadPairs.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 84, 94, 14, 80,
                              90, 10),
        ScsiMolTest::ScsiTest("DnaTest.mol", false,
                              SCSRBaseHbondOptions::UseSapAll, 254, 300, 14,
                              250, 296, 10),
        ScsiMolTest::ScsiTest("wobblePairs2.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 169, 196, 26,
                              165, 192, 22),
        ScsiMolTest::ScsiTest("wobblePairs.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 169, 196, 26,
                              165, 192, 22),
        ScsiMolTest::ScsiTest("KellanRNA.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 299, 353, 44,
                              295, 349, 40),
        ScsiMolTest::ScsiTest("DnaTest2.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 83, 97, 14, 79,
                              93, 10),
        ScsiMolTest::ScsiTest("DnaTest3.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 165, 194, 26,
                              161, 190, 22),
        ScsiMolTest::ScsiTest("KellanError.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 244, 263, 39,
                              236, 255, 31),

        ScsiMolTest::ScsiTest("TestRNA2_fixed.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 106, 117, 17,
                              104, 115, 15),

        ScsiMolTest::ScsiTest("DnaTest.mol", true, SCSRBaseHbondOptions::Auto,
                              254, 300, 38, 250, 296, 34),
        ScsiMolTest::ScsiTest("wobblePairs2.mol", true,
                              SCSRBaseHbondOptions::Auto, 169, 196, 26, 165,
                              192, 22),
        ScsiMolTest::ScsiTest("wobblePairs.mol", true,
                              SCSRBaseHbondOptions::Auto, 169, 196, 26, 165,
                              192, 22),
        ScsiMolTest::ScsiTest("DnaBadPairs.mol", true,
                              SCSRBaseHbondOptions::Auto, 84, 94, 14, 80, 90,
                              10),
        ScsiMolTest::ScsiTest("DnaTest2.mol", true, SCSRBaseHbondOptions::Auto,
                              83, 97, 14, 79, 93, 10),
        ScsiMolTest::ScsiTest("DnaTest3.mol", true, SCSRBaseHbondOptions::Auto,
                              165, 194, 26, 161, 190, 22),
        ScsiMolTest::ScsiTest("KellanError.mol", true,
                              SCSRBaseHbondOptions::Auto, 244, 263, 39, 236,
                              255, 31),
        ScsiMolTest::ScsiTest("TestRNA2_fixed.mol", true,
                              SCSRBaseHbondOptions::Auto, 106, 117, 17, 104,
                              115, 15),

        ScsiMolTest::ScsiTest("TrastuzumabMaxPlus3Register.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 7606, 7793, 1451,
                              7080, 7267, 925),
        ScsiMolTest::ScsiTest("TrastuzumabMaxRegister.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 7576, 7763, 1445,
                              7053, 7240, 922),
        ScsiMolTest::ScsiTest("Mixed.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 51, 54, 8, 51,
                              54, 8),
        ScsiMolTest::ScsiTest("CrossLink.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 47, 48, 10, 45,
                              46, 8),
        ScsiMolTest::ScsiTest("cyclic.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 45, 47, 8, 45,
                              47, 8),

        ScsiMolTest::ScsiTest("Triplet.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 30, 30, 6, 27,
                              27, 3),
        ScsiMolTest::ScsiTest("FromBioviaDoc.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 27, 26, 7, 25,
                              24, 5),
        ScsiMolTest::ScsiTest("testSCSR.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 64, 66, 15, 57,
                              59, 8),
        ScsiMolTest::ScsiTest("badAtomName.mol", false,
                              SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0, 0,
                              0),
        ScsiMolTest::ScsiTest("badClass.mol", false,
                              SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0, 0,
                              0),
        ScsiMolTest::ScsiTest("badClassTemplate.mol", false,
                              SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0, 0,
                              0),
        ScsiMolTest::ScsiTest("badMissingTemplate.mol", false,
                              SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0, 0,
                              0),
        ScsiMolTest::ScsiTest("obj3dTest.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 27, 26, 7, 25,
                              24, 5),
        ScsiMolTest::ScsiTest("obj3dTest2.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 27, 26, 7, 25,
                              24, 5),
        ScsiMolTest::ScsiTest("obj3dFoundTwice.mol", false,
                              SCSRBaseHbondOptions::UseSapAll, 27, 26, 0, 27,
                              26, 0),
        ScsiMolTest::ScsiTest("SgroupFoundTwice.mol", false,
                              SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0, 0,
                              0),
    };
    ScsiMolTest scsiMolTest;
    for (auto scsiTest : scsiTests) {
      BOOST_LOG(rdInfoLog) << "Test: " << scsiTest.fileName << std::endl;

      scsiMolTest.testScsiFiles(&scsiTest);
    }
  }
}
TEST_CASE("nestedParens", "nestedParens") {
  SECTION("basics") {
    BOOST_LOG(rdInfoLog) << "testing names with parens" << std::endl;
    std::string filename = "Ugly.mol";
    INFO(filename);

    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/macromols/" + filename;

    RDKit::v2::FileParsers::MolFileParserParams pp;
    pp.sanitize = true;
    pp.removeHs = false;
    pp.strictParsing = true;

    RDKit::v2::FileParsers::MolFromSCSRParams molFromSCSRParams;
    molFromSCSRParams.includeLeavingGroups = true;
    molFromSCSRParams.scsrBaseHbondOptions = SCSRBaseHbondOptions::Auto;

    std::unique_ptr<RDKit::RWMol> mol;
    REQUIRE_NOTHROW(mol = MolFromSCSRFile(fName, pp, molFromSCSRParams));

    CHECK(mol != nullptr);
    CHECK(mol->getNumAtoms() == 140);
    CHECK(mol->getNumBonds() == 151);
    CHECK(getSubstanceGroups(*mol).size() == 11);

    // check that the macro atoms parsed have the correct names.  In the
    // output mol these names will not appear in the SGROUP that defined the
    // macroatom. we willcheck just a couple of them

    std::string sgroupName;
    getSubstanceGroups(*mol)[1].getProp("LABEL", sgroupName);
    std::string expected = "((cPr)O(2S-Me)Et)NGly";
    CHECK(sgroupName.substr(0, expected.length()) == expected);

    getSubstanceGroups(*mol)[5].getProp("LABEL", sgroupName);
    expected = "Phe(b-Me2)";
    CHECK(sgroupName.substr(0, expected.length()) == expected);

    getSubstanceGroups(*mol)[7].getProp("LABEL", sgroupName);
    expected = "(PhO(2S-Me)Et)NGly";
    CHECK(sgroupName.substr(0, expected.length()) == expected);
  }
}

TEST_CASE("threeLetterCodeTest", "threeLetterCodeTest") {
  SECTION("basics") {
    std::list<ScsiMolTest::ScsiTest> scsiTests{
        ScsiMolTest::ScsiTest("PepTla.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 26, 25, 7, 26,
                              25, 7),

    };
    ScsiMolTest scsiMolTest;

    for (auto scsiTest : scsiTests) {
      BOOST_LOG(rdInfoLog) << "Test: " << scsiTest.fileName << std::endl;

      scsiMolTest.threeLetterCodeTest(&scsiTest);
    }
  }
}
