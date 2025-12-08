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
#include <GraphMol/FileParsers/FileWriters.h>
#include <GraphMol/Atropisomers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileWriters.h>

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

  enum class ExpectedStatus {
    Success,
    FailsScsrParsing,
    FailsMolConversion
  };

  class ScsiTest {
   public:
    std::string fileName;
    unsigned int totalAtomCount;
    unsigned int totalBondCount;
    unsigned int sgroupCount;
    unsigned int totalQueryAtomCount;
    unsigned int totalQueryBondCount;
    unsigned int querySgroupCount;
    ExpectedStatus expectedStatus;
    SCSRBaseHbondOptions scsrBaseHbondOptions;
    ScsiTest(std::string fileNameInit, ExpectedStatus expectedStatusInit,
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
          expectedStatus(expectedStatusInit),
          scsrBaseHbondOptions(scsrBaseHbondOptions) {};
  };

  class ScsiMakeTest {
   public:
    std::string fileName;
    std::string templateFileName;
    unsigned int atomCount;
    unsigned int bondCount;
    unsigned int templateCount;
    bool expectedStatus;

    ScsiMakeTest(std::string fileNameInit, std::string templateFileNameInit,
                 unsigned int atomCountInit = 0, unsigned int bondCountInit = 0,
                 unsigned int templateCountInit = 0,
                 bool expectedStatusInit = true)
        : fileName(fileNameInit),
          templateFileName(templateFileNameInit),
          atomCount(atomCountInit),
          bondCount(bondCountInit),
          templateCount(templateCountInit),
          expectedStatus(expectedStatusInit) {};
  };

 public:
  std::list<ScsiMolTest::ScsiTest> scsiTests;
  std::list<ScsiMolTest::ScsiMakeTest> scsiMakeTests;

  ScsiMolTest() {
    scsiTests = {
        ScsiMolTest::ScsiTest("ValenceErrorScsr.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 38, 39, 6, 35, 36, 3),
        ScsiMolTest::ScsiTest("ValenceErrorScsr2.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 28, 28, 6, 25, 25, 3),
        ScsiMolTest::ScsiTest("RiboseFullname.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 45, 49, 8, 43, 47, 6),
        ScsiMolTest::ScsiTest("testSCSR.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 64, 66, 15, 57,
                              59, 8),

        ScsiMolTest::ScsiTest("Conjugate.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 91, 91, 14, 87, 87,
                              10),
        ScsiMolTest::ScsiTest("ModifiedPeptide2.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 438, 444, 81, 407,
                              413, 50),
        ScsiMolTest::ScsiTest("DnaBadPairs_NoCh.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 84, 94, 14, 80, 90,
                              10),
        ScsiMolTest::ScsiTest("DnaBadPairs.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 84, 94, 14, 80,
                              90, 10),
        ScsiMolTest::ScsiTest("DnaTest.mol", ExpectedStatus::FailsMolConversion,
                              SCSRBaseHbondOptions::UseSapAll, 254, 300, 14,
                              250, 296, 10),
        ScsiMolTest::ScsiTest("wobblePairs2.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 169, 196, 26,
                              165, 192, 22),
        ScsiMolTest::ScsiTest("wobblePairs.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 169, 196, 26,
                              165, 192, 22),
        ScsiMolTest::ScsiTest("KellanRNA.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 299, 353, 44,
                              295, 349, 40),
        ScsiMolTest::ScsiTest("DnaTest2.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 83, 97, 14, 79,
                              93, 10),
        ScsiMolTest::ScsiTest("DnaTest3.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 165, 194, 26,
                              161, 190, 22),
        ScsiMolTest::ScsiTest("KellanError.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 244, 263, 39,
                              236, 255, 31),

        ScsiMolTest::ScsiTest("TestRNA2_fixed.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 106, 117, 17,
                              104, 115, 15),

        ScsiMolTest::ScsiTest("DnaTest.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 254, 300, 38, 250,
                              296, 34),
        ScsiMolTest::ScsiTest("wobblePairs2.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 169, 196, 26, 165,
                              192, 22),
        ScsiMolTest::ScsiTest("wobblePairs.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 169, 196, 26, 165,
                              192, 22),
        ScsiMolTest::ScsiTest("DnaBadPairs.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 84, 94, 14, 80, 90,
                              10),
        ScsiMolTest::ScsiTest("DnaTest2.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 83, 97, 14, 79, 93,
                              10),
        ScsiMolTest::ScsiTest("DnaTest3.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 165, 194, 26, 161,
                              190, 22),
        ScsiMolTest::ScsiTest("KellanError.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 244, 263, 39, 236,
                              255, 31),
        ScsiMolTest::ScsiTest("TestRNA2_fixed.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 106, 117, 17, 104,
                              115, 15),

        ScsiMolTest::ScsiTest(
            "TrastuzumabMaxPlus3Register.mol", ExpectedStatus::Success,
            SCSRBaseHbondOptions::UseSapAll, 7606, 7793, 1451, 7080, 7267, 925),
        ScsiMolTest::ScsiTest(
            "TrastuzumabMaxRegister.mol", ExpectedStatus::Success,
            SCSRBaseHbondOptions::UseSapAll, 7576, 7763, 1445, 7053, 7240, 922),
        ScsiMolTest::ScsiTest("Mixed.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 51, 54, 8, 51,
                              54, 8),
        ScsiMolTest::ScsiTest("CrossLink.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 47, 48, 10, 45,
                              46, 8),
        ScsiMolTest::ScsiTest("cyclic.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 45, 47, 8, 45,
                              47, 8),

        ScsiMolTest::ScsiTest("Triplet.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 30, 30, 6, 27,
                              27, 3),
        ScsiMolTest::ScsiTest("FromBioviaDoc.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 27, 26, 7, 25,
                              24, 5),
        ScsiMolTest::ScsiTest("testSCSR.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 64, 66, 15, 57,
                              59, 8),
        ScsiMolTest::ScsiTest(
            "badAtomName.mol", ExpectedStatus::FailsScsrParsing,
            SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0, 0, 0),
        ScsiMolTest::ScsiTest("badClass.mol", ExpectedStatus::FailsScsrParsing,
                              SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0, 0,
                              0),
        ScsiMolTest::ScsiTest(
            "badClassTemplate.mol", ExpectedStatus::FailsScsrParsing,
            SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0, 0, 0),
        ScsiMolTest::ScsiTest(
            "badMissingTemplate.mol", ExpectedStatus::FailsScsrParsing,
            SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0, 0, 0),
        ScsiMolTest::ScsiTest("obj3dTest.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 27, 26, 7, 25,
                              24, 5),
        ScsiMolTest::ScsiTest("obj3dTest2.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 27, 26, 7, 25,
                              24, 5),
        ScsiMolTest::ScsiTest(
            "obj3dFoundTwice.mol", ExpectedStatus::FailsScsrParsing,
            SCSRBaseHbondOptions::UseSapAll, 27, 26, 0, 27, 26, 0),
        ScsiMolTest::ScsiTest(
            "SgroupFoundTwice.mol", ExpectedStatus::FailsScsrParsing,
            SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0, 0, 0),
    };
    scsiMakeTests = {
        ScsiMolTest::ScsiMakeTest("DnaTest2Full.mol",
                                  "DnaTest2FullTemplates.mol", 10, 10, 6),
        ScsiMolTest::ScsiMakeTest("testSCSRFull.mol",
                                  "testSCSRFullTemplates.mol", 8, 7, 6),
    };
  }

  void testScsiFiles(const ScsiTest *scsiTest) {
    BOOST_LOG(rdInfoLog) << "testing scsr  files" << std::endl;

    INFO(scsiTest->fileName);

    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase +
                        "/Code/GraphMol/FileParsers/test_data/macromols/" +
                        scsiTest->fileName;
    std::string fOutName = rdbase +
                           "/Code/GraphMol/FileParsers/test_data/macromols/" +
                           scsiTest->fileName + ".out.mol";

    RDKit::v2::FileParsers::MolFileParserParams pp;
    pp.sanitize = true;
    pp.removeHs = false;
    pp.strictParsing = true;

    RDKit::v2::FileParsers::MolFromSCSRParams molFromSCSRParams;
    molFromSCSRParams.includeLeavingGroups = true;
    molFromSCSRParams.scsrBaseHbondOptions = scsiTest->scsrBaseHbondOptions;

    std::unique_ptr<RDKit::ROMol> mol;
    if (scsiTest->expectedStatus == ExpectedStatus::Success) {
      REQUIRE_NOTHROW(mol = MolFromSCSRFile(fName, pp, molFromSCSRParams));
    } else {
      REQUIRE_THROWS(mol = MolFromSCSRFile(fName, pp, molFromSCSRParams));
      return;
    }

    // auto scsrMol = SCSRMolFromSCSRFile(fName, pp);
    // auto outScsrMol = MolToScsrMol(*(mol.get()), *(scsrMol.get()));

    RDKit::Chirality::removeNonExplicit3DChirality(*(mol.get()));

    CHECK(mol != nullptr);
    CHECK(mol->getNumAtoms() == scsiTest->totalAtomCount);
    CHECK(mol->getNumBonds() == scsiTest->totalBondCount);
    CHECK(getSubstanceGroups(*mol).size() == scsiTest->sgroupCount);
    MolWriterParams params;
    RDKit::MolToMolFile(*mol, fOutName, params, -1);

    // now make the expanded mol in "query" mode - not including any leaving
    // groups

    molFromSCSRParams.includeLeavingGroups = false;
    std::unique_ptr<RDKit::RWMol> molNoLeavingGroups;
    if (scsiTest->expectedStatus == ExpectedStatus::Success) {
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

    return;
  }

  void testMakeScsiFiles(const ScsiMakeTest *scsiMakeTest) {
    BOOST_LOG(rdInfoLog) << "testing scsr  files" << std::endl;

    INFO(scsiMakeTest->fileName);

    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase +
                        "/Code/GraphMol/FileParsers/test_data/macromols/" +
                        scsiMakeTest->fileName;
    std::string tfName = rdbase +
                         "/Code/GraphMol/FileParsers/test_data/macromols/" +
                         scsiMakeTest->templateFileName;
    std::string fOutName = rdbase +
                           "/Code/GraphMol/FileParsers/test_data/macromols/" +
                           scsiMakeTest->fileName + ".createout.mol";

    RDKit::v2::FileParsers::MolFileParserParams pp;
    pp.sanitize = true;
    pp.removeHs = false;
    pp.strictParsing = true;

    std::unique_ptr<RDKit::ROMol> mol;
    mol = MolFromMolFile(fName, pp);
    std::unique_ptr<RDKit::SCSRMol> templateMol;
    templateMol = SCSRMolFromSCSRFile(tfName, pp);

    std::unique_ptr<RDKit::SCSRMol> outScsrMol;
    if (scsiMakeTest->expectedStatus) {
      REQUIRE_NOTHROW(outScsrMol =
                          MolToScsrMol(*(mol.get()), *(templateMol.get())));
    } else {
      REQUIRE_THROWS(outScsrMol =
                         MolToScsrMol(*(mol.get()), *(templateMol.get())));
      return;
    }

    CHECK(outScsrMol != nullptr);
    CHECK(outScsrMol->getMol()->getNumAtoms() == scsiMakeTest->atomCount);
    CHECK(outScsrMol->getMol()->getNumBonds() == scsiMakeTest->bondCount);
    CHECK(outScsrMol->getTemplateCount() == scsiMakeTest->templateCount);

    SCSRMolToSCSRMolFile(*(outScsrMol.get()), fOutName);

    return;
  }

  void testScsiMols(const ScsiTest *scsiTest) {
    BOOST_LOG(rdInfoLog)
        << "testing scsr mol generation and generation of scsr files from scsrMols"
        << std::endl;

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

    std::unique_ptr<RDKit::SCSRMol> scsrMol;
    std::string outScsrName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/macromols/" +
        scsiTest->fileName + ".scsrOut.mol";
    if (scsiTest->expectedStatus == ExpectedStatus::FailsScsrParsing) {
      REQUIRE_THROWS(scsrMol = SCSRMolFromSCSRFile(fName, pp));
      return;
    } else {
      REQUIRE_NOTHROW(scsrMol = SCSRMolFromSCSRFile(fName, pp));
    }

    SCSRMolToSCSRMolFile(*(scsrMol.get()), outScsrName);

    if (scsiTest->expectedStatus == ExpectedStatus::Success) {
      REQUIRE_NOTHROW(mol =
                          MolFromSCSRFile(outScsrName, pp, molFromSCSRParams));
    } else {
      REQUIRE_THROWS(mol = MolFromSCSRFile(outScsrName, pp, molFromSCSRParams));
      return;
    }

    RDKit::Chirality::removeNonExplicit3DChirality(*(mol.get()));

    CHECK(mol != nullptr);
    CHECK(mol->getNumAtoms() == scsiTest->totalAtomCount);
    CHECK(mol->getNumBonds() == scsiTest->totalBondCount);
    CHECK(getSubstanceGroups(*mol).size() == scsiTest->sgroupCount);

    return;
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
    if (scsiTest->expectedStatus == ExpectedStatus::Success) {
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
    if (scsiTest->expectedStatus == ExpectedStatus::Success) {
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
    ScsiMolTest scsiMolTest;
    for (auto scsiTest : scsiMolTest.scsiTests) {
      BOOST_LOG(rdInfoLog) << "Test: " << scsiTest.fileName << std::endl;

      scsiMolTest.testScsiFiles(&scsiTest);
    }
  }
}
TEST_CASE("scsiMols", "scsiMols") {
  SECTION("basics") {
    ScsiMolTest scsiMolTest;
    for (auto scsiTest : scsiMolTest.scsiTests) {
      BOOST_LOG(rdInfoLog) << "Test: " << scsiTest.fileName << std::endl;

      scsiMolTest.testScsiMols(&scsiTest);
    }
  }
}

TEST_CASE("makeScsrMols", "makeScsrMols") {
  SECTION("basics") {
    ScsiMolTest scsiMolTest;
    for (auto scsiMakeTest : scsiMolTest.scsiMakeTests) {
      BOOST_LOG(rdInfoLog) << "Test: " << scsiMakeTest.fileName << std::endl;

      scsiMolTest.testMakeScsiFiles(&scsiMakeTest);
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
    std::string expected = "AA_2_((cPr)O(2S-Me)Et)NGly+((cPr)O(2S-Me)Et)NGly";
    CHECK(sgroupName.substr(0, expected.length()) == expected);

    getSubstanceGroups(*mol)[5].getProp("LABEL", sgroupName);
    expected = "AA_6_Phe(b-Me2)+Phe(b-Me2)";
    CHECK(sgroupName.substr(0, expected.length()) == expected);

    getSubstanceGroups(*mol)[7].getProp("LABEL", sgroupName);
    expected = "AA_8_(PhO(2S-Me)Et)NGly+(PhO(2S-Me)Et)NGly";
    CHECK(sgroupName.substr(0, expected.length()) == expected);
  }
}

TEST_CASE("threeLetterCodeTest", "threeLetterCodeTest") {
  SECTION("basics") {
    std::list<ScsiMolTest::ScsiTest> scsiTests{
        ScsiMolTest::ScsiTest(
            "PepTla.mol", ScsiMolTest::ExpectedStatus::Success,
            SCSRBaseHbondOptions::UseSapAll, 26, 25, 7, 26, 25, 7),

    };
    ScsiMolTest scsiMolTest;

    for (auto scsiTest : scsiTests) {
      BOOST_LOG(rdInfoLog) << "Test: " << scsiTest.fileName << std::endl;

      scsiMolTest.threeLetterCodeTest(&scsiTest);
    }
  }
}
