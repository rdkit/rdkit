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

class ScsrMolTest {
 public:
 public:
  bool generateExpectedFiles = false;

  enum class ExpectedStatus {
    Success,
    FailsScsrParsing,
    FailsMolConversion
  };

  class ScsrTest {
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
    ScsrTest(std::string fileNameInit, ExpectedStatus expectedStatusInit,
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

  class ScsrMakeTest {
   public:
    std::string fileName;
    std::string templateFileName;
    unsigned int atomCount;
    unsigned int bondCount;
    unsigned int templateCount;
    bool expectedStatus;

    ScsrMakeTest(std::string fileNameInit, std::string templateFileNameInit,
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
  std::list<ScsrMolTest::ScsrTest> ScsrTests;

  ScsrMolTest() {
    ScsrTests = {
        ScsrMolTest::ScsrTest("ValenceErrorScsr.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 38, 39, 6, 35, 36, 3),
        ScsrMolTest::ScsrTest("ValenceErrorScsr2.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 28, 28, 6, 25, 25, 3),
        ScsrMolTest::ScsrTest("RiboseFullname.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 45, 49, 8, 43, 47, 6),
        ScsrMolTest::ScsrTest("testSCSR.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 64, 66, 15, 57,
                              59, 8),

        ScsrMolTest::ScsrTest("Conjugate.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 91, 91, 14, 87, 87,
                              10),
        ScsrMolTest::ScsrTest("ModifiedPeptide2.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 438, 444, 81, 407,
                              413, 50),
        ScsrMolTest::ScsrTest("DnaBadPairs_NoCh.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 84, 94, 14, 80, 90,
                              10),
        ScsrMolTest::ScsrTest("DnaBadPairs.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 84, 94, 14, 80,
                              90, 10),
        ScsrMolTest::ScsrTest("DnaTest.mol", ExpectedStatus::FailsMolConversion,
                              SCSRBaseHbondOptions::UseSapAll, 254, 300, 14,
                              250, 296, 10),
        ScsrMolTest::ScsrTest("wobblePairs2.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 169, 196, 26,
                              165, 192, 22),
        ScsrMolTest::ScsrTest("wobblePairs.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 169, 196, 26,
                              165, 192, 22),
        ScsrMolTest::ScsrTest("KellanRNA.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 299, 353, 44,
                              295, 349, 40),
        ScsrMolTest::ScsrTest("DnaTest2.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 83, 97, 14, 79,
                              93, 10),
        ScsrMolTest::ScsrTest("DnaTest3.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 165, 194, 26,
                              161, 190, 22),
        ScsrMolTest::ScsrTest("KellanError.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 244, 263, 39,
                              236, 255, 31),

        ScsrMolTest::ScsrTest("TestRNA2_fixed.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 106, 117, 17,
                              104, 115, 15),

        ScsrMolTest::ScsrTest("DnaTest.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 254, 300, 38, 250,
                              296, 34),
        ScsrMolTest::ScsrTest("wobblePairs2.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 169, 196, 26, 165,
                              192, 22),
        ScsrMolTest::ScsrTest("wobblePairs.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 169, 196, 26, 165,
                              192, 22),
        ScsrMolTest::ScsrTest("DnaBadPairs.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 84, 94, 14, 80, 90,
                              10),
        ScsrMolTest::ScsrTest("DnaTest2.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 83, 97, 14, 79, 93,
                              10),
        ScsrMolTest::ScsrTest("DnaTest3.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 165, 194, 26, 161,
                              190, 22),
        ScsrMolTest::ScsrTest("KellanError.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 244, 263, 39, 236,
                              255, 31),
        ScsrMolTest::ScsrTest("TestRNA2_fixed.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::Auto, 106, 117, 17, 104,
                              115, 15),

        ScsrMolTest::ScsrTest(
            "TrastuzumabMaxPlus3Register.mol", ExpectedStatus::Success,
            SCSRBaseHbondOptions::UseSapAll, 7606, 7793, 1451, 7080, 7267, 925),
        ScsrMolTest::ScsrTest(
            "TrastuzumabMaxRegister.mol", ExpectedStatus::Success,
            SCSRBaseHbondOptions::UseSapAll, 7576, 7763, 1445, 7053, 7240, 922),
        ScsrMolTest::ScsrTest("Mixed.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 51, 54, 8, 51,
                              54, 8),
        ScsrMolTest::ScsrTest("CrossLink.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 47, 48, 10, 45,
                              46, 8),
        ScsrMolTest::ScsrTest("cyclic.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 45, 47, 8, 45,
                              47, 8),

        ScsrMolTest::ScsrTest("Triplet.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 30, 30, 6, 27,
                              27, 3),
        ScsrMolTest::ScsrTest("FromBioviaDoc.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 27, 26, 7, 25,
                              24, 5),
        ScsrMolTest::ScsrTest("testSCSR.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 64, 66, 15, 57,
                              59, 8),
        ScsrMolTest::ScsrTest(
            "badAtomName.mol", ExpectedStatus::FailsScsrParsing,
            SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0, 0, 0),
        ScsrMolTest::ScsrTest("badClass.mol", ExpectedStatus::FailsScsrParsing,
                              SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0, 0,
                              0),
        ScsrMolTest::ScsrTest(
            "badClassTemplate.mol", ExpectedStatus::FailsScsrParsing,
            SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0, 0, 0),
        ScsrMolTest::ScsrTest(
            "badMissingTemplate.mol", ExpectedStatus::FailsScsrParsing,
            SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0, 0, 0),
        ScsrMolTest::ScsrTest("obj3dTest.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 27, 26, 7, 25,
                              24, 5),
        ScsrMolTest::ScsrTest("obj3dTest2.mol", ExpectedStatus::Success,
                              SCSRBaseHbondOptions::UseSapAll, 27, 26, 7, 25,
                              24, 5),
        ScsrMolTest::ScsrTest(
            "obj3dFoundTwice.mol", ExpectedStatus::FailsScsrParsing,
            SCSRBaseHbondOptions::UseSapAll, 27, 26, 0, 27, 26, 0),
        ScsrMolTest::ScsrTest(
            "SgroupFoundTwice.mol", ExpectedStatus::FailsScsrParsing,
            SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0, 0, 0),
    };
  }

  void testScsrFiles(const ScsrTest *ScsrTest) {
    BOOST_LOG(rdInfoLog) << "testing scsr  files" << std::endl;

    INFO(ScsrTest->fileName);

    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase +
                        "/Code/GraphMol/FileParsers/test_data/macromols/" +
                        ScsrTest->fileName;
    std::string fOutName = rdbase +
                           "/Code/GraphMol/FileParsers/test_data/macromols/" +
                           ScsrTest->fileName + ".out.mol";
    std::string fOutName2 = rdbase +
                            "/Code/GraphMol/FileParsers/test_data/macromols/" +
                            ScsrTest->fileName + ".scsrout.mol";

    RDKit::v2::FileParsers::MolFileParserParams pp;
    pp.sanitize = true;
    pp.removeHs = false;
    pp.strictParsing = true;

    RDKit::v2::FileParsers::MolFromSCSRParams molFromSCSRParams;
    molFromSCSRParams.includeLeavingGroups = true;
    molFromSCSRParams.scsrBaseHbondOptions = ScsrTest->scsrBaseHbondOptions;

    std::unique_ptr<RDKit::ROMol> mol;
    if (ScsrTest->expectedStatus == ExpectedStatus::Success) {
      REQUIRE_NOTHROW(mol = MolFromSCSRFile(fName, pp, molFromSCSRParams));
    } else {
      REQUIRE_THROWS(mol = MolFromSCSRFile(fName, pp, molFromSCSRParams));
      return;
    }

    RDKit::Chirality::removeNonExplicit3DChirality(*(mol.get()));

    CHECK(mol != nullptr);
    CHECK(mol->getNumAtoms() == ScsrTest->totalAtomCount);
    CHECK(mol->getNumBonds() == ScsrTest->totalBondCount);
    CHECK(getSubstanceGroups(*mol).size() == ScsrTest->sgroupCount);
    MolWriterParams params;
    RDKit::MolToMolFile(*mol, fOutName, params, -1);

    auto scsrMol = SCSRMolFromSCSRFile(fName, pp);

    std::unique_ptr<RDKit::SCSRMol> outScsrMol;

    REQUIRE_NOTHROW(outScsrMol = MolToScsrMol(*(mol.get()), *(scsrMol.get())));

    CHECK(outScsrMol != nullptr);
    CHECK(outScsrMol->getMol()->getNumAtoms() ==
          scsrMol->getMol()->getNumAtoms());
    CHECK(outScsrMol->getMol()->getNumBonds() ==
          scsrMol->getMol()->getNumBonds());
    CHECK(outScsrMol->getTemplateCount() == scsrMol->getTemplateCount());

    SCSRMolToSCSRMolFile(*(outScsrMol.get()), fOutName2);

    std::unique_ptr<RDKit::ROMol> molReadBackIn;
    REQUIRE_NOTHROW(molReadBackIn =
                        MolFromSCSRFile(fOutName2, pp, molFromSCSRParams));

    CHECK(molReadBackIn != nullptr);
    CHECK(molReadBackIn->getNumAtoms() == mol->getNumAtoms());
    CHECK(molReadBackIn->getNumBonds() == mol->getNumBonds());

    // now make the expanded mol in "query" mode - not including any leaving
    // groups

    molFromSCSRParams.includeLeavingGroups = false;
    std::unique_ptr<RDKit::RWMol> molNoLeavingGroups;
    if (ScsrTest->expectedStatus == ExpectedStatus::Success) {
      REQUIRE_NOTHROW(molNoLeavingGroups =
                          MolFromSCSRFile(fName, pp, molFromSCSRParams));
    } else {
      REQUIRE_THROWS(molNoLeavingGroups =
                         MolFromSCSRFile(fName, pp, molFromSCSRParams));
      return;
    }

    CHECK(molNoLeavingGroups != nullptr);
    CHECK(molNoLeavingGroups->getNumAtoms() == ScsrTest->totalQueryAtomCount);
    CHECK(molNoLeavingGroups->getNumBonds() == ScsrTest->totalQueryBondCount);
    CHECK(getSubstanceGroups(*molNoLeavingGroups).size() ==
          ScsrTest->querySgroupCount);

    return;
  }

  void threeLetterCodeTest(const ScsrTest *ScsrTest) {
    BOOST_LOG(rdInfoLog) << "testing scsr  files with three letter codes"
                         << std::endl;

    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase +
                        "/Code/GraphMol/FileParsers/test_data/macromols/" +
                        ScsrTest->fileName;

    RDKit::v2::FileParsers::MolFileParserParams pp;
    pp.sanitize = false;
    pp.removeHs = false;
    pp.strictParsing = false;

    RDKit::v2::FileParsers::MolFromSCSRParams molFromSCSRParams;
    molFromSCSRParams.includeLeavingGroups = true;
    molFromSCSRParams.scsrTemplateNames =
        RDKit::v2::FileParsers::SCSRTemplateNames::AsEntered;
    molFromSCSRParams.scsrBaseHbondOptions = ScsrTest->scsrBaseHbondOptions;

    std::unique_ptr<RDKit::RWMol> mol;
    if (ScsrTest->expectedStatus == ExpectedStatus::Success) {
      REQUIRE_NOTHROW(mol = MolFromSCSRFile(fName, pp, molFromSCSRParams));
    } else {
      REQUIRE_THROWS(mol = MolFromSCSRFile(fName, pp, molFromSCSRParams));
      return;
    }

    RDKit::Chirality::removeNonExplicit3DChirality(*(mol.get()));

    CHECK(mol);
    CHECK(mol->getNumAtoms() == ScsrTest->totalAtomCount);
    CHECK(mol->getNumBonds() == ScsrTest->totalBondCount);
    CHECK(getSubstanceGroups(*mol).size() == ScsrTest->sgroupCount);

    molFromSCSRParams.includeLeavingGroups = true;
    molFromSCSRParams.scsrTemplateNames =
        RDKit::v2::FileParsers::SCSRTemplateNames::UseFirstName;
    std::unique_ptr<RWMol> mol2;
    if (ScsrTest->expectedStatus == ExpectedStatus::Success) {
      REQUIRE_NOTHROW(mol2 = MolFromSCSRFile(fName, pp, molFromSCSRParams));
    } else {
      REQUIRE_THROWS(mol2 = MolFromSCSRFile(fName, pp, molFromSCSRParams));
    }

    // const std::unique_ptr<RDKit::RWMol> mol;
    CHECK(mol2);
    CHECK(mol2->getNumAtoms() == ScsrTest->totalQueryAtomCount);
    CHECK(mol2->getNumBonds() == ScsrTest->totalQueryBondCount);
    CHECK(getSubstanceGroups(*mol2).size() == ScsrTest->querySgroupCount);
  };
};

TEST_CASE("scsrTests", "scsrTests") {
  SECTION("basics") {
    ScsrMolTest ScsrMolTest;
    for (auto ScsrTest : ScsrMolTest.ScsrTests) {
      BOOST_LOG(rdInfoLog) << "Test: " << ScsrTest.fileName << std::endl;

      ScsrMolTest.testScsrFiles(&ScsrTest);
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
    std::list<ScsrMolTest::ScsrTest> scsrTests{
        ScsrMolTest::ScsrTest(
            "PepTla.mol", ScsrMolTest::ExpectedStatus::Success,
            SCSRBaseHbondOptions::UseSapAll, 26, 25, 7, 26, 25, 7),

    };
    ScsrMolTest scsrMolTest;

    for (auto scsrTest : scsrTests) {
      BOOST_LOG(rdInfoLog) << "Test: " << scsrTest.fileName << std::endl;

      scsrMolTest.threeLetterCodeTest(&scsrTest);
    }
  }
}
