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
  std::string testToRun;
  bool generateExpectedFiles;

  ScsiMolTest() {
    testToRun = "";
    generateExpectedFiles = false;
  }

  class ScsiTest {
   public:
    unsigned int atomCount;
    unsigned int bondCount;
    std::string fileName;
    unsigned int templateCount;
    unsigned int totalAtomCount;
    unsigned int totalBondCount;
    unsigned int totalQueryAtomCount;
    unsigned int totalQueryBondCount;
    bool scsrParseResult;
    bool scsrExpandResult;
    ScsrBaseHbondOptions scsrBaseHbondOptions;

    ScsiTest(std::string fileNameInit, bool scsrParseResult,
             bool scsrExpandResult, ScsrBaseHbondOptions scsrBaseHbondOptions,
             unsigned int atomCountInit, unsigned int bondCountInit,
             unsigned int templateCountInit, unsigned int totalAtomCountInit,
             unsigned int totalBondCountInit,
             unsigned int totalQueryAtomCountInit,
             unsigned int totalQueryBondCountInit)
        : atomCount(atomCountInit),
          bondCount(bondCountInit),
          fileName(fileNameInit),
          templateCount(templateCountInit),
          totalAtomCount(totalAtomCountInit),
          totalBondCount(totalBondCountInit),
          totalQueryAtomCount(totalQueryAtomCountInit),
          totalQueryBondCount(totalQueryBondCountInit),
          scsrParseResult(scsrParseResult),
          scsrExpandResult(scsrExpandResult),
          scsrBaseHbondOptions(scsrBaseHbondOptions) {};
  };

  void generateNewExpectedFilesIfSoSpecified(std::string filename,
                                             std::string dataToWrite) {
    if (generateExpectedFiles) {
      std::ofstream out;
      out.open(filename);
      out << dataToWrite;
    }
  }

  std::string GetExpectedValue(std::string expectedFileName) {
    std::stringstream expectedMolStr;
    std::ifstream in;
    in.open(expectedFileName);
    expectedMolStr << in.rdbuf();
    return expectedMolStr.str();
  }

  void testScsiFiles(const ScsiTest *scsiTest) {
    BOOST_LOG(rdInfoLog) << "testing scsr  files" << std::endl;

    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase +
                        "/Code/GraphMol/FileParsers/test_data/macromols/" +
                        scsiTest->fileName;

    RDKit::v2::FileParsers::MolFileParserParams pp;
    pp.sanitize = false;
    pp.removeHs = false;
    pp.strictParsing = false;
    std::unique_ptr<RDKit::SCSRMol> scsrMol, scsrMolFromStream;
    if (scsiTest->scsrParseResult) {
      REQUIRE_NOTHROW(scsrMol = SCSRMolFromScsrFile(fName, pp));
    } else {
      REQUIRE_THROWS(scsrMol = SCSRMolFromScsrFile(fName, pp));
      return;
    }

    RDKit::Chirality::removeNonExplicit3DChirality(*scsrMol->getMol());

    CHECK(scsrMol);
    CHECK(scsrMol->getMol()->getNumAtoms() == scsiTest->atomCount);
    CHECK(scsrMol->getMol()->getNumBonds() == scsiTest->bondCount);
    CHECK(scsrMol->getTemplateCount() == scsiTest->templateCount);

    RDKit::v2::FileParsers::MolFromScsrParams molFromScsrParams;
    molFromScsrParams.includeLeavingGroups = true;
    molFromScsrParams.scsrBaseHbondOptions = scsiTest->scsrBaseHbondOptions;

    std::unique_ptr<RWMol> mol;
    if (scsiTest->scsrExpandResult) {
      REQUIRE_NOTHROW(mol = MolFromSCSRMol(scsrMol.get(), molFromScsrParams));
    } else {
      REQUIRE_THROWS(mol = MolFromSCSRMol(scsrMol.get(), molFromScsrParams));
      return;
    }

    CHECK(mol != nullptr);
    CHECK(mol->getNumAtoms() == scsiTest->totalAtomCount);
    CHECK(mol->getNumBonds() == scsiTest->totalBondCount);

    {
      std::string expectedMolName = fName + ".expected.sdf";
      std::string outMolStr = "";
      try {
        outMolStr = MolToMolBlock(*mol, true, 0, true, true);
      } catch (const RDKit::KekulizeException &) {
        outMolStr = "";
      }
      if (outMolStr == "") {
        outMolStr = MolToMolBlock(*mol, true, 0, false,
                                  true);  // try without kekule'ing
      }

      generateNewExpectedFilesIfSoSpecified(fName + ".NEW.sdf", outMolStr);

      CHECK(GetExpectedValue(expectedMolName) == outMolStr);
    }
    // now make the expanded mol in "query" mode - not including any leaving
    // groups

    std::unique_ptr<RWMol> mol2;
    molFromScsrParams.includeLeavingGroups = false;
    if (scsiTest->scsrExpandResult) {
      REQUIRE_NOTHROW(mol2 = MolFromSCSRMol(scsrMol.get(), molFromScsrParams));
    } else {
      REQUIRE_THROWS(mol2 = MolFromSCSRMol(scsrMol.get(), molFromScsrParams));
      return;
    }

    // const std::unique_ptr<RDKit::RWMol> mol;
    CHECK(mol2);
    CHECK(mol2->getNumAtoms() == scsiTest->totalQueryAtomCount);
    CHECK(mol2->getNumBonds() == scsiTest->totalQueryBondCount);

    {
      std::string expectedMolName = fName + ".expected2.sdf";
      std::string outMolStr = "";
      try {
        outMolStr = MolToMolBlock(*mol2, true, 0, true, true);
      } catch (const RDKit::KekulizeException &) {
        outMolStr = "";
      }
      if (outMolStr == "") {
        outMolStr = MolToMolBlock(*mol2, true, 0, false,
                                  true);  // try without kekule'ing
      }

      generateNewExpectedFilesIfSoSpecified(fName + ".NEW2.sdf", outMolStr);

      CHECK(GetExpectedValue(expectedMolName) == outMolStr);
    }
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

    std::unique_ptr<RDKit::SCSRMol> scsrMol;
    if (scsiTest->scsrParseResult) {
      CHECK_NOTHROW(scsrMol = SCSRMolFromScsrFile(fName, pp));
    } else {
      CHECK_THROWS(scsrMol = SCSRMolFromScsrFile(fName, pp));
    }
    RDKit::Chirality::removeNonExplicit3DChirality(*scsrMol->getMol());

    REQUIRE(scsrMol);
    CHECK(scsrMol->getMol()->getNumAtoms() == scsiTest->atomCount);
    CHECK(scsrMol->getMol()->getNumBonds() == scsiTest->bondCount);
    CHECK(scsrMol->getTemplateCount() == scsiTest->templateCount);

    RDKit::v2::FileParsers::MolFromScsrParams molFromScsrParams;
    molFromScsrParams.includeLeavingGroups = true;
    molFromScsrParams.scsrTemplateNames =
        RDKit::v2::FileParsers::ScsrTemplateNamesAsEntered;
    molFromScsrParams.scsrBaseHbondOptions = scsiTest->scsrBaseHbondOptions;
    std::unique_ptr<RWMol> mol;
    if (scsiTest->scsrExpandResult) {
      REQUIRE_NOTHROW(mol = MolFromSCSRMol(scsrMol.get(), molFromScsrParams));
    } else {
      REQUIRE_THROWS(mol = MolFromSCSRMol(scsrMol.get(), molFromScsrParams));
      return;
    }

    CHECK(mol);
    CHECK(mol->getNumAtoms() == scsiTest->totalAtomCount);
    CHECK(mol->getNumBonds() == scsiTest->totalBondCount);

    {
      std::string expectedMolName = fName + ".expected.sdf";
      std::string outMolStr = "";
      try {
        outMolStr = MolToMolBlock(*mol, true, 0, true, true);
      } catch (const RDKit::KekulizeException &) {
        outMolStr = "";
      }
      if (outMolStr == "") {
        outMolStr = MolToMolBlock(*mol, true, 0, false,
                                  true);  // try without kekule'ing
      }

      generateNewExpectedFilesIfSoSpecified(fName + ".NEW.sdf", outMolStr);

      CHECK(GetExpectedValue(expectedMolName) == outMolStr);
    }
    // now make the expanded mol in "query" mode - not including any leaving
    // groups

    molFromScsrParams.includeLeavingGroups = true;
    molFromScsrParams.scsrTemplateNames =
        RDKit::v2::FileParsers::ScsrTemplateNamesUseFirstName;
    std::unique_ptr<RWMol> mol2;
    if (scsiTest->scsrExpandResult) {
      REQUIRE_NOTHROW(mol2 = MolFromSCSRMol(scsrMol.get(), molFromScsrParams));
    } else {
      REQUIRE_THROWS(mol2 = MolFromSCSRMol(scsrMol.get(), molFromScsrParams));
    }

    // const std::unique_ptr<RDKit::RWMol> mol;
    CHECK(mol2);
    CHECK(mol2->getNumAtoms() == scsiTest->totalQueryAtomCount);
    CHECK(mol2->getNumBonds() == scsiTest->totalQueryBondCount);

    {
      std::string expectedMolName = fName + ".expected2.sdf";
      std::string outMolStr = "";
      try {
        outMolStr = MolToMolBlock(*mol2, true, 0, true, true);
      } catch (const RDKit::KekulizeException &) {
        outMolStr = "";
      }
      if (outMolStr == "") {
        outMolStr = MolToMolBlock(*mol2, true, 0, false,
                                  true);  // try without kekule'ing
      }

      generateNewExpectedFilesIfSoSpecified(fName + ".NEW2.sdf", outMolStr);

      CHECK(GetExpectedValue(expectedMolName) == outMolStr);
    }
  };
};

TEST_CASE("scsiTests", "scsiTests") {
  SECTION("basics") {
    std::list<ScsiMolTest::ScsiTest> scsiTests{
        ScsiMolTest::ScsiTest("DnaBadPairs_NoCh.mol", true, true,
                              ScsrBaseHbondOptionsAuto, 10, 10, 7, 84, 94, 80,
                              90),
        ScsiMolTest::ScsiTest("DnaBadPairs.mol", true, true,
                              ScsrBaseHbondOptionsUseSapAll, 10, 10, 7, 84, 94,
                              80, 90),
        ScsiMolTest::ScsiTest("DnaTest.mol", true, false,
                              ScsrBaseHbondOptionsUseSapAll, 34, 38, 7, 254,
                              300, 250, 296),
        ScsiMolTest::ScsiTest("wobblePairs2.mol", true, true,
                              ScsrBaseHbondOptionsUseSapAll, 22, 24, 7, 169,
                              196, 165, 192),
        ScsiMolTest::ScsiTest("wobblePairs.mol", true, true,
                              ScsrBaseHbondOptionsUseSapAll, 22, 24, 7, 169,
                              196, 165, 192),
        ScsiMolTest::ScsiTest("KellanRNA.mol", true, true,
                              ScsrBaseHbondOptionsUseSapAll, 40, 45, 6, 299,
                              353, 295, 349),
        ScsiMolTest::ScsiTest("DnaTest2.mol", true, true,
                              ScsrBaseHbondOptionsUseSapAll, 10, 10, 6, 83, 97,
                              79, 93),
        ScsiMolTest::ScsiTest("DnaTest3.mol", true, true,
                              ScsrBaseHbondOptionsUseSapAll, 22, 24, 8, 165,
                              194, 161, 190),
        ScsiMolTest::ScsiTest("KellanError.mol", true, true,
                              ScsrBaseHbondOptionsUseSapAll, 31, 30, 15, 244,
                              263, 236, 255),
        ScsiMolTest::ScsiTest("TestRNA2_fixed.mol", true, true,
                              ScsrBaseHbondOptionsUseSapAll, 15, 14, 5, 106,
                              117, 104, 115),

        ScsiMolTest::ScsiTest("DnaTest.mol", true, true,
                              ScsrBaseHbondOptionsAuto, 34, 38, 7, 254, 300,
                              250, 296),
        ScsiMolTest::ScsiTest("wobblePairs2.mol", true, true,
                              ScsrBaseHbondOptionsAuto, 22, 24, 7, 169, 196,
                              165, 192),
        ScsiMolTest::ScsiTest("wobblePairs.mol", true, true,
                              ScsrBaseHbondOptionsAuto, 22, 24, 7, 169, 196,
                              165, 192),
        ScsiMolTest::ScsiTest("KellanRNA.mol", true, true,
                              ScsrBaseHbondOptionsAuto, 40, 45, 6, 299, 353,
                              295, 349),
        ScsiMolTest::ScsiTest("DnaBadPairs.mol", true, true,
                              ScsrBaseHbondOptionsAuto, 10, 10, 7, 84, 94, 80,
                              90),
        ScsiMolTest::ScsiTest("DnaTest2.mol", true, true,
                              ScsrBaseHbondOptionsAuto, 10, 10, 6, 83, 97, 79,
                              93),
        ScsiMolTest::ScsiTest("DnaTest3.mol", true, true,
                              ScsrBaseHbondOptionsAuto, 22, 24, 8, 165, 194,
                              161, 190),
        ScsiMolTest::ScsiTest("KellanError.mol", true, true,
                              ScsrBaseHbondOptionsAuto, 31, 30, 15, 244, 263,
                              236, 255),
        ScsiMolTest::ScsiTest("TestRNA2_fixed.mol", true, true,
                              ScsrBaseHbondOptionsAuto, 15, 14, 5, 106, 117,
                              104, 115),

        ScsiMolTest::ScsiTest("TrastuzumabMaxPlus3Register.mol", true, true,
                              ScsrBaseHbondOptionsUseSapAll, 925, 933, 21, 7606,
                              7793, 7080, 7267),
        ScsiMolTest::ScsiTest("TrastuzumabMaxRegister.mol", true, true,
                              ScsrBaseHbondOptionsUseSapAll, 922, 930, 21, 7576,
                              7763, 7053, 7240),
        ScsiMolTest::ScsiTest("Mixed.mol", true, true,
                              ScsrBaseHbondOptionsUseSapAll, 14, 16, 5, 51, 54,
                              51, 54),
        ScsiMolTest::ScsiTest("CrossLink.mol", true, true,
                              ScsrBaseHbondOptionsUseSapAll, 8, 8, 5, 47, 48,
                              45, 46),
        ScsiMolTest::ScsiTest("cyclic.mol", true, true,
                              ScsrBaseHbondOptionsUseSapAll, 8, 9, 5, 45, 47,
                              45, 47),

        ScsiMolTest::ScsiTest("Triplet.mol", true, true,
                              ScsrBaseHbondOptionsUseSapAll, 3, 2, 3, 30, 30,
                              27, 27),
        ScsiMolTest::ScsiTest("FromBioviaDoc.mol", true, true,
                              ScsrBaseHbondOptionsUseSapAll, 5, 4, 3, 27, 26,
                              25, 24),
        ScsiMolTest::ScsiTest("testScsr.mol", true, true,
                              ScsrBaseHbondOptionsUseSapAll, 8, 7, 6, 64, 66,
                              57, 59),
        ScsiMolTest::ScsiTest("badAtomName.mol", false, false,
                              ScsrBaseHbondOptionsUseSapAll, 8, 7, 6, 0, 0, 0,
                              0),
        ScsiMolTest::ScsiTest("badClass.mol", false, false,
                              ScsrBaseHbondOptionsUseSapAll, 8, 7, 6, 0, 0, 0,
                              0),
        ScsiMolTest::ScsiTest("badClassTemplate.mol", false, false,
                              ScsrBaseHbondOptionsUseSapAll, 8, 7, 6, 0, 0, 0,
                              0),
        ScsiMolTest::ScsiTest("badMissingTemplate.mol", false, false,
                              ScsrBaseHbondOptionsUseSapAll, 8, 7, 6, 0, 0, 0,
                              0),
        ScsiMolTest::ScsiTest("obj3dTest.mol", true, true,
                              ScsrBaseHbondOptionsUseSapAll, 5, 4, 3, 27, 26,
                              25, 24),
        ScsiMolTest::ScsiTest("obj3dTest2.mol", true, true,
                              ScsrBaseHbondOptionsUseSapAll, 5, 4, 3, 27, 26,
                              25, 24),
        ScsiMolTest::ScsiTest("obj3dFoundTwice.mol", false, false,
                              ScsrBaseHbondOptionsUseSapAll, 5, 4, 3, 27, 26,
                              27, 26),
        ScsiMolTest::ScsiTest("SgroupFoundTwice.mol", false, false,
                              ScsrBaseHbondOptionsUseSapAll, 5, 4, 3, 0, 0, 0,
                              0),
    };
    ScsiMolTest scsiMolTest;
    for (auto scsiTest : scsiTests) {
      BOOST_LOG(rdInfoLog) << "Test: " << scsiTest.fileName << std::endl;

      scsiMolTest.testScsiFiles(&scsiTest);
    }
  }
}

TEST_CASE("threeLetterCodeTest", "threeLetterCodeTest") {
  SECTION("basics") {
    std::list<ScsiMolTest::ScsiTest> scsiTests{
        ScsiMolTest::ScsiTest("PepTla.mol", true, true,
                              ScsrBaseHbondOptionsUseSapAll, 4, 3, 4, 26, 25,
                              26, 25),

    };
    ScsiMolTest scsiMolTest;

    for (auto scsiTest : scsiTests) {
      BOOST_LOG(rdInfoLog) << "Test: " << scsiTest.fileName << std::endl;

      scsiMolTest.threeLetterCodeTest(&scsiTest);
    }
  }
}
