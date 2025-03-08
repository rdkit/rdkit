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
    std::string fileName;
    unsigned int totalAtomCount;
    unsigned int totalBondCount;
    unsigned int totalQueryAtomCount;
    unsigned int totalQueryBondCount;
    bool scsrExpandResult;
    SCSRBaseHbondOptions scsrBaseHbondOptions;

    ScsiTest(std::string fileNameInit, bool scsrExpandResult,
             SCSRBaseHbondOptions scsrBaseHbondOptions,
             unsigned int totalAtomCountInit, unsigned int totalBondCountInit,
             unsigned int totalQueryAtomCountInit,
             unsigned int totalQueryBondCountInit)
        : fileName(fileNameInit),

          totalAtomCount(totalAtomCountInit),
          totalBondCount(totalBondCountInit),
          totalQueryAtomCount(totalQueryAtomCountInit),
          totalQueryBondCount(totalQueryBondCountInit),
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

    INFO(scsiTest->fileName);

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

    molFromSCSRParams.includeLeavingGroups = false;
    std::unique_ptr<RDKit::RWMol> molFromStream;
    if (scsiTest->scsrExpandResult) {
      REQUIRE_NOTHROW(molFromStream =
                          MolFromSCSRFile(fName, pp, molFromSCSRParams));
    } else {
      REQUIRE_THROWS(molFromStream =
                         MolFromSCSRFile(fName, pp, molFromSCSRParams));
      return;
    }

    {
      std::string expectedMolName = fName + ".expected2.sdf";
      std::string outMolStr = "";
      try {
        outMolStr = MolToMolBlock(*molFromStream, true, 0, true, true);
      } catch (const RDKit::KekulizeException &) {
        outMolStr = "";
      }
      if (outMolStr == "") {
        outMolStr = MolToMolBlock(*molFromStream, true, 0, false,
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
        ScsiMolTest::ScsiTest("DnaBadPairs_NoCh.mol", true,
                              SCSRBaseHbondOptions::Auto, 84, 94, 80, 90),
        ScsiMolTest::ScsiTest("DnaBadPairs.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 84, 94, 80, 90),
        ScsiMolTest::ScsiTest("DnaTest.mol", false,
                              SCSRBaseHbondOptions::UseSapAll, 254, 300, 250,
                              296),
        ScsiMolTest::ScsiTest("wobblePairs2.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 169, 196, 165,
                              192),
        ScsiMolTest::ScsiTest("wobblePairs.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 169, 196, 165,
                              192),
        ScsiMolTest::ScsiTest("KellanRNA.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 299, 353, 295,
                              349),
        ScsiMolTest::ScsiTest("DnaTest2.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 83, 97, 79, 93),
        ScsiMolTest::ScsiTest("DnaTest3.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 165, 194, 161,
                              190),
        ScsiMolTest::ScsiTest("KellanError.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 244, 263, 236,
                              255),
        ScsiMolTest::ScsiTest("TestRNA2_fixed.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 106, 117, 104,
                              115),

        ScsiMolTest::ScsiTest("DnaTest.mol", true, SCSRBaseHbondOptions::Auto,
                              254, 300, 250, 296),
        ScsiMolTest::ScsiTest("wobblePairs2.mol", true,
                              SCSRBaseHbondOptions::Auto, 169, 196, 165, 192),
        ScsiMolTest::ScsiTest("wobblePairs.mol", true,
                              SCSRBaseHbondOptions::Auto, 169, 196, 165, 192),
        ScsiMolTest::ScsiTest("KellanRNA.mol", true, SCSRBaseHbondOptions::Auto,
                              299, 353, 295, 349),
        ScsiMolTest::ScsiTest("DnaBadPairs.mol", true,
                              SCSRBaseHbondOptions::Auto, 84, 94, 80, 90),
        ScsiMolTest::ScsiTest("DnaTest2.mol", true, SCSRBaseHbondOptions::Auto,
                              83, 97, 79, 93),
        ScsiMolTest::ScsiTest("DnaTest3.mol", true, SCSRBaseHbondOptions::Auto,
                              165, 194, 161, 190),
        ScsiMolTest::ScsiTest("KellanError.mol", true,
                              SCSRBaseHbondOptions::Auto, 244, 263, 236, 255),
        ScsiMolTest::ScsiTest("TestRNA2_fixed.mol", true,
                              SCSRBaseHbondOptions::Auto, 106, 117, 104, 115),

        ScsiMolTest::ScsiTest("TrastuzumabMaxPlus3Register.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 7606, 7793, 7080,
                              7267),
        ScsiMolTest::ScsiTest("TrastuzumabMaxRegister.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 7576, 7763, 7053,
                              7240),
        ScsiMolTest::ScsiTest("Mixed.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 51, 54, 51, 54),
        ScsiMolTest::ScsiTest("CrossLink.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 47, 48, 45, 46),
        ScsiMolTest::ScsiTest("cyclic.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 45, 47, 45, 47),

        ScsiMolTest::ScsiTest("Triplet.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 30, 30, 27, 27),
        ScsiMolTest::ScsiTest("FromBioviaDoc.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 27, 26, 25, 24),
        ScsiMolTest::ScsiTest("testSCSR.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 64, 66, 57, 59),
        ScsiMolTest::ScsiTest("badAtomName.mol", false,
                              SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0),
        ScsiMolTest::ScsiTest("badClass.mol", false,
                              SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0),
        ScsiMolTest::ScsiTest("badClassTemplate.mol", false,
                              SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0),
        ScsiMolTest::ScsiTest("badMissingTemplate.mol", false,
                              SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0),
        ScsiMolTest::ScsiTest("obj3dTest.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 27, 26, 25, 24),
        ScsiMolTest::ScsiTest("obj3dTest2.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 27, 26, 25, 24),
        ScsiMolTest::ScsiTest("obj3dFoundTwice.mol", false,
                              SCSRBaseHbondOptions::UseSapAll, 27, 26, 27, 26),
        ScsiMolTest::ScsiTest("SgroupFoundTwice.mol", false,
                              SCSRBaseHbondOptions::UseSapAll, 0, 0, 0, 0),
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
        ScsiMolTest::ScsiTest("PepTla.mol", true,
                              SCSRBaseHbondOptions::UseSapAll, 26, 25, 26, 25),

    };
    ScsiMolTest scsiMolTest;

    for (auto scsiTest : scsiTests) {
      BOOST_LOG(rdInfoLog) << "Test: " << scsiTest.fileName << std::endl;

      scsiMolTest.threeLetterCodeTest(&scsiTest);
    }
  }
}
