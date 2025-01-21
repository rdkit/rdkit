//
//  Copyright (C) 2002-2021 Collaboartive Drug Discovery and other RDKit
//  contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Atropisomers.h>

#include <string>
#include <fstream>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/lexical_cast.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <filesystem>
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

 public:
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
    bool expectedResult;

    ScsiTest(std::string fileNameInit, bool expectedResultInit,
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
          expectedResult(expectedResultInit) {};
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

    try {
      RDKit::v2::FileParsers::MolFileParserParams pp;
      pp.sanitize = false;
      pp.removeHs = false;
      pp.strictParsing = false;

      std::unique_ptr<RDKit::SCSRMol> scsrMol = SCSRMolFromScsrFile(fName, pp);
      RDKit::Chirality::removeNonExplicit3DChirality(*scsrMol->getMol());

      TEST_ASSERT(scsrMol != nullptr);
      TEST_ASSERT(scsrMol->getMol()->getNumAtoms() == scsiTest->atomCount)
      TEST_ASSERT(scsrMol->getMol()->getNumBonds() == scsiTest->bondCount)
      TEST_ASSERT(scsrMol->getTemplateCount() == scsiTest->templateCount)

      std::fstream inStream;
      inStream.open(fName);

      unsigned int line = 0;
      auto scsrMolFromStream = SCSRMolFromScsrDataStream(inStream, line, pp);

      inStream.close();

      RDKit::Chirality::removeNonExplicit3DChirality(
          *scsrMolFromStream->getMol());

      TEST_ASSERT(scsrMolFromStream != nullptr);
      TEST_ASSERT(scsrMolFromStream->getMol()->getNumAtoms() ==
                  scsiTest->atomCount)
      TEST_ASSERT(scsrMolFromStream->getMol()->getNumBonds() ==
                  scsiTest->bondCount)
      TEST_ASSERT(scsrMolFromStream->getTemplateCount() ==
                  scsiTest->templateCount)

      RDKit::v2::FileParsers::MolFromScsrParams molFromScsrParams;
      molFromScsrParams.includeLeavingGroups = true;
      auto mol = MolFromSCSRMol(*scsrMol, molFromScsrParams);

      TEST_ASSERT(mol != nullptr)
      TEST_ASSERT(mol->getNumAtoms() == scsiTest->totalAtomCount)
      TEST_ASSERT(mol->getNumBonds() == scsiTest->totalBondCount)

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

        TEST_ASSERT(GetExpectedValue(expectedMolName) == outMolStr);
      }
      // now make the expanded mol in "query" mode - not including any leaving
      // groups

      molFromScsrParams.includeLeavingGroups = false;
      auto mol2 = MolFromSCSRMol(*scsrMol, molFromScsrParams);

      // const std::unique_ptr<RDKit::RWMol> mol;
      TEST_ASSERT(mol2 != nullptr)
      TEST_ASSERT(mol2->getNumAtoms() == scsiTest->totalQueryAtomCount)
      TEST_ASSERT(mol2->getNumBonds() == scsiTest->totalQueryBondCount)

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

        TEST_ASSERT(GetExpectedValue(expectedMolName) == outMolStr);
      }
    } catch (const std::exception &e) {
      if (scsiTest->expectedResult != false) {
        throw;
      }
      return;
    }
    TEST_ASSERT(scsiTest->expectedResult == true);
  }

  void threeLetterCodeTest(const ScsiTest *scsiTest) {
    BOOST_LOG(rdInfoLog) << "testing scsr  files with three letter codes"
                         << std::endl;

    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase +
                        "/Code/GraphMol/FileParsers/test_data/macromols/" +
                        scsiTest->fileName;

    try {
      RDKit::v2::FileParsers::MolFileParserParams pp;
      pp.sanitize = false;
      pp.removeHs = false;
      pp.strictParsing = false;

      std::unique_ptr<RDKit::SCSRMol> scsrMol = SCSRMolFromScsrFile(fName, pp);
      RDKit::Chirality::removeNonExplicit3DChirality(*scsrMol->getMol());

      TEST_ASSERT(scsrMol != nullptr);
      TEST_ASSERT(scsrMol->getMol()->getNumAtoms() == scsiTest->atomCount)
      TEST_ASSERT(scsrMol->getMol()->getNumBonds() == scsiTest->bondCount)
      TEST_ASSERT(scsrMol->getTemplateCount() == scsiTest->templateCount)

      RDKit::v2::FileParsers::MolFromScsrParams molFromScsrParams;
      molFromScsrParams.includeLeavingGroups = true;
      molFromScsrParams.scsrTemplateNames =
          RDKit::v2::FileParsers::ScsrTemplateNamesAsEntered;

      auto mol = MolFromSCSRMol(*scsrMol, molFromScsrParams);

      TEST_ASSERT(mol != nullptr)
      TEST_ASSERT(mol->getNumAtoms() == scsiTest->totalAtomCount)
      TEST_ASSERT(mol->getNumBonds() == scsiTest->totalBondCount)

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

        TEST_ASSERT(GetExpectedValue(expectedMolName) == outMolStr);
      }
      // now make the expanded mol in "query" mode - not including any leaving
      // groups

      molFromScsrParams.includeLeavingGroups = true;
      molFromScsrParams.scsrTemplateNames =
          RDKit::v2::FileParsers::ScsrTemplateNamesUseFirstName;
      auto mol2 = MolFromSCSRMol(*scsrMol, molFromScsrParams);

      // const std::unique_ptr<RDKit::RWMol> mol;
      TEST_ASSERT(mol2 != nullptr)
      TEST_ASSERT(mol2->getNumAtoms() == scsiTest->totalQueryAtomCount)
      TEST_ASSERT(mol2->getNumBonds() == scsiTest->totalQueryBondCount)

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

        TEST_ASSERT(GetExpectedValue(expectedMolName) == outMolStr);
      }
    } catch (const std::exception &e) {
      if (scsiTest->expectedResult != false) {
        throw;
      }
      return;
    }
    TEST_ASSERT(scsiTest->expectedResult == true);
  }

  void RunTests() {
    // the molecule tests

    if (testToRun == "" || testToRun == "scsiTests") {
      std::list<ScsiTest> scsiTests{
          ScsiTest("TestRNA2_fixed.mol", true, 15, 14, 5, 106, 117, 104, 115),
          ScsiTest("Mixed.mol", true, 14, 16, 5, 51, 54, 51, 54),
          ScsiTest("CrossLink.mol", true, 8, 8, 5, 47, 48, 45, 46),
          ScsiTest("cyclic.mol", true, 8, 9, 5, 45, 47, 45, 47),

          ScsiTest("Triplet.mol", true, 3, 2, 3, 30, 30, 27, 27),
          ScsiTest("FromBioviaDoc.mol", true, 5, 4, 3, 27, 26, 25, 24),
          ScsiTest("testScsr.mol", true, 8, 7, 6, 64, 66, 57, 59),
          ScsiTest("testScsr.mol", false, 1, 7, 6, 64, 66, 64, 66),
          ScsiTest("testScsr.mol", false, 8, 1, 6, 64, 66, 64, 66),
          ScsiTest("badAtomName.mol", false, 8, 7, 6, 0, 0, 0, 0),
          ScsiTest("badClass.mol", false, 8, 7, 6, 0, 0, 0, 0),
          ScsiTest("badClassTemplate.mol", false, 8, 7, 6, 0, 0, 0, 0),
          ScsiTest("badMissingTemplate.mol", false, 8, 7, 6, 0, 0, 0, 0),
          ScsiTest("obj3dTest.mol", true, 5, 4, 3, 27, 26, 25, 24),
          ScsiTest("obj3dTest2.mol", true, 5, 4, 3, 27, 26, 25, 24),
          ScsiTest("obj3dFoundTwice.mol", false, 5, 4, 3, 27, 26, 27, 26),
          ScsiTest("SgroupFoundTwice.mol", false, 5, 4, 3, 0, 0, 0, 0),
      };

      for (auto scsiTest : scsiTests) {
        BOOST_LOG(rdInfoLog) << "Test: " << scsiTest.fileName << std::endl;

        testScsiFiles(&scsiTest);
      }
    }

    if (testToRun == "" || testToRun == "threeLetterCodeTest") {
      std::list<ScsiTest> scsiTests{
          // ScsiTest("KellanError.mol", true, 31, 30, 15, 244, 263, 236, 255),
          ScsiTest("PepTla.mol", true, 4, 3, 4, 26, 25, 26, 25),

      };

      for (auto scsiTest : scsiTests) {
        BOOST_LOG(rdInfoLog) << "Test: " << scsiTest.fileName << std::endl;

        threeLetterCodeTest(&scsiTest);
      }
    }
  }
};

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  RDLog::InitLogs();
  boost::logging::enable_logs("rdApp.info");
  BOOST_LOG(rdInfoLog) << " ---- Running with POSIX locale ----- " << std::endl;

  ScsiMolTest scsiMolTest;

  if (argc > 1) {
    scsiMolTest.testToRun = argv[1];
  }

  if (argc > 2 && std::string(argv[2]) == "generate") {
    scsiMolTest.generateExpectedFiles = true;
  }

  RDLog::InitLogs();
  BOOST_LOG(rdInfoLog) << " ---- Running with POSIX locale ----- " << std::endl;

  scsiMolTest.RunTests();

  return 0;
}
