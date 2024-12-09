//
//  Copyright (C) 2002-2021 Collaboartive Drug Discovery and other RDKit
//  contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/SequenceParsers.h>
#include <GraphMol/FileParsers/SequenceWriters.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/Depictor/RDDepictor.h>

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Atropisomers.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/test_fixtures.h>

#include <string>
#include <fstream>
#include <boost/lexical_cast.hpp>
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
    bool expectedResult;

    ScsiTest(std::string fileNameInit, bool expectedResultInit,
             int atomCountInit, int bondCountInit, int templateCountInit)
        : atomCount(atomCountInit),
          bondCount(bondCountInit),
          fileName(fileNameInit),
          templateCount(templateCountInit),
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
    BOOST_LOG(rdInfoLog) << "testing mol files with atropisomers" << std::endl;

    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase +
                        "/Code/GraphMol/FileParsers/test_data/macromols/" +
                        scsiTest->fileName;

    try {
      RDKit::v2::FileParsers::MolFileParserParams pp;
      pp.sanitize = false;
      pp.removeHs = false;
      pp.strictParsing = false;

      std::unique_ptr<SCSRMol> scsiMol(ScsrMolFromScsrFile(fName, pp));
      RDKit::Chirality::removeNonExplicit3DChirality(*scsiMol->getMol());

      TEST_ASSERT(scsiMol != nullptr);
      TEST_ASSERT(scsiMol->getMol()->getNumAtoms() == scsiTest->atomCount)
      TEST_ASSERT(scsiMol->getMol()->getNumBonds() == scsiTest->bondCount)
      TEST_ASSERT(scsiMol->getTemplateCount() == scsiTest->templateCount)

      {
        // std::string outMolStr = MolToMolBlock(*mol, true, 0, true, true);

        // generateNewExpectedFilesIfSoSpecified(fName + ".NEW.sdf", outMolStr);

        // TEST_ASSERT(GetExpectedValue(fName + ".expected.sdf") == outMolStr);
      }

      return;
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
          ScsiTest("Triplet.mol", true, 3, 2, 3),
          ScsiTest("FromBioviaDoc.mol", true, 5, 4, 3),
          ScsiTest("testScsr.mol", true, 8, 7, 6),
          ScsiTest("testScsr.mol", false, 1, 7, 6),
          ScsiTest("testScsr.mol", false, 8, 1, 6),
          ScsiTest("badAtomName.mol", false, 8, 7, 6),
          ScsiTest("badClass.mol", false, 8, 7, 6),
          ScsiTest("badClassTemplate.mol", false, 8, 7, 6),
          ScsiTest("badMissingTemplate.mol", false, 8, 7, 6),
          ScsiTest("obj3dTest.mol", true, 5, 4, 3),
          ScsiTest("obj3dTest2.mol", true, 5, 4, 3),
          ScsiTest("obj3dFoundTwice.mol", false, 5, 4, 3),
          ScsiTest("SgroupFoundTwice.mol", false, 5, 4, 3),
      };

      for (auto scsiTest : scsiTests) {
        BOOST_LOG(rdInfoLog) << "Test: " << scsiTest.fileName << std::endl;

        testScsiFiles(&scsiTest);
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
