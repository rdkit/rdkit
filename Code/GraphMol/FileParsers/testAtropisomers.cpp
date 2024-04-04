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

#include <string>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <filesystem>
using namespace RDKit;

class MolAtropTest {
 public:
 public:
  std::string testToRun;
  bool generateExpectedFiles;

  MolAtropTest() {
    testToRun = "";
    generateExpectedFiles = false;
  }

 public:
  class MolTest {
   public:
    unsigned int atomCount;
    unsigned int bondCount;
    std::string fileName;
    bool expectedResult;

    MolTest(std::string fileNameInit, bool expectedResultInit,
            int atomCountInit, int bondCountInit)
        : atomCount(atomCountInit),
          bondCount(bondCountInit),
          fileName(fileNameInit),
          expectedResult(expectedResultInit){};
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

  void testMolFiles(const MolTest *molFileTest) {
    BOOST_LOG(rdInfoLog) << "testing mol files with atropisomers" << std::endl;

    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase +
                        "/Code/GraphMol/FileParsers/test_data/atropisomers/" +
                        molFileTest->fileName;

    // first pass - no sanitiaton, but reapplyMolBlockWedging is ON

    try {
      std::unique_ptr<RWMol> mol(MolFileToMol(fName, false, false, false));
      RDKit::Chirality::removeNonExplicit3DChirality(*mol);

      TEST_ASSERT(mol != nullptr);
      TEST_ASSERT(mol->getNumAtoms() == molFileTest->atomCount)
      TEST_ASSERT(mol->getNumBonds() == molFileTest->bondCount)

      {
        MolOps::Kekulize(*mol);
        RDKit::Chirality::reapplyMolBlockWedging(*mol);
        std::string outMolStr = MolToMolBlock(*mol, true, 0, true, true);

        generateNewExpectedFilesIfSoSpecified(fName + ".NEW.sdf", outMolStr);

        TEST_ASSERT(GetExpectedValue(fName + ".expected.sdf") == outMolStr);
      }

      // 2nd pass without reapplying the mol block wedging -
      //     the itropisomers will be marked automatically
      // SANITIZATION IS ON

      mol = std::unique_ptr<RWMol>(MolFileToMol(fName, true, false, false));
      RDKit::Chirality::removeNonExplicit3DChirality(*mol);

      {
        MolOps::Kekulize(*mol);
        std::string outMolStr = MolToMolBlock(*mol, true, 0, true, true);

        generateNewExpectedFilesIfSoSpecified(fName + ".NEW2.sdf", outMolStr);

        TEST_ASSERT(GetExpectedValue(fName + ".expected2.sdf") == outMolStr);
      }

      // CXSMILES

      mol = std::unique_ptr<RWMol>(MolFileToMol(fName, false, false, false));
      RDKit::Chirality::removeNonExplicit3DChirality(*mol);

      {
        std::string expectedFileName = fName + ".expected.cxsmi";
        SmilesWriteParams ps;
        ps.canonical = false;
        ps.doKekule = true;

        unsigned int flags = SmilesWrite::CXSmilesFields::CX_COORDS |
                             SmilesWrite::CXSmilesFields::CX_MOLFILE_VALUES |
                             SmilesWrite::CXSmilesFields::CX_ATOM_PROPS |
                             SmilesWrite::CXSmilesFields::CX_BOND_CFG |
                             SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO
            //| SmilesWrite::CXSmilesFields::CX_ALL
            ;

        std::string smilesOut =
            MolToCXSmiles(*mol, ps, flags, RestoreBondDirOptionTrue);

        generateNewExpectedFilesIfSoSpecified(fName + ".NEW.cxsmi", smilesOut);

        TEST_ASSERT(GetExpectedValue(fName + ".expected.cxsmi") == smilesOut);
      }

      BOOST_LOG(rdInfoLog) << "done" << std::endl;
    } catch (const std::exception &e) {
      if (molFileTest->expectedResult != false) {
        throw;
      }
      return;
    }
    TEST_ASSERT(molFileTest->expectedResult == true);

    return;
  }

  void RunTests() {
    // the molecule tests

    if (testToRun == "" || testToRun == "sdfTests") {
      std::list<MolTest> sdfTests{
          MolTest("atropWedgeTest.sdf", true, 16, 16),
          MolTest("AtropTest.sdf", true, 38, 41),
          MolTest("AtropManyChiralsEnhanced.sdf", true, 20, 20),
          MolTest("AtropManyChiralsEnhanced2.sdf", true, 20, 20),
          MolTest("AtropManyChirals.sdf", true, 20, 20),
          MolTest("BMS-986142.sdf", true, 42, 47),
          MolTest("BMS-986142_3d_chiral.sdf", true, 72, 77),
          MolTest("BMS-986142_3d.sdf", true, 72, 77),
          MolTest("BMS-986142_atrop1.sdf", true, 42, 47),
          MolTest("BMS-986142_atrop2.sdf", true, 42, 47),
          MolTest("BMS-986142_atrop3.sdf", true, 42, 47),
          MolTest("BMS-986142_atrop4.sdf", true, 42, 47),
          MolTest("BMS-986142_atrop5.sdf", true, 42, 47),
          MolTest("BMS-986142_atrop6.sdf", true, 42, 47),
          MolTest("BMS-986142_atrop7.sdf", true, 42, 47),
          MolTest("BMS-986142_atrop8.sdf", true, 42, 47),
          MolTest("BMS-986142_atropBad2.sdf", true, 42, 47),
          MolTest("JDQ443.sdf", true, 38, 44),
          MolTest("JDQ443_3d.sdf", true, 66, 72),
          MolTest("JDQ443_atrop1.sdf", true, 38, 44),
          MolTest("JDQ443_atrop2.sdf", true, 38, 44),
          MolTest("JDQ443_atrop3.sdf", true, 38, 44),
          MolTest("JDQ443_atropBad1.sdf", true, 38, 44),
          MolTest("RP-6306.sdf", true, 24, 26),
          MolTest("RP-6306_atrop1.sdf", true, 24, 26),
          MolTest("RP-6306_atrop2.sdf", true, 24, 26),
          MolTest("RP-6306_atrop3.sdf", true, 24, 26),
          MolTest("RP-6306_atrop4.sdf", true, 24, 26),
          MolTest("RP-6306_atrop5.sdf", true, 24, 26),
          MolTest("RP-6306_atropBad1.sdf", true, 24, 26),
          MolTest("RP-6306_atropBad2.sdf", true, 24, 26),
          // note the rp-6306_3d.sdf is backwards from the 2D versions
          // the 2D version were based on images from drug hunter
          // the 3D version came from PUBCHEM
          MolTest("RP-6306_3d.sdf", true, 44, 46),
          MolTest("Sotorasib.sdf", true, 41, 45),
          MolTest("Sotorasib_atrop1.sdf", true, 41, 45),
          MolTest("Sotorasib_atrop2.sdf", true, 41, 45),
          MolTest("Sotorasib_atrop3.sdf", true, 41, 45),
          MolTest("Sotorasib_atrop4.sdf", true, 41, 45),
          MolTest("Sotorasib_atrop5.sdf", true, 41, 45),
          MolTest("Sotorasib_atropBad1.sdf", true, 41, 45),
          MolTest("Sotorasib_atropBad2.sdf", true, 41, 45),
          // note the sotorasib_3d.sdf is backwards from the 2D versions
          // the 2D version were based on images from drug hunter
          // the 3D version came from PUBCHEM
          MolTest("Sotorasib_3d.sdf", true, 71, 75),
          MolTest("ZM374979.sdf", true, 45, 49),
          MolTest("ZM374979_atrop1.sdf", true, 45, 49),
          MolTest("ZM374979_atrop2.sdf", true, 45, 49),
          MolTest("ZM374979_atrop3.sdf", true, 45, 49),
          MolTest("ZM374979_atropBad1.sdf", true, 45, 49),
          MolTest("Mrtx1719.sdf", true, 33, 37),
          // note the Mrtx1719_3d.sdf is backwards from the 2D versions
          // the 2D version were based on images from drug hunter
          // the 3D version came from PUBCHEM
          MolTest("Mrtx1719_3d.sdf", true, 51, 55),
          MolTest("Mrtx1719_atrop1.sdf", true, 33, 37),
          MolTest("Mrtx1719_atrop2.sdf", true, 33, 37),
          MolTest("Mrtx1719_atrop3.sdf", true, 33, 37),
          MolTest("Mrtx1719_atropBad1.sdf", true, 33, 37),
      };

      for (auto sdfTest : sdfTests) {
        BOOST_LOG(rdInfoLog) << "Test: " << sdfTest.fileName << std::endl;

        testMolFiles(&sdfTest);
      }
    }
  }
};

void testLookForAtropisomersInSDdfFiles(std::string fileName,
                                        unsigned int expectedHits,
                                        unsigned int expectedMisses) {
  BOOST_LOG(rdInfoLog) << "Looking for atropisomers in " << fileName
                       << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/atropisomers/" + fileName;

  std::ifstream in;
  in.open(fName);
  std::string line;
  unsigned int foundCount = 0;
  unsigned int notFoundCount = 0;
  while (!in.eof()) {
    std::string molBlock = "";
    while (std::getline(in, line)) {
      if (line.find("$$$$") != std::string::npos) {
        break;
      }

      molBlock += line + "\n";
    }

    if (molBlock.length() < 10) {
      continue;  // try for another;
    }

    std::unique_ptr<RWMol> mol(MolBlockToMol(molBlock, false, false, false));
    TEST_ASSERT(mol != nullptr);

    auto hasAtropisomers = RDKit::Atropisomers::doesMolHaveAtropisomers(*mol);

    if (hasAtropisomers) {
      BOOST_LOG(rdInfoLog) << "Found atropisomers in " << fileName << std::endl;
      foundCount++;
      printf("Atropisomers- %d hits   %d misses\r", foundCount, notFoundCount);
      std::flush(std::cout);
      std::ofstream out;
      out.open(fName + "_" + std::to_string(foundCount) + ".sdf");
      out << molBlock << std::endl;
    } else {
      notFoundCount++;
      if (notFoundCount % 100 == 0) {
        printf("Atropisomers- %d hits   %d misses\r", foundCount,
               notFoundCount);
        std::flush(std::cout);
      }
    }
  }
  printf("\nFinal results:\nFound atropisomers in %s - %d hits   %d misses\n",
         fileName.c_str(), foundCount, notFoundCount);

  TEST_ASSERT(foundCount == expectedHits);
  TEST_ASSERT(notFoundCount == expectedMisses);
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  RDLog::InitLogs();
  boost::logging::enable_logs("rdApp.info");
  BOOST_LOG(rdInfoLog) << " ---- Running with POSIX locale ----- " << std::endl;

  MolAtropTest molAtropTest;

  if (argc > 1) {
    molAtropTest.testToRun = argv[1];
  }

  if (argc > 2 && std::string(argv[2]) == "generate") {
    molAtropTest.generateExpectedFiles = true;
  }

  RDLog::InitLogs();
  BOOST_LOG(rdInfoLog) << " ---- Running with POSIX locale ----- " << std::endl;

  molAtropTest.RunTests();
  testLookForAtropisomersInSDdfFiles("TestMultInSDF.sdf", 1, 4);

  return 0;
}
