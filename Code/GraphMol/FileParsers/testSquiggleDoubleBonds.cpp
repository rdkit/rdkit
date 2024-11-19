//
//  Copyright (C) 2002-2021 Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RDGeneral/test.h"
#include <catch2/catch_all.hpp>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <RDGeneral/FileParseException.h>
#include <boost/algorithm/string.hpp>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/Chirality.h>

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>

#include <string>

using namespace RDKit;

class MolTest {
 public:
  std::string fileName;
  bool expectedResult;

  unsigned int atomCount;
  unsigned int bondCount;

  MolTest(std::string fileNameInit, bool expectedResultInit, int atomCountInit,
          int bondCountInit)
      : fileName(fileNameInit),
        expectedResult(expectedResultInit),
        atomCount(atomCountInit),
        bondCount(bondCountInit) {};
};

void testMolFiles(const MolTest *molFileTest) {
  BOOST_LOG(rdInfoLog) << "testing marvin writing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/" + molFileTest->fileName;

  try {
    std::unique_ptr<RWMol> mol(MolFileToMol(fName, true, false, false));
    RDKit::Chirality::reapplyMolBlockWedging(*mol);

    CHECK(mol != nullptr);
    CHECK(mol->getNumAtoms() == molFileTest->atomCount);
    CHECK(mol->getNumBonds() == molFileTest->bondCount);

    {
      std::string expectedMrvName = fName + ".expected.sdf";
      std::string outMolStr = "";
      try {
        outMolStr = MolToMolBlock(*mol, true, 0, true, true);
      } catch (const RDKit::KekulizeException &e) {
        outMolStr = "";
      } catch (...) {
        throw;  // re-trhow the error if not a kekule error
      }
      if (outMolStr == "") {
        outMolStr = MolToMolBlock(*mol, true, 0, false,
                                  true);  // try without kekule'ing
      }

      // code to create the expected files for new or changed tests

      // {
      //   std::ofstream out;
      //   out.open(fName + ".NEW.sdf");
      //   out << outMolStr;
      // }

      std::stringstream expectedMolStr;
      std::ifstream in;
      in.open(expectedMrvName);
      expectedMolStr << in.rdbuf();
      std::string expectedStr = expectedMolStr.str();

      CHECK(expectedStr == outMolStr);
    }

    BOOST_LOG(rdInfoLog) << "done" << std::endl;
  } catch (const std::exception &e) {
    if (molFileTest->expectedResult != false) {
      throw;
    }
    return;
  }

  CHECK(molFileTest->expectedResult == true);

  return;
}

TEST_CASE("CrossedDoubleBondWithChiralNbr") {
  SECTION("SIMPLE") {
    MolTest molTest("CrossedDoubleBondWithChiralNbr.sdf", true, 10, 9);
    testMolFiles(&molTest);
  }
}

TEST_CASE("CrossedDoubleBondWithChiralNbr2") {
  SECTION("SIMPLE") {
    MolTest molTest("CrossedDoubleBondWithChiralNbr2.sdf", true, 10, 9);
    testMolFiles(&molTest);
  }
}

TEST_CASE("SimpleWiggleDoubleBond") {
  SECTION("SIMPLE") {
    MolTest molTest("SimpleWiggleDoubleBond.sdf", true, 6, 5);
    testMolFiles(&molTest);
  }
}
