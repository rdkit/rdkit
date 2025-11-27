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
#include <sstream>
#include <vector>
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
          expectedResult(expectedResultInit) {};
  };

  class KekuleTest {
   public:
    unsigned int atomCount;
    unsigned int bondCount;
    std::string smiles;
    std::string expectedOutput;
    bool expectedResult;

    KekuleTest(std::string smilesInit, std::string nameInit,
               bool expectedResultInit, int atomCountInit, int bondCountInit)
        : atomCount(atomCountInit),
          bondCount(bondCountInit),
          smiles(smilesInit),
          expectedOutput(nameInit),
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

  void testAromAtropMolFile(const MolTest *molFileTest) {
    BOOST_LOG(rdInfoLog) << "testing aromtizing mol files with atropisomers"
                         << std::endl;

    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase +
                        "/Code/GraphMol/FileParsers/test_data/atropisomers/" +
                        molFileTest->fileName;

    UseLegacyStereoPerceptionFixture useLegacy(false);

    std::unique_ptr<RWMol> mol =
        std::unique_ptr<RWMol>(MolFileToMol(fName, false, false, false));

    TEST_ASSERT(mol->getNumAtoms() == molFileTest->atomCount)
    TEST_ASSERT(mol->getNumBonds() == molFileTest->bondCount)

    testAromAtropMol(mol.get(), molFileTest->expectedResult, fName);
  }

  void testAromAtropMol(RWMol *mol, bool expectedResult, std::string fName) {
    try {
      RDKit::Chirality::removeNonExplicit3DChirality(*mol);

      auto sanitizeOps =
          MolOps::SANITIZE_SETAROMATICITY + MolOps::SANITIZE_CLEANUPCHIRALITY +
          MolOps::SANITIZE_SETHYBRIDIZATION + MolOps::SANITIZE_SETCONJUGATION;

      unsigned int operationThatFailed = 0;
      RDKit::MolOps::sanitizeMol(*mol, operationThatFailed, sanitizeOps);
      for (const auto &atom : mol->atoms()) {
        if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
          BOOST_LOG(rdInfoLog) << "atom with stereo" << std::endl;
        }
      }
      for (const auto &bond : mol->bonds()) {
        if (bond->getBondDir() != Bond::BondDir::NONE) {
          BOOST_LOG(rdInfoLog) << "bond with wedging" << std::endl;
        }
      }
      {
        SmilesWriteParams ps;
        ps.canonical = true;
        ps.doKekule = false;
        ps.allBondsExplicit = false;
        ps.allHsExplicit = false;
        ps.doIsomericSmiles = true;
        ps.doKekule = false;
        ps.doRandom = false;
        ps.rootedAtAtom = -1;

        unsigned int flags = SmilesWrite::CXSmilesFields::CX_COORDS |
                             SmilesWrite::CXSmilesFields::CX_MOLFILE_VALUES |
                             SmilesWrite::CXSmilesFields::CX_ATOM_PROPS |
                             SmilesWrite::CXSmilesFields::CX_BOND_CFG |
                             SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO;

        std::string smilesOut =
            MolToCXSmiles(*mol, ps, flags, RestoreBondDirOptionClear);
        std::string expectedFileName = fName + ".expectedArom.cxsmi";

        generateNewExpectedFilesIfSoSpecified(fName + ".NEWArom.cxsmi",
                                              smilesOut);

        TEST_ASSERT(GetExpectedValue(fName + ".expectedArom.cxsmi") ==
                    smilesOut);
      }

      BOOST_LOG(rdInfoLog) << "done" << std::endl;
    } catch (const std::exception &e) {
      if (expectedResult != false) {
        throw;
      }
      return;
    }
    TEST_ASSERT(expectedResult == true);
  }

  void testKekuleWedgeError(RWMol *mol, std::string expectedSmi,
                            bool expectedResult, unsigned int expectedAtomCount,
                            unsigned int expectedBondCount) {
    BOOST_LOG(rdInfoLog) << "testing aromatic atropisomers" << std::endl;

    TEST_ASSERT(mol != nullptr);
    TEST_ASSERT(mol->getNumAtoms() == expectedAtomCount)
    TEST_ASSERT(mol->getNumBonds() == expectedBondCount)

    try {
      SmilesWriteParams ps;
      ps.canonical = true;
      ps.doKekule = true;
      ps.allBondsExplicit = false;
      ps.allHsExplicit = false;
      ps.doIsomericSmiles = true;
      ps.doRandom = false;
      ps.rootedAtAtom = -1;

      unsigned int flags = SmilesWrite::CXSmilesFields::CX_MOLFILE_VALUES |
                           SmilesWrite::CXSmilesFields::CX_ATOM_PROPS |
                           SmilesWrite::CXSmilesFields::CX_BOND_CFG |
                           SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO;

      std::string smilesOut =
          MolToCXSmiles(*mol, ps, flags, RestoreBondDirOptionClear);
      TEST_ASSERT(expectedSmi == smilesOut);

    } catch (const std::exception &e) {
      if (expectedResult != false) {
        throw;
      }
      return;
    }
    TEST_ASSERT(expectedResult == true);

    return;
  }

  void testKekuleWedgeErrorSmiles(const KekuleTest *kekuleTest) {
    BOOST_LOG(rdInfoLog) << "testing aromatic atropisomers" << std::endl;

    UseLegacyStereoPerceptionFixture useLegacy(false);

    RDKit::SmilesParserParams params;
    params.allowCXSMILES = true;    // recognize and parse CXSMILES
    params.debugParse = false;      // disable debugging in the SMILES parser
    params.parseName = true;        // parse (and set) the molecule name as well
    params.removeHs = false;        // do not remove Hs
    params.replacements = nullptr;  // no SMILES "macros"
    params.sanitize = true;         // do not sanitize the molecule
    params.skipCleanup =
        false;  // do not skip the final cleanup stage (for internal RDKit use)
    params.strictCXSMILES =
        false;  // do not throw an exception if the CXSMILES parsing fails

    std::unique_ptr<RWMol> mol =
        std::unique_ptr<RWMol>(SmilesToMol(kekuleTest->smiles, params));

    testKekuleWedgeError(mol.get(), kekuleTest->expectedOutput,
                         kekuleTest->expectedResult, kekuleTest->atomCount,
                         kekuleTest->bondCount);
  }

  void testKekuleWedgeErrorMol(const MolTest *kekuleTest) {
    BOOST_LOG(rdInfoLog) << "testing aromatic atropisomers" << std::endl;

    UseLegacyStereoPerceptionFixture useLegacy(false);
    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase +
                        "/Code/GraphMol/FileParsers/test_data/atropisomers/" +
                        kekuleTest->fileName;

    std::unique_ptr<RWMol> mol(MolFileToMol(fName, true, false, false));

    testKekuleWedgeError(mol.get(), GetExpectedValue(fName + ".expected.cxsmi"),
                         kekuleTest->expectedResult, kekuleTest->atomCount,
                         kekuleTest->bondCount);
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
          // macrocycles
          MolTest("macrocycle-9-meta-wedge.sdf", true, 24, 26),
          MolTest("macrocycle-9-ortho-wedge.sdf", true, 24, 26),
          MolTest("macrocycle-8-meta-wedge.sdf", true, 23, 25),
          MolTest("macrocycle-8-ortho-wedge.sdf", true, 23, 25),
          MolTest("macrocycle-7-meta-wedge.sdf", true, 22, 24),
          MolTest("macrocycle-7-ortho-wedge.sdf", true, 22, 24),
          MolTest("macrocycle-6-meta-wedge.sdf", true, 21, 23),
          MolTest("macrocycle-6-ortho-wedge.sdf", true, 21, 23),
          MolTest("macrocycle-5-meta-wedge.sdf", true, 20, 22),
          MolTest("macrocycle-5-ortho-wedge.sdf", true, 20, 22),
          MolTest("macrocycle-9-meta-Cl-ortho-wedge.sdf", true, 25, 27),
          MolTest("macrocycle-8-meta-Cl-ortho-wedge.sdf", true, 24, 26),
          MolTest("macrocycle-7-meta-Cl-ortho-wedge.sdf", true, 23, 25),
          MolTest("macrocycle-6-meta-Cl-ortho-wedge.sdf", true, 22, 24),
          MolTest("macrocycle-5-meta-Cl-ortho-wedge.sdf", true, 21, 23),
          MolTest("macrocycle-9-meta-broken-wedge.sdf", true, 24, 25),
          MolTest("macrocycle-9-ortho-broken-wedge.sdf", true, 24, 25),
          MolTest("macrocycle-8-meta-broken-wedge.sdf", true, 23, 24),
          MolTest("macrocycle-8-ortho-broken-wedge.sdf", true, 23, 24),
          MolTest("macrocycle-7-meta-broken-wedge.sdf", true, 22, 23),
          MolTest("macrocycle-7-ortho-broken-wedge.sdf", true, 22, 23),
          MolTest("macrocycle-6-meta-broken-wedge.sdf", true, 21, 22),
          MolTest("macrocycle-6-ortho-broken-wedge.sdf", true, 21, 22),
          MolTest("macrocycle-5-meta-broken-wedge.sdf", true, 20, 21),
          MolTest("macrocycle-5-ortho-broken-wedge.sdf", true, 20, 21),
          MolTest("macrocycle-9-meta-hash.sdf", true, 24, 26),
          MolTest("macrocycle-9-ortho-hash.sdf", true, 24, 26),
          MolTest("macrocycle-8-meta-hash.sdf", true, 23, 25),
          MolTest("macrocycle-8-ortho-hash.sdf", true, 23, 25),
          MolTest("macrocycle-7-meta-hash.sdf", true, 22, 24),
          MolTest("macrocycle-7-ortho-hash.sdf", true, 22, 24),
          MolTest("macrocycle-6-meta-hash.sdf", true, 21, 23),
          MolTest("macrocycle-6-ortho-hash.sdf", true, 21, 23),
          MolTest("macrocycle-5-meta-hash.sdf", true, 20, 22),
          MolTest("macrocycle-5-ortho-hash.sdf", true, 20, 22),
          MolTest("macrocycle-9-meta-Cl-ortho-hash.sdf", true, 25, 27),
          MolTest("macrocycle-8-meta-Cl-ortho-hash.sdf", true, 24, 26),
          MolTest("macrocycle-7-meta-Cl-ortho-hash.sdf", true, 23, 25),
          MolTest("macrocycle-6-meta-Cl-ortho-hash.sdf", true, 22, 24),
          MolTest("macrocycle-5-meta-Cl-ortho-hash.sdf", true, 21, 23),
          MolTest("macrocycle-9-meta-broken-hash.sdf", true, 24, 25),
          MolTest("macrocycle-9-ortho-broken-hash.sdf", true, 24, 25),
          MolTest("macrocycle-8-meta-broken-hash.sdf", true, 23, 24),
          MolTest("macrocycle-8-ortho-broken-hash.sdf", true, 23, 24),
          MolTest("macrocycle-7-meta-broken-hash.sdf", true, 22, 23),
          MolTest("macrocycle-7-ortho-broken-hash.sdf", true, 22, 23),
          MolTest("macrocycle-6-meta-broken-hash.sdf", true, 21, 22),
          MolTest("macrocycle-6-ortho-broken-hash.sdf", true, 21, 22),
          MolTest("macrocycle-5-meta-broken-hash.sdf", true, 20, 21),
          MolTest("macrocycle-5-ortho-broken-hash.sdf", true, 20, 21),
      };

      for (auto sdfTest : sdfTests) {
        BOOST_LOG(rdInfoLog) << "Test: " << sdfTest.fileName << std::endl;

        testMolFiles(&sdfTest);
      }
    }

    if (testToRun == "" || testToRun == "AromAtropMol") {
      std::list<MolTest> sdfTests{
          MolTest("BMS-986142_atrop1.sdf", true, 42, 47),
          MolTest("BMS-986142_3d_chiral.sdf", true, 72, 77),
          MolTest("BMS-986142_3d.sdf", true, 72, 77),
      };

      for (auto sdfTest : sdfTests) {
        BOOST_LOG(rdInfoLog) << "Test: " << sdfTest.fileName << std::endl;

        testAromAtropMolFile(&sdfTest);
      }
    }

    if (testToRun == "" || testToRun == "KekuleWedgeError") {
      std::list<KekuleTest> kekuleTests{
          KekuleTest(
              "CC1C(C2C(Cl)=CC=CC=2C)=C(Cl)C=CC=1 |wU:3.3,(17.57,-4.18,;17.57,-5.18,;16.71,-5.67,;15.84,-5.18,;15.84,-4.18,;16.71,-3.67,;14.98,-3.67,;14.11,-4.18,;14.11,-5.18,;14.98,-5.67,;14.98,-6.67,;16.71,-6.67,;15.84,-7.18,;17.57,-7.18,;18.44,-6.67,;18.44,-5.67,)|",
              "CC1=C(C2=C(C)C=CC=C2Cl)C(Cl)=CC=C1 |wU:2.11|", true, 16, 17),
          KekuleTest("CC1C(C2C(Cl)=CC=CC=2C)=C(Cl)C=CC=1 |wD:3.3|",
                     "CC1=C(C2=C(C)C=CC=C2Cl)C(Cl)=CC=C1 |wU:2.11|", true, 16,
                     17),
          KekuleTest(
              "CC1C=CC(O)=C(C)C=1N1C2C(=NC=CC=2OC2C=NC=CC=2)C(C(=O)N)=C1N |wD:8.7|",
              "CC1=C(O)C=CC(C)=C1N1C(N)=C(C(N)=O)C2=NC=CC(OC3=CN=CC=C3)=C21 |wD:9.10|",
              true, 29, 32),
      };

      for (auto kekuleTest : kekuleTests) {
        BOOST_LOG(rdInfoLog) << "Test: " << kekuleTest.smiles << std::endl;

        testKekuleWedgeErrorSmiles(&kekuleTest);
      }
    }

    if (testToRun == "" || testToRun == "KekuleWedgeErrorMol") {
      std::list<MolTest> kekuleTests{
          MolTest("atropWedgeError2.mol", true, 16, 17),
          MolTest("atropWedgeError.mol", true, 29, 32),
      };

      for (auto kekuleTest : kekuleTests) {
        BOOST_LOG(rdInfoLog) << "Test: " << kekuleTest.fileName << std::endl;

        testKekuleWedgeErrorMol(&kekuleTest);
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

void testSulfinamideExamplesHaveNoAtropisomers() {
  const std::vector<std::string> controlFiles = {
      "sulfinamide-double-bond-O-R.mol",
      "sulfinamide-single-bond-O-R.mol",
      "sulfinamide-single-bond-O-S.mol",
  };
  std::string rdbase = getenv("RDBASE");
  for (const auto &file : controlFiles) {
    auto fName = rdbase +
                 "/Code/GraphMol/FileParsers/test_data/atropisomers/" +
                 file;
    BOOST_LOG(rdInfoLog) << "Validating absence of atropisomers in " << file
                         << std::endl;
    auto mol = std::unique_ptr<RWMol>(MolFileToMol(fName, true, false, false));
    TEST_ASSERT(mol);
    const Conformer *conf =
        mol->getNumConformers() ? &mol->getConformer() : nullptr;
    Atropisomers::detectAtropisomerChirality(*mol, conf);
    TEST_ASSERT(!Atropisomers::doesMolHaveAtropisomers(*mol));
  }
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
  testSulfinamideExamplesHaveNoAtropisomers();

  return 0;
}
