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

  class SmiTest {
   public:
    unsigned int atomCount;
    unsigned int bondCount;
    std::string smiles;
    std::string name;
    bool expectedResult;

    SmiTest(std::string smilesInit, std::string nameInit,
            bool expectedResultInit, int atomCountInit, int bondCountInit)
        : atomCount(atomCountInit),
          bondCount(bondCountInit),
          smiles(smilesInit),
          name(nameInit),
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

  void testAromAtropSmiles(const SmiTest *smiTest) {
    BOOST_LOG(rdInfoLog) << "testing aromtizing smiles with atropisomers"
                         << std::endl;

    std::string rdbase = getenv("RDBASE");
    std::string fBase = rdbase +
                        "/Code/GraphMol/FileParsers/test_data/atropisomers/" +
                        smiTest->name;

    UseLegacyStereoPerceptionFixture useLegacy(false);

    RDKit::SmilesParserParams params;
    params.allowCXSMILES = true;    // recognize and parse CXSMILES
    params.debugParse = false;      // disable debugging in the SMILES parser
    params.parseName = true;        // parse (and set) the molecule name as well
    params.removeHs = false;        // do not remove Hs
    params.replacements = nullptr;  // no SMILES "macros"
    params.sanitize = false;        // do not sanitize the molecule
    params.skipCleanup =
        false;  // do not skip the final cleanup stage (for internal RDKit use)
    params.strictCXSMILES =
        false;  // do not throw an exception if the CXSMILES parsing fails

    std::unique_ptr<RWMol> mol =
        std::unique_ptr<RWMol>(SmilesToMol(smiTest->smiles, params));
    TEST_ASSERT(mol->getNumAtoms() == smiTest->atomCount)
    TEST_ASSERT(mol->getNumBonds() == smiTest->bondCount)

    testAromAtropMol(mol.get(), smiTest->expectedResult, fBase);
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

      RDKit::Chirality::removeNonExplicit3DChirality(*mol);

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
      auto sanitizeOps = MolOps::SANITIZE_SETAROMATICITY +
                         MolOps::SANITIZE_CLEANUPCHIRALITY +
                         MolOps::SANITIZE_SETHYBRIDIZATIONFORCE +
                         MolOps::SANITIZE_SETCONJUGATION;

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

    if (testToRun == "" || testToRun == "AromAtropSmi") {
      std::list<SmiTest> sdfTests{
          SmiTest(
              "[H:0][O:0][C:0]([C:0]([H:0])([H:0])[H:0])([C:0]([H:0])([H:0])[H:0])[C:0]1([H:0])[C:0]([H:0])([H:0])[c:0]2[c:0]([c:0]3[c:0](-[c:0]4[c:0]([H:0])[c:0]([H:0])[c:0]([H:0])[c:0](-[n:0]5[c:0](=[O:0])[c:0]6[c:0]([H:0])[c:0]([H:0])[c:0]([H:0])[c:0]([F:0])[c:0]6[n:0]([C:0]([H:0])([H:0])[H:0])[c:0]5=[O:0])[c:0]4[C:0]([H:0])([H:0])[H:0])[c:0]([F:0])[c:0]([H:0])[c:0]([C:0](=[O:0])[N:0]([H:0])[H:0])[c:0]3[n:0]2[H:0])[C:0]([H:0])([H:0])[C:0]1([H:0])[H:0] |(-1.784,2.6885,-15353.5;-1.7706,2.244,-15353.7;-2.0401,2.2783,-15354.4;-1.6036,2.7864,-15354.8;-1.726,2.7756,-15355.4;-1.0502,2.7191,-15354.7;-1.6944,3.3158,-15354.6;-2.7841,2.5509,-15354.3;-3.0456,2.5184,-15354.8;-2.8004,3.0803,-15354.1;-3.0848,2.2449,-15353.9;-2.0126,1.5618,-15354.7;-2.2375,1.5925,-15355.2;-2.4426,1.041,-15354.3;-2.3286,1.1012,-15353.7;-2.9948,1.1437,-15354.3;-2.2836,0.3326,-15354.5;-1.7256,0.1306,-15354.9;-1.7725,-0.5951,-15355;-1.3692,-1.0929,-15355.3;-0.7571,-0.8884,-15355.6;-0.8155,-0.7282,-15356.4;-1.3113,-0.7619,-15356.7;-0.2325,-0.5275,-15356.7;-0.2779,-0.4082,-15357.3;0.4089,-0.514,-15356.4;0.8596,-0.3516,-15356.7;0.4671,-0.6745,-15355.7;1.1391,-0.6579,-15355.4;1.3334,-0.0217,-15355.1;1.0062,0.5319,-15355.1;2.0064,-0.0037,-15354.7;2.2121,0.5921,-15354.4;1.8803,1.0542,-15354.4;2.8519,0.6365,-15354.1;3.0094,1.0843,-15353.8;3.2878,0.0607,-15354.1;3.7881,0.0748,-15353.9;3.0825,-0.5334,-15354.4;3.5292,-1.064,-15354.4;2.437,-0.5836,-15354.8;2.2153,-1.1841,-15355.1;2.6813,-1.7797,-15355.1;2.7755,-1.9594,-15354.6;3.161,-1.6393,-15355.4;2.4579,-2.2139,-15355.4;1.5572,-1.2626,-15355.4;1.3554,-1.8171,-15355.7;-0.1157,-0.875,-15355.3;-0.0488,-1.067,-15354.6;-0.2399,-1.5804,-15354.5;0.4797,-1.0378,-15354.4;-0.341,-0.6899,-15354.3;-1.5711,-1.7842,-15355.3;-1.1945,-2.279,-15355.6;-2.1638,-1.9781,-15354.9;-2.2805,-2.5457,-15354.9;-2.5729,-1.4909,-15354.6;-3.1903,-1.7006,-15354.2;-3.5657,-1.2943,-15353.9;-3.3571,-2.413,-15354.2;-3.0833,-2.7864,-15354.4;-3.7881,-2.5608,-15353.9;-2.3657,-0.7991,-15354.6;-2.6659,-0.2442,-15354.4;-3.1024,-0.231,-15354.1;-1.2069,0.6208,-15355.2;-0.6821,0.4199,-15355.1;-1.2966,0.7025,-15355.7;-1.2568,1.3334,-15354.8;-1.0143,1.2664,-15354.3;-0.9475,1.6827,-15355.1),wD:11.11,28.47|",
              "BMS-986142_3d_chiral", true, 72, 77),
          SmiTest(
              "[H:0][O:0][C:0]([C:0]([H:0])([H:0])[H:0])([C:0]([H:0])([H:0])[H:0])[C:0]1([H:0])[C:0]([H:0])([H:0])[c:0]2[c:0]([c:0]3[c:0](-[c:0]4[c:0]([H:0])[c:0]([H:0])[c:0]([H:0])[c:0](-[n:0]5[c:0](=[O:0])[c:0]6[c:0]([H:0])[c:0]([H:0])[c:0]([H:0])[c:0]([F:0])[c:0]6[n:0]([C:0]([H:0])([H:0])[H:0])[c:0]5=[O:0])[c:0]4[C:0]([H:0])([H:0])[H:0])[c:0]([F:0])[c:0]([H:0])[c:0]([C:0](=[O:0])[N:0]([H:0])[H:0])[c:0]3[n:0]2[H:0])[C:0]([H:0])([H:0])[C:0]1([H:0])[H:0] |(-1.784,2.6885,-15353.5;-1.7706,2.244,-15353.7;-2.0401,2.2783,-15354.4;-1.6036,2.7864,-15354.8;-1.726,2.7756,-15355.4;-1.0502,2.7191,-15354.7;-1.6944,3.3158,-15354.6;-2.7841,2.5509,-15354.3;-3.0456,2.5184,-15354.8;-2.8004,3.0803,-15354.1;-3.0848,2.2449,-15353.9;-2.0126,1.5618,-15354.7;-2.2375,1.5925,-15355.2;-2.4426,1.041,-15354.3;-2.3286,1.1012,-15353.7;-2.9948,1.1437,-15354.3;-2.2836,0.3326,-15354.5;-1.7256,0.1306,-15354.9;-1.7725,-0.5951,-15355;-1.3692,-1.0929,-15355.3;-0.7571,-0.8884,-15355.6;-0.8155,-0.7282,-15356.4;-1.3113,-0.7619,-15356.7;-0.2325,-0.5275,-15356.7;-0.2779,-0.4082,-15357.3;0.4089,-0.514,-15356.4;0.8596,-0.3516,-15356.7;0.4671,-0.6745,-15355.7;1.1391,-0.6579,-15355.4;1.3334,-0.0217,-15355.1;1.0062,0.5319,-15355.1;2.0064,-0.0037,-15354.7;2.2121,0.5921,-15354.4;1.8803,1.0542,-15354.4;2.8519,0.6365,-15354.1;3.0094,1.0843,-15353.8;3.2878,0.0607,-15354.1;3.7881,0.0748,-15353.9;3.0825,-0.5334,-15354.4;3.5292,-1.064,-15354.4;2.437,-0.5836,-15354.8;2.2153,-1.1841,-15355.1;2.6813,-1.7797,-15355.1;2.7755,-1.9594,-15354.6;3.161,-1.6393,-15355.4;2.4579,-2.2139,-15355.4;1.5572,-1.2626,-15355.4;1.3554,-1.8171,-15355.7;-0.1157,-0.875,-15355.3;-0.0488,-1.067,-15354.6;-0.2399,-1.5804,-15354.5;0.4797,-1.0378,-15354.4;-0.341,-0.6899,-15354.3;-1.5711,-1.7842,-15355.3;-1.1945,-2.279,-15355.6;-2.1638,-1.9781,-15354.9;-2.2805,-2.5457,-15354.9;-2.5729,-1.4909,-15354.6;-3.1903,-1.7006,-15354.2;-3.5657,-1.2943,-15353.9;-3.3571,-2.413,-15354.2;-3.0833,-2.7864,-15354.4;-3.7881,-2.5608,-15353.9;-2.3657,-0.7991,-15354.6;-2.6659,-0.2442,-15354.4;-3.1024,-0.231,-15354.1;-1.2069,0.6208,-15355.2;-0.6821,0.4199,-15355.1;-1.2966,0.7025,-15355.7;-1.2568,1.3334,-15354.8;-1.0143,1.2664,-15354.3;-0.9475,1.6827,-15355.1),wD:28.47|", "BMS-986142_3d", true,
              72, 77),
      };

      for (auto sdfTest : sdfTests) {
        BOOST_LOG(rdInfoLog) << "Test: " << sdfTest.name << std::endl;

        testAromAtropSmiles(&sdfTest);
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
