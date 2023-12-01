//
//  Copyright (C) 2022-2023 Tad Hurst, Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/SequenceParsers.h>
#include <GraphMol/FileParsers/SequenceWriters.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include "MarvinParser.h"

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <string>
#include <fstream>
#include <filesystem>

using namespace RDKit;

class MrvTests {
 public:
  std::string testToRun;
  bool generateExpectedFiles;

  MrvTests() {
    testToRun = "";
    generateExpectedFiles = false;
  }

  class MolTest {
   public:
    unsigned int atomCount;
    unsigned int bondCount;
    std::string fileName;
    bool expectedResult;
    bool sanitizeFlag;
    bool reapplyMolBlockWedging;

    MolTest(std::string fileNameInit, bool expectedResultInit,
            int atomCountInit, int bondCountInit, bool sanitizeFlagInit = true,
            bool reapplyMolBlockWedgingInit = true)
        : atomCount(atomCountInit),
          bondCount(bondCountInit),
          fileName(fileNameInit),
          expectedResult(expectedResultInit),
          sanitizeFlag(sanitizeFlagInit),
          reapplyMolBlockWedging(reapplyMolBlockWedgingInit){};
  };

  class RxnTest {
   public:
    std::string fileName;
    bool expectedResult;
    unsigned int reactantCount;
    unsigned int agentCount;
    unsigned int productCount;
    unsigned int warnings;
    unsigned int errors;

    RxnTest(std::string fileNameInit, bool expectedResultInit,
            int reactantCountInit, int agentCountInit, int productCountInit,
            int warnInit, int errorInit)
        : fileName(fileNameInit),
          expectedResult(expectedResultInit),
          reactantCount(reactantCountInit),
          agentCount(agentCountInit),
          productCount(productCountInit),
          warnings(warnInit),
          errors(errorInit){};
  };

  class SmilesTest {
   public:
    std::string name;
    std::string smiles;
    bool expectedResult;
    bool sanitizeFlag;
    unsigned int atomCount;
    unsigned int bondCount;

    SmilesTest(std::string nameInit, std::string smilesInit,
               bool expectedResultInit, int atomCountInit, int bondCountInit,
               bool sanitizeFlagInit)
        : name(nameInit),
          smiles(smilesInit),
          expectedResult(expectedResultInit),
          sanitizeFlag(sanitizeFlagInit),
          atomCount(atomCountInit),
          bondCount(bondCountInit){};

    SmilesTest(std::string nameInit, std::string smilesInit,
               bool expectedResultInit, int atomCountInit, int bondCountInit)
        : name(nameInit),
          smiles(smilesInit),
          expectedResult(expectedResultInit),
          sanitizeFlag(true),
          atomCount(atomCountInit),
          bondCount(bondCountInit){};
  };

  RWMol *GetMol(const MolTest *molTest) {
    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/MarvinParse/test_data/" + molTest->fileName;

    try {
      return MrvFileToMol(fName, molTest->sanitizeFlag, false);
    } catch (const std::exception &e) {
      std::cerr << e.what() << '\n';
      throw BadFileException("Could not parse the MRV block");
    }
  }

  ChemicalReaction *GetReaction(const RxnTest *rxnTest) {
    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/MarvinParse/test_data/" + rxnTest->fileName;

    try {
      return MrvFileToChemicalReaction(fName, true, false);
    } catch (const std::exception &e) {
      std::cerr << e.what() << '\n';
    }
    throw BadFileException("Could not parse the MRV block");
  }

  std::string GetExpectedValue(std::string expectedFileName) {
    std::stringstream expectedMolStr;
    std::ifstream in;
    in.open(expectedFileName);
    expectedMolStr << in.rdbuf();
    return expectedMolStr.str();
  }

  void generateNewExpectedFilesIfSoSpecified(std::string filename,
                                             std::string dataToWrite) {
    if (generateExpectedFiles) {
      std::ofstream out;
      out.open(filename);
      out << dataToWrite;
    }
  }

  void testSmilesToMarvin(const SmilesTest *smilesTest) {
    BOOST_LOG(rdInfoLog) << "testing smiles to marvin " << std::endl;
    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/MarvinParse/test_data/" + smilesTest->name;

    class LocalVars  // protect against mem leak on error
    {
     public:
      RWMol *smilesMol = nullptr;

      LocalVars(){};

      ~LocalVars() { delete smilesMol; }
    } localVars;

    try {
      SmilesParserParams smilesParserParams;
      smilesParserParams.sanitize = smilesTest->sanitizeFlag;
      smilesParserParams.removeHs = smilesTest->sanitizeFlag;

      localVars.smilesMol = SmilesToMol(smilesTest->smiles, smilesParserParams);
      Chirality::reapplyMolBlockWedging(*localVars.smilesMol);

      TEST_ASSERT(localVars.smilesMol->getNumAtoms() == smilesTest->atomCount);
      TEST_ASSERT(localVars.smilesMol->getNumBonds() == smilesTest->bondCount);

      // test round trip back to smiles
      {
        std::string expectedMrvName =
            fName + (smilesTest->sanitizeFlag ? "" : ".nosan") +
            ".expected.smi";

        SmilesWriteParams ps;
        ps.canonical = false;

        std::string smilesOut = MolToSmiles(*localVars.smilesMol, ps);

        generateNewExpectedFilesIfSoSpecified(
            fName + (smilesTest->sanitizeFlag ? "" : ".nosan") + ".NEW.smi",
            smilesOut);

        TEST_ASSERT(GetExpectedValue(expectedMrvName) == smilesOut);
      }

      {
        std::string expectedMrvName =
            fName + (smilesTest->sanitizeFlag ? "" : ".nosan") +
            ".expected.sdf";
        std::string outMolStr = "";
        try {
          outMolStr = MolToMolBlock(*localVars.smilesMol, true, 0, true, true);
        } catch (const RDKit::KekulizeException &e) {
          outMolStr = "";
        }
        if (outMolStr == "") {
          outMolStr = MolToMolBlock(*localVars.smilesMol, true, 0, false,
                                    true);  // try without kekule'ing
        }

        generateNewExpectedFilesIfSoSpecified(
            fName + (smilesTest->sanitizeFlag ? "" : ".nosan") + ".NEW.sdf",
            outMolStr);

        TEST_ASSERT(GetExpectedValue(expectedMrvName) == outMolStr);
      }
      {
        std::string mrvBlock;
        std::string expectedMrvName =
            fName + (smilesTest->sanitizeFlag ? "" : ".nosan") +
            ".expected.mrv";
        std::string outMolStr = "";
        try {
          outMolStr =
              MolToMrvBlock(*localVars.smilesMol, true, -1, true, false);
        } catch (const RDKit::KekulizeException &e) {
          outMolStr = "";
        } catch (...) {
          throw;  // re-throw the error if not a kekule error
        }
        if (outMolStr == "") {
          outMolStr = MolToMrvBlock(*localVars.smilesMol, true, -1, false,
                                    false);  // try without kekule'ing
        }

        generateNewExpectedFilesIfSoSpecified(
            fName + (smilesTest->sanitizeFlag ? "" : ".nosan") + ".NEW.mrv",
            outMolStr);

        TEST_ASSERT(GetExpectedValue(expectedMrvName) == outMolStr);
      }
      BOOST_LOG(rdInfoLog) << "done" << std::endl;
    } catch (const std::exception &e) {
      if (smilesTest->expectedResult != false) {
        throw;
      }
      return;
    }

    TEST_ASSERT(smilesTest->expectedResult == true);
  }

  void testMarvinMol(const MolTest *molTest) {
    BOOST_LOG(rdInfoLog) << "testing marvin parsing" << std::endl;

    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/MarvinParse/test_data/" + molTest->fileName;

    class LocalVars  // protect against mem leak on error
    {
     public:
      RWMol *mol = nullptr;

      LocalVars(){};

      ~LocalVars() { delete (RWMol *)mol; }
    } localVars;

    try {
      if (MrvFileIsReaction(fName)) {
        TEST_ASSERT(molTest->expectedResult == false);
        return;
      }

      localVars.mol = GetMol(molTest);

      if (molTest->reapplyMolBlockWedging) {
        Chirality::reapplyMolBlockWedging(*localVars.mol);
      }

      TEST_ASSERT(localVars.mol != nullptr);

      TEST_ASSERT(localVars.mol->getNumAtoms() == molTest->atomCount)
      TEST_ASSERT(localVars.mol->getNumBonds() == molTest->bondCount)

      {
        std::string expectedMrvName =
            fName + (molTest->sanitizeFlag ? "" : ".nosan") + ".expected.sdf";
        std::string outMolStr = "";
        try {
          outMolStr = MolToMolBlock(*localVars.mol, true, 0, true, true);
        } catch (const RDKit::KekulizeException &e) {
          outMolStr = "";
        }
        if (outMolStr == "") {
          outMolStr = MolToMolBlock(*localVars.mol, true, 0, false,
                                    true);  // try without kekule'ing
        }

        generateNewExpectedFilesIfSoSpecified(
            fName + (molTest->sanitizeFlag ? "" : ".nosan") + ".NEW.sdf",
            outMolStr);

        TEST_ASSERT(GetExpectedValue(expectedMrvName) == outMolStr);
      }

      {
        std::string expectedMrvName =
            fName + (molTest->sanitizeFlag ? "" : ".nosan") + ".expected.mrv";

        std::string outMolStr = "";
        try {
          outMolStr = MolToMrvBlock(*localVars.mol, true, -1, true, false);
        } catch (const RDKit::KekulizeException &e) {
          outMolStr = "";
        }
        if (outMolStr == "") {
          outMolStr = MolToMrvBlock(*localVars.mol, true, -1, false,
                                    false);  // try without kekule'ing
        }

        generateNewExpectedFilesIfSoSpecified(
            fName + (molTest->sanitizeFlag ? "" : ".nosan") + ".NEW.mrv",
            outMolStr);

        TEST_ASSERT(GetExpectedValue(expectedMrvName) == outMolStr);

        BOOST_LOG(rdInfoLog) << "done" << std::endl;
      }
    } catch (const std::exception &e) {
      if (molTest->expectedResult != false) {
        throw;
      }
      return;
    }
    TEST_ASSERT(molTest->expectedResult == true);

    return;
  }

  void testMarvinRxn(const RxnTest *rxnTest) {
    BOOST_LOG(rdInfoLog) << "testing marvin parsing" << std::endl;

    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/MarvinParse/test_data/" + rxnTest->fileName;

    class LocalVars  // protect against mem leak on error
    {
     public:
      ChemicalReaction *rxn = nullptr;

      LocalVars(){};

      ~LocalVars() { delete (ChemicalReaction *)rxn; }
    } localVars;

    try {
      if (!MrvFileIsReaction(fName)) {
        TEST_ASSERT(rxnTest->expectedResult == false);
        return;
      }

      localVars.rxn = GetReaction(rxnTest);

      // check for errors

      unsigned int nWarn = 0, nError = 0;

      TEST_ASSERT(localVars.rxn != nullptr);

      TEST_ASSERT(localVars.rxn->getNumReactantTemplates() ==
                  rxnTest->reactantCount);
      TEST_ASSERT(localVars.rxn->getNumProductTemplates() ==
                  rxnTest->productCount);
      TEST_ASSERT(localVars.rxn->getNumAgentTemplates() == rxnTest->agentCount);
      localVars.rxn->initReactantMatchers(true);

      if (localVars.rxn->getNumReactantTemplates() > 0 &&
          localVars.rxn->getNumProductTemplates() > 0) {
        TEST_ASSERT(localVars.rxn->validate(nWarn, nError, true));
      } else {
        nWarn = 0;
        nError = 0;
      }

      TEST_ASSERT(nWarn == rxnTest->warnings);
      TEST_ASSERT(nError == rxnTest->errors);

      // make sure the Rxn is kekule'ed

      for (auto mol : localVars.rxn->getReactants()) {
        auto rwMol = (RWMol *)mol.get();
        if (rwMol->needsUpdatePropertyCache()) {
          rwMol->updatePropertyCache(false);
        }
        MolOps::Kekulize(*rwMol);
      }
      for (auto mol : localVars.rxn->getAgents()) {
        auto rwMol = (RWMol *)mol.get();
        if (rwMol->needsUpdatePropertyCache()) {
          rwMol->updatePropertyCache(false);
        }
        MolOps::Kekulize(*rwMol);
      }
      for (auto mol : localVars.rxn->getProducts()) {
        auto rwMol = (RWMol *)mol.get();
        if (rwMol->needsUpdatePropertyCache()) {
          rwMol->updatePropertyCache(false);
        }
        MolOps::Kekulize(*rwMol);
      }

      {
        std::string outMolStr =
            ChemicalReactionToRxnBlock(*localVars.rxn, false, true);

        std::string expectedRxnName = fName + ".expected.rxn";

        generateNewExpectedFilesIfSoSpecified(fName + ".NEW.rxn", outMolStr);

        TEST_ASSERT(GetExpectedValue(expectedRxnName) == outMolStr);
      }

      {
        std::string outMolStr =
            ChemicalReactionToMrvBlock(*localVars.rxn, false);

        std::string expectedRxnName = fName + ".expected.mrv";

        generateNewExpectedFilesIfSoSpecified(fName + ".NEW.mrv", outMolStr);

        TEST_ASSERT(GetExpectedValue(expectedRxnName) == outMolStr);
      }
      BOOST_LOG(rdInfoLog) << "done" << std::endl;
      BOOST_LOG(rdInfoLog) << "done" << std::endl;
    } catch (const std::exception &e) {
      if (rxnTest->expectedResult != false) {
        throw;
      }
      return;
    }

    TEST_ASSERT(rxnTest->expectedResult == true);

    return;
  }

  void testMolFiles(const MolTest *molFileTest) {
    BOOST_LOG(rdInfoLog) << "testing marvin writing" << std::endl;

    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase + "/Code/GraphMol/MarvinParse/test_data/" +
                        molFileTest->fileName;

    class LocalVars  // protext against mem leak on error
    {
     public:
      RWMol *mol = nullptr;

      LocalVars(){};

      ~LocalVars() { delete mol; }
    } localVars;

    try {
      localVars.mol =
          MolFileToMol(fName, molFileTest->sanitizeFlag, false, false);
      if (molFileTest->reapplyMolBlockWedging) {
        reapplyMolBlockWedging(*localVars.mol);
      }

      TEST_ASSERT(localVars.mol != nullptr);
      TEST_ASSERT(localVars.mol->getNumAtoms() == molFileTest->atomCount)
      TEST_ASSERT(localVars.mol->getNumBonds() == molFileTest->bondCount)

      {
        std::string expectedMrvName =
            fName + (molFileTest->sanitizeFlag ? "" : ".nosan") +
            ".expected.sdf";
        std::string outMolStr = "";
        try {
          outMolStr = MolToMolBlock(*localVars.mol, true, 0, true, true);
        } catch (const RDKit::KekulizeException &e) {
          outMolStr = "";
        }
        if (outMolStr == "") {
          outMolStr = MolToMolBlock(*localVars.mol, true, 0, false,
                                    true);  // try without kekule'ing
        }

        generateNewExpectedFilesIfSoSpecified(
            fName + (molFileTest->sanitizeFlag ? "" : ".nosan") + ".NEW.sdf",
            outMolStr);
        TEST_ASSERT(GetExpectedValue(expectedMrvName) == outMolStr);
      }

      {
        std::string expectedMrvName =
            fName + (molFileTest->sanitizeFlag ? "" : ".nosan") +
            ".expected.mrv";

        std::string outMolStr = "";
        try {
          outMolStr = MolToMrvBlock(*localVars.mol, true, -1, true, false);
        } catch (const RDKit::KekulizeException &e) {
          outMolStr = "";
        }
        if (outMolStr == "") {
          // try without kekule'ing
          outMolStr = MolToMrvBlock(*localVars.mol, true, -1, false, false);
        }

        generateNewExpectedFilesIfSoSpecified(
            fName + (molFileTest->sanitizeFlag ? "" : ".nosan") + ".NEW.mrv",
            outMolStr);
        TEST_ASSERT(GetExpectedValue(expectedMrvName) == outMolStr);
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
    RDKit::Chirality::setUseLegacyStereoPerception(false);
    printf("Using new chirality perception\n");

    // the molecule tests - starting with molfiles/sdf
    if (testToRun == "" || testToRun == "sdfTests") {
      std::list<MolTest> sdfTests{
          MolTest("NewChiralTest.sdf", true, 13, 14, true,
                  false),  // wedges NOT reapplied
          MolTest("NewChiralTest.sdf", true, 13, 14, false,
                  false),  // not sanitized, wedges NOT reapplied
          MolTest("NewChiralTestNoChiral.sdf", true, 13, 14, true,
                  false),  // wedges NOT reapplied
          MolTest("NewChiralTestNoChiral.sdf", true, 13, 14, false,
                  false),  // not sanitized, wedges NOT reapplied
          MolTest("NewChiralTestAllChiral.sdf", true, 14, 15, true,
                  false),  // wedges NOT reapplied
          MolTest("NewChiralTestAllChiral.sdf", true, 14, 15, false,
                  false),  // not sanitized, wedges NOT reapplied
          MolTest("ProblemShort.mol", true, 993, 992),
          MolTest("lostStereoAnd.sdf", true, 6, 5),
          MolTest("DoubleBondChain.sdf", true, 22, 22),
          MolTest("UnitsError.sdf", true, 17, 18),
          MolTest("StarAtom.sdf", true, 17, 16),
          MolTest("nonProprietary.mol", true, 17, 16),
          MolTest("ChiralTest.sdf", true, 8, 7),
          MolTest("TestBond1.mol", true, 10, 10),
          MolTest("Sgroup_MUL_ParentInMiddle.sdf", true, 17, 16)};

      for (auto sdfTest : sdfTests) {
        BOOST_LOG(rdInfoLog) << "Test: " << sdfTest.fileName << std::endl;

        printf("Test\n\n %s\n\n", sdfTest.fileName.c_str());
        testMolFiles(&sdfTest);
      }
    }

    if (testToRun == "" || testToRun == "molFileTests") {
      std::list<MolTest> molFileTests{
          MolTest("Cubane.mrv", true, 16, 20),
          MolTest("NewChiralTest.mrv", true, 13, 14, true,
                  false),  // wedges NOT reapplied
          MolTest("NewChiralTest.mrv", true, 13, 14, false,
                  false),  // not sanitized, wedges NOT reapplied
          MolTest("NewChiralTestNoChiral.mrv", true, 13, 14, true,
                  false),  // wedges NOT reapplied
          MolTest("NewChiralTestNoChiral.mrv", true, 13, 14, false,
                  false),  // not sanitized, wedges NOT reapplied
          MolTest("NewChiralTestAllChiral.mrv", true, 14, 15, true,
                  false),  // wedges NOT reapplied
          MolTest("NewChiralTestAllChiral.mrv", true, 14, 15, false,
                  false),  // not sanitized, wedges NOT reapplied
          MolTest("lostStereoAnd.mrv", true, 6, 5),
          MolTest("DoubleBondChain.mrv", true, 22, 22),
          MolTest("WigglyAndCrossed.mrv", true, 8, 7),
          MolTest("BondTypes.mrv", true, 26, 25),
          MolTest("EmbeddedSGroupSUP_MUL.mrv", true, 17, 17),
          MolTest("EmbeddedSgroupCOP_SUP.mrv", true, 10, 10),
          MolTest("EmbeddedSgroupDAT_SUP.mrv", true, 10, 10),
          MolTest("EmbeddedSgroupMULTICENTER_SUP.mrv", true, 11, 10),
          MolTest("EmbeddedSgroupMUL_MUL.mrv", true, 141, 140),
          MolTest("EmbeddedSgroupMUL_MUL2.mrv", true, 23, 22),
          MolTest("EmbeddedSgroupMUL_SUP.mrv", true, 129, 128),
          MolTest("EmbeddedSgroupSRU_SUP.mrv", true, 10, 10),
          MolTest("EmbeddedSgroupSUPEXP_SUP.mrv", true, 10, 10),
          MolTest("EmbeddedSgroupSUPEXP_SUP2.mrv", true, 10, 10),
          MolTest("EmbeddedSgroupSUP_MULTICENTER.mrv", true, 10, 8),
          MolTest("EmbeddedSgroupSUP_SUP.mrv", true, 12, 11,
                  false),  // not sanitized
          MolTest("EmbeddedSgroupSUP_SUP2.mrv", true, 12, 12),
          MolTest("RgroupBad.mrv", true, 9, 9),
          MolTest("valenceLessThanDrawn.mrv", true, 14, 14),
          MolTest("data_sgroup_no_fieldname.mrv", true, 4, 3),
          MolTest("data_sgroup_empty_field_data.mrv", true, 2, 1),
          MolTest("radical_value.mrv", true, 3, 2),
          MolTest("emptyOneLineAtomList.mrv", true, 0, 0),
          MolTest("mrvValence_value.mrv", true, 3, 2),
          MolTest("ChiralTest2.mrv", true, 46, 47),
          MolTest("ChiralTest.mrv", true, 8, 7),
          MolTest("SnCl2.mrv", true, 3, 2),
          MolTest("SnH2Cl2.mrv", true, 3, 2),
          MolTest("marvin01.mrv", true, 11, 11),
          MolTest("marvin02.mrv", true, 9, 9),
          MolTest("marvin07.mrv", true, 12, 11),
          MolTest("marvin10.mrv", true, 10, 10),
          MolTest("marvin06.mrv", true, 11, 11),
          MolTest("marvin12.mrv", true, 31, 33),
          MolTest("EmptyMol.mrv", true, 0, 0),
          MolTest("Sparse.mrv", true, 0, 0),
          MolTest("Sparse2.mrv", true, 0, 0),
          MolTest("Sparse3.mrv", true, 0, 0),
          MolTest("MarvinNoCoords.mrv", true, 6, 6),
          MolTest("aspirin.mrv", true, 13, 13),
          MolTest("MarvinStereoGroupsZeros.mrv", true, 8, 8),
          MolTest("MarvinStereoGroupsAbs.mrv", true, 8, 8),
          MolTest("triphenylphosphine.mrv", true, 19, 21),
          MolTest("MarvinOldSuperGroupTest.mrv", true, 89, 93),
          MolTest("RadicalTests.mrv", true, 8, 7),
          MolTest("AnyBond.mrv", true, 4, 3),
          MolTest("cisBenzene.mrv", true, 6, 6),
          MolTest("DativeBond.mrv", true, 6, 5),
          MolTest("MultipleSgroup.mrv", true, 123, 122),
          MolTest("SgroupExpanded.mrv", true, 5, 4),
          MolTest("SgroupMultAttach.mrv", true, 44, 45),
          MolTest("MarvinMissingX2.mrv", true, 12, 11),
          MolTest("MarvinMissingY2.mrv", true, 12, 11),
          MolTest("DataSgroup.mrv", true, 7, 6),
          MolTest("MulticenterSgroup.mrv", true, 17, 16,
                  false),  // not sanitized),
          MolTest("GenericSgroup.mrv", true, 13, 13),
          MolTest("MonomerSgroup.mrv", true, 4, 3),
          MolTest("modification_sgroup.mrv", true, 54, 40),
          MolTest("copolymer_sgroup.mrv", true, 19, 18),
          MolTest("MultipleSgroupParentInMiddleOfAtomBlock.mrv", true, 23, 22),
          MolTest("EmbeddedSgroups.mrv", false, 14, 14),
          MolTest("marvin03.mrv", false, 31, 33),
          MolTest("MarvinBadMissingMolID.mrv", false, 12, 11),
          MolTest("MarvinBadMissingAtomID.mrv", false, 12, 11),
          MolTest("MarvinBadX2.mrv", false, 12, 11),
          MolTest("MarvinBadY2.mrv", false, 12, 11),
          MolTest("MarvinBadStereoGroupsAbs.mrv", false, 8, 8),
          MolTest("MarvinBadElementType.mrv", false, 12, 11),
          MolTest("MarvinBadMissingBondID.mrv", false, 12, 11),
          MolTest("MarvinBadMissingBondAtomRefs", false, 12, 11),
          MolTest("MarvinBadMissingBondOrder.mrv", false, 12, 11),
          MolTest("MarvinBadMissingSruMolID.mrv", false, 12, 11),
          MolTest("MarvinBadMissingSruID.mrv", false, 12, 11),
          MolTest("MarvinBadMissingSruRole.mrv", false, 12, 11),
          MolTest("MarvinBadMissingSruAtomRef.mrv", false, 12, 11),
          MolTest("MarvinBadMissingSruTitle.mrv", false, 12, 11),
          MolTest("MarvinBadSruAtomRef.mrv", false, 12, 11),
          MolTest("MarvinBadSruID.mrv", false, 12, 11),
          MolTest("MarvinBadSruRole.mrv", false, 12, 11),
          MolTest("MarvinBadSruAtomRef.mrv", false, 12, 11),
          MolTest("MarvinBadSruAtomRef.mrv", false, 12, 11),
          MolTest("MarvinBadSruConnect.mrv", false, 12, 11),
          MolTest("MarvinBadSupAttachAtom.mrv", false, 9, 9),
          MolTest("MarvinBadSupAttachBond.mrv", false, 9, 9),
          MolTest("MarvinBadSupAttachOrder.mrv", false, 9, 9),
          MolTest("MarvinBadSupAttachAtom.mrv", false, 9, 9),
          MolTest("MarvinBadSupAttachAtom.mrv", false, 9, 9),
          MolTest("MarvinBadSupAttachAtom.mrv", false, 9, 9),
          MolTest("MarvinBadSupMissingAttachBond.mrv", false, 9, 9),
          MolTest("MarvinBadSupMissingAttachOrder.mrv", false, 9, 9),
          MolTest("marvin03.mrv", false, 1, 1)};  // should fail

      for (auto molFileTest : molFileTests) {
        BOOST_LOG(rdInfoLog) << "Test: " << molFileTest.fileName << std::endl;

        printf("Test\n\n %s\n\n", molFileTest.fileName.c_str());
        testMarvinMol(&molFileTest);
      }
    }
    // now the reactions

    if (testToRun == "" || testToRun == "rxnFileTests") {
      std::list<RxnTest> rxnFileTests{

          RxnTest("AlexRxn.mrv", true, 1, 0, 1, 2, 0),
          RxnTest("BadReactionSign.mrv", true, 2, 0, 1, 3, 0),
          RxnTest("bondArray_node.mrv", true, 2, 4, 1, 3, 0),
          RxnTest("marvin03.mrv", true, 1, 1, 1, 2, 0),
          RxnTest("marvin04.mrv", true, 2, 1, 2, 4, 0),
          RxnTest("marvin08.mrv", true, 2, 3, 2, 4, 0),
          RxnTest("marvin09.mrv", true, 2, 3, 2, 4, 0),
          RxnTest("marvin11.mrv", true, 2, 0, 1, 0, 0),
          RxnTest("marvin05.mrv", true, 2, 1, 1, 3, 0),
          RxnTest("EmptyRxn.mrv", true, 0, 0, 0, 0, 0),
          RxnTest("RxnNoCoords.mrv", true, 2, 0, 1, 3, 0),
          RxnTest("mrvValenceZero.mrv", true, 3, 0, 1, 4, 0),
          RxnTest("condition_coordinates_mpoint.mrv", true, 1, 0, 1, 0, 0),
          RxnTest("marvin01.mrv", false, 2, 1, 1, 3, 0),
          RxnTest("aspirineSynthesisWithAttributes.mrv", true, 2, 0, 1, 3, 0),
          RxnTest("marvin01.mrv", false, 2, 1, 1, 3, 0)};  // should fail

      for (auto rxnFileTest : rxnFileTests) {
        printf("Test\n\n %s\n\n", rxnFileTest.fileName.c_str());
        testMarvinRxn(&rxnFileTest);
      }
    }

    // now smiles tests

    if (testToRun == "" || testToRun == "smiTests") {
      std::list<SmilesTest> smiTests{
          SmilesTest("NewChiralTest", R"(C[C@H]1CC[C@@]2(CC[C@H](Cl)CC2)CC1)",
                     true, 13, 14, true),
          SmilesTest("NewChiralTest", R"(C[C@H]1CC[C@@]2(CC[C@H](Cl)CC2)CC1)",
                     true, 13, 14, false),
          SmilesTest("NewChiralTestNoChiral",
                     R"(C[C@H]1CCC2(CC[C@H](Cl)CC2)CC1)", true, 13, 14, true),
          SmilesTest("NewChiralTestNoChiral",
                     R"(C[C@H]1CCC2(CC[C@H](Cl)CC2)CC1)", true, 13, 14, false),
          SmilesTest("NewChiralTest2", R"(C[C@H]1CCC2(CC1)CC[C@H](C)C(C)C2)",
                     true, 14, 15, true),
          SmilesTest("NewChiralTest2", R"(C[C@H]1CCC2(CC1)CC[C@H](C)C(C)C2)",
                     true, 14, 15, false),
          SmilesTest("NewChiralTest2AllChiral",
                     R"(C[C@H]1CC[C@@]2(CC1)CC[C@H](C)C(C)C2)", true, 14, 15,
                     true),
          SmilesTest("NewChiralTest2AllChiral",
                     R"(C[C@H]1CC[C@@]2(CC1)CC[C@H](C)C(C)C2)", true, 14, 15,
                     false),

          SmilesTest("DoubleBondChain",
                     R"(CC1=C(\C=C\C(C)=C\C=C\C(C)=C/C(O)=O)C(C)(C)CCC1)", true,
                     22, 22),
          // SmilesTest(
          //     "Macrocycle2",
          //     R"(CC1OC(=O)CC(O)CC(O)CC(O)CCC(O)C(O)CC2(O)CC(O)C(C(CC(O[C@@H]3O[C@H](C)[C@@H](O)[C@H](N)[C@@H]3O)\C=C\C=C\C=C\C=C\CC\C=C\C=C\C(C)C(O)C1C)O2)C(O)=O
          //     |t:42,44,46,48,52,54|)",
          //     true, 65, 67),
          SmilesTest(
              "Na_Mg_Al_OH",
              "[OH-].[OH-].[OH-].[O--].[Na+].[Mg++].[Al+3].[Si].OC([O-])=O",
              true, 12, 3),
          SmilesTest("Pb", "[Pb]", true, 1, 0),
          SmilesTest("O_Mg_Si", "[O].[Mg].[Si]", true, 3, 0),
          SmilesTest("SquiggleBond", "CN1N=C(SC1=NC(C)=O)S(N)(=O)=O |c:2|",
                     true, 14, 14),
          SmilesTest(
              "BigMacrocycle",
              "C[C@@H]1CCCCCCCCC(=O)OCCN[C@H](C)CCCCCCCCC(=O)OCCN[C@H](C)CCCCCCCCC(=O)OCCN1",
              true, 48, 48),
          SmilesTest("Smiles1", "N[C@@H]([O-])c1cc[13c]cc1", true, 9, 9)};

      for (auto smiTest : smiTests) {
        printf("Test\n\n %s\n\n", smiTest.name.c_str());
        testSmilesToMarvin(&smiTest);
      }
    }
  }
};

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  MrvTests mrvTests;

  if (argc > 1) {
    mrvTests.testToRun = argv[1];
  }

  if (argc > 2 && std::string(argv[2]) == "generate") {
    mrvTests.generateExpectedFiles = true;
  }

  RDLog::InitLogs();
  BOOST_LOG(rdInfoLog) << " ---- Running with POSIX locale ----- " << std::endl;

  mrvTests.RunTests();  // run with C locale

  return 0;
}
