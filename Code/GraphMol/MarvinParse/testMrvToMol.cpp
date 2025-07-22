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
#include <GraphMol/test_fixtures.h>
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
using namespace v2::FileParsers;

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
          reapplyMolBlockWedging(reapplyMolBlockWedgingInit) {};
  };

  class ScsrMolTest {
   public:
    unsigned int atomCount;
    unsigned int bondCount;
    unsigned int sGroupCount;
    std::string fileName;
    bool expectedResult;

    ScsrMolTest(std::string fileNameInit, bool expectedResultInit,
                int atomCountInit, int bondCountInit, int sGroupCountInit)
        : atomCount(atomCountInit),
          bondCount(bondCountInit),
          sGroupCount(sGroupCountInit),
          fileName(fileNameInit),
          expectedResult(expectedResultInit) {};
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
          errors(errorInit) {};
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
          bondCount(bondCountInit) {};

    SmilesTest(std::string nameInit, std::string smilesInit,
               bool expectedResultInit, int atomCountInit, int bondCountInit)
        : name(nameInit),
          smiles(smilesInit),
          expectedResult(expectedResultInit),
          sanitizeFlag(true),
          atomCount(atomCountInit),
          bondCount(bondCountInit) {};
  };

  RWMol *GetMolv1(const MolTest *molTest) {
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

  ChemicalReaction *GetReactionv1(const RxnTest *rxnTest) {
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

  std::unique_ptr<RWMol> GetMolv2(const MolTest *molTest) {
    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/MarvinParse/test_data/" + molTest->fileName;

    v2::MarvinParser::MrvParserParams params;
    params.removeHs = false;
    params.sanitize = molTest->sanitizeFlag;
    try {
      return v2::MarvinParser::MolFromMrvFile(fName, params);
    } catch (const std::exception &e) {
      std::cerr << e.what() << '\n';
      throw BadFileException("Could not parse the MRV block");
    }
  }

  std::unique_ptr<ChemicalReaction> GetReactionv2(const RxnTest *rxnTest) {
    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/MarvinParse/test_data/" + rxnTest->fileName;

    v2::MarvinParser::MrvParserParams params;
    params.removeHs = false;
    params.sanitize = true;
    try {
      return ReactionFromMrvFile(fName, params);
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

    try {
      SmilesParserParams smilesParserParams;
      smilesParserParams.sanitize = smilesTest->sanitizeFlag;
      smilesParserParams.allowCXSMILES = true;
      smilesParserParams.removeHs = false;  // do not remove Hs

      std::unique_ptr<RWMol> smilesMol(
          SmilesToMol(smilesTest->smiles, smilesParserParams));

      Chirality::reapplyMolBlockWedging(*smilesMol);

      TEST_ASSERT(smilesMol->getNumAtoms() == smilesTest->atomCount);
      TEST_ASSERT(smilesMol->getNumBonds() == smilesTest->bondCount);

      // test round trip back to smiles
      {
        std::string expectedMrvName =
            fName + (smilesTest->sanitizeFlag ? "" : ".nosan") +
            ".expected.smi";

        SmilesWriteParams ps;
        ps.canonical = false;

        std::string smilesOut = MolToSmiles(*smilesMol, ps);

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
        outMolStr = MolToMolBlock(*smilesMol, true, 0, true, true);

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
          outMolStr = MolToMrvBlock(*smilesMol, true, -1, true, false);
        } catch (const RDKit::KekulizeException &) {
          outMolStr = "";

          if (outMolStr == "") {
            outMolStr = MolToMrvBlock(*smilesMol, true, -1, false,
                                      false);  // try without kekule'ing
          }

          generateNewExpectedFilesIfSoSpecified(
              fName + (smilesTest->sanitizeFlag ? "" : ".nosan") + ".NEW.mrv",
              outMolStr);

          TEST_ASSERT(GetExpectedValue(expectedMrvName) == outMolStr);
        }
        BOOST_LOG(rdInfoLog) << "done" << std::endl;
      }
    } catch (const std::exception &) {
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

    try {
      if (MrvFileIsReaction(fName)) {
        TEST_ASSERT(molTest->expectedResult == false);
        return;
      }

      {
        // make sure the v1 API still works, this calls the v2 API, so we don't
        // need to do serious correctness testing here, we'll let the v2 tests
        // below handle that
        std::unique_ptr<RWMol> mol2(GetMolv1(molTest));
        TEST_ASSERT(mol2 != nullptr);
        TEST_ASSERT(mol2->getNumAtoms() == molTest->atomCount)
        TEST_ASSERT(mol2->getNumBonds() == molTest->bondCount)
      }

      auto mol = GetMolv2(molTest);
      TEST_ASSERT(mol != nullptr);
      RDKit::Chirality::removeNonExplicit3DChirality(*mol);

      if (molTest->reapplyMolBlockWedging) {
        Chirality::reapplyMolBlockWedging(*mol);
      }

      TEST_ASSERT(mol->getNumAtoms() == molTest->atomCount);
      TEST_ASSERT(mol->getNumBonds() == molTest->bondCount);

      {
        std::string expectedMrvName =
            fName + (molTest->sanitizeFlag ? "" : ".nosan") +
            (molTest->reapplyMolBlockWedging ? "" : ".noreapply") +
            ".expected.sdf";
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

        generateNewExpectedFilesIfSoSpecified(
            fName + (molTest->sanitizeFlag ? "" : ".nosan") +
                (molTest->reapplyMolBlockWedging ? "" : ".noreapply") +
                ".NEW.sdf",
            outMolStr);

        TEST_ASSERT(GetExpectedValue(expectedMrvName) == outMolStr);
      }

      {
        std::string expectedMrvName =
            fName + (molTest->sanitizeFlag ? "" : ".nosan") +
            (molTest->reapplyMolBlockWedging ? "" : ".noreapply") +
            ".expected.mrv";

        std::string outMolStr = "";
        try {
          outMolStr = MolToMrvBlock(*mol, true, -1, true, false);
        } catch (const RDKit::KekulizeException &) {
          outMolStr = "";
        }
        if (outMolStr == "") {
          outMolStr = MolToMrvBlock(*mol, true, -1, false,
                                    false);  // try without kekule'ing
        }
        generateNewExpectedFilesIfSoSpecified(
            fName + (molTest->sanitizeFlag ? "" : ".nosan") +
                (molTest->reapplyMolBlockWedging ? "" : ".noreapply") +
                ".NEW.mrv",
            outMolStr);

        TEST_ASSERT(GetExpectedValue(expectedMrvName) == outMolStr);

        BOOST_LOG(rdInfoLog) << "done" << std::endl;
      }
    } catch (const std::exception &) {
      if (molTest->expectedResult != false) {
        throw;
      }
      return;
    }
    TEST_ASSERT(molTest->expectedResult == true);

    return;
  }

  void testMarvinFromScsr(const ScsrMolTest *scsrMolTest) {
    BOOST_LOG(rdInfoLog) << "testing marvin generate from an SCSR mol file"
                         << std::endl;

    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase + "/Code/GraphMol/MarvinParse/test_data/" +
                        scsrMolTest->fileName;

    try {
      RDKit::v2::FileParsers::MolFileParserParams pp;
      pp.sanitize = false;
      pp.removeHs = false;
      pp.strictParsing = true;

      RDKit::v2::FileParsers::MolFromSCSRParams molFromSCSRParams;
      molFromSCSRParams.includeLeavingGroups = true;
      molFromSCSRParams.scsrBaseHbondOptions = SCSRBaseHbondOptions::Auto;

      std::unique_ptr<RDKit::RWMol> mol;
      mol = MolFromSCSRFile(fName, pp, molFromSCSRParams);

      TEST_ASSERT(mol != nullptr);
      TEST_ASSERT(mol->getNumAtoms() == scsrMolTest->atomCount);
      TEST_ASSERT(mol->getNumBonds() == scsrMolTest->bondCount);
      TEST_ASSERT(getSubstanceGroups(*mol).size() == scsrMolTest->sGroupCount);

      std::string outMolStr = "";
      outMolStr = MolToMrvBlock(*mol, true, -1, true, false);

      RDKit::v2::MarvinParser::MrvParserParams mpp;
      mol = RDKit::v2::MarvinParser::MolFromMrvBlock(outMolStr, mpp);

      TEST_ASSERT(mol != nullptr);
      TEST_ASSERT(mol->getNumAtoms() == scsrMolTest->atomCount);
      TEST_ASSERT(mol->getNumBonds() == scsrMolTest->bondCount);
      TEST_ASSERT(getSubstanceGroups(*mol).size() == scsrMolTest->sGroupCount);

      BOOST_LOG(rdInfoLog) << "done" << std::endl;
    } catch (const std::exception &) {
      if (scsrMolTest->expectedResult != false) {
        throw;
      }
      return;
    }
    TEST_ASSERT(scsrMolTest->expectedResult == true);

    return;
  }

  void testMarvinRxn(const RxnTest *rxnTest) {
    BOOST_LOG(rdInfoLog) << "testing marvin parsing" << std::endl;

    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/MarvinParse/test_data/" + rxnTest->fileName;

    try {
      if (!MrvFileIsReaction(fName)) {
        TEST_ASSERT(rxnTest->expectedResult == false);
        return;
      }

      // reaction test
      std::unique_ptr<ChemicalReaction> rxn(GetReactionv1(rxnTest));

      // check for errors

      unsigned int nWarn = 0, nError = 0;

      TEST_ASSERT(rxn != nullptr);

      TEST_ASSERT(rxn->getNumReactantTemplates() == rxnTest->reactantCount);
      TEST_ASSERT(rxn->getNumProductTemplates() == rxnTest->productCount);
      TEST_ASSERT(rxn->getNumAgentTemplates() == rxnTest->agentCount);
      rxn->initReactantMatchers(true);

      if (rxn->getNumReactantTemplates() > 0 &&
          rxn->getNumProductTemplates() > 0) {
        TEST_ASSERT(rxn->validate(nWarn, nError, true));
      } else {
        nWarn = 0;
        nError = 0;
      }

      TEST_ASSERT(nWarn == rxnTest->warnings);
      TEST_ASSERT(nError == rxnTest->errors);

      // make sure the Rxn is kekule'ed

      for (auto mol : rxn->getReactants()) {
        auto rwMol = (RWMol *)mol.get();
        if (rwMol->needsUpdatePropertyCache()) {
          rwMol->updatePropertyCache(false);
        }
        MolOps::Kekulize(*rwMol);
      }
      for (auto mol : rxn->getAgents()) {
        auto rwMol = (RWMol *)mol.get();
        if (rwMol->needsUpdatePropertyCache()) {
          rwMol->updatePropertyCache(false);
        }
        MolOps::Kekulize(*rwMol);
      }
      for (auto mol : rxn->getProducts()) {
        auto rwMol = (RWMol *)mol.get();
        if (rwMol->needsUpdatePropertyCache()) {
          rwMol->updatePropertyCache(false);
        }
        MolOps::Kekulize(*rwMol);
      }

      {
        std::string outMolStr = ChemicalReactionToRxnBlock(*rxn, false, true);

        generateNewExpectedFilesIfSoSpecified(fName + ".NEW.rxn", outMolStr);

        std::string expectedRxnName = fName + ".expected.rxn";

        TEST_ASSERT(GetExpectedValue(expectedRxnName) == outMolStr);
      }

      {
        std::string outMolStr = ChemicalReactionToMrvBlock(*rxn, false);

        generateNewExpectedFilesIfSoSpecified(fName + ".NEW.mrv", outMolStr);

        std::string expectedRxnName = fName + ".expected.mrv";

        TEST_ASSERT(GetExpectedValue(expectedRxnName) == outMolStr);
      }

      BOOST_LOG(rdInfoLog) << "done" << std::endl;
    } catch (const std::exception &) {
      if (rxnTest->expectedResult != false) {
        throw;
      }
      return;
    }

    TEST_ASSERT(rxnTest->expectedResult == true);

    return;
  }

  void testMarvinRxnMols(const RxnTest *rxnTest) {
    BOOST_LOG(rdInfoLog) << "testing marvin parsing" << std::endl;

    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/MarvinParse/test_data/" + rxnTest->fileName;

    try {
      std::unique_ptr<ChemicalReaction> rxn(
          MrvFileToChemicalReaction(fName, false, false));

      std::unique_ptr<ROMol> oneMol(ChemicalReactionToRxnMol(*rxn));

      auto rwMol = (RWMol *)oneMol.get();
      if (rwMol->needsUpdatePropertyCache()) {
        rwMol->updatePropertyCache(false);
      }
      MolOps::Kekulize(*rwMol);
      RDKit::Chirality::reapplyMolBlockWedging(*rwMol);

      {
        std::string outMolStr = MolToMrvBlock(*rwMol, false, -1, false, false);

        std::string expectedRxnName = fName + ".expected.mrv";

        generateNewExpectedFilesIfSoSpecified(fName + ".NEW.mrv", outMolStr);

        TEST_ASSERT(GetExpectedValue(expectedRxnName) == outMolStr);
      }
      BOOST_LOG(rdInfoLog) << "done" << std::endl;
      BOOST_LOG(rdInfoLog) << "done" << std::endl;
    } catch (const std::exception &) {
      if (rxnTest->expectedResult != false) {
        throw;
      }
      return;
    }

    TEST_ASSERT(rxnTest->expectedResult == true);

    return;
  }

  void testMarvin3dChiral(const MolTest *molTest) {
    BOOST_LOG(rdInfoLog) << "testing marvin parsing for chirality from 3d"
                         << std::endl;

    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/MarvinParse/test_data/" + molTest->fileName;

    try {
      if (MrvFileIsReaction(fName)) {
        TEST_ASSERT(molTest->expectedResult == false);
        return;
      }

      auto mol = GetMolv2(molTest);
      TEST_ASSERT(mol != nullptr);
      RDKit::Chirality::removeNonExplicit3DChirality(*mol);

      Chirality::reapplyMolBlockWedging(*mol);

      TEST_ASSERT(mol->getNumAtoms() == molTest->atomCount)
      TEST_ASSERT(mol->getNumBonds() == molTest->bondCount)

      {
        std::string expectedMrvName = fName + ".expected.sdf";
        std::string outMolStr = "";
        try {
          outMolStr = MolToMolBlock(*mol, true, 0, true, true);
        } catch (const RDKit::KekulizeException &) {
          outMolStr = "";
        } catch (...) {
          throw;  // re-trhow the error if not a kekule error
        }
        if (outMolStr == "") {
          outMolStr = MolToMolBlock(*mol, true, 0, false,
                                    true);  // try without kekule'ing
        }

        generateNewExpectedFilesIfSoSpecified(fName + ".NEW.sdf", outMolStr);

        TEST_ASSERT(GetExpectedValue(expectedMrvName) == outMolStr);
      }

      {
        std::string expectedMrvName = fName + ".expected.mrv";

        std::string outMolStr = "";
        try {
          outMolStr = MolToMrvBlock(*mol, true, -1, true, false);
          RDKit::Chirality::removeNonExplicit3DChirality(*mol);

        } catch (const RDKit::KekulizeException &) {
          outMolStr = "";
        } catch (...) {
          throw;  // re-throw the error if not a kekule error
        }
        if (outMolStr == "") {
          outMolStr = MolToMrvBlock(*mol, true, -1, false,
                                    false);  // try without kekule'ing
        }
        generateNewExpectedFilesIfSoSpecified(fName + ".NEW.mrv", outMolStr);

        TEST_ASSERT(GetExpectedValue(expectedMrvName) == outMolStr);

        BOOST_LOG(rdInfoLog) << "done" << std::endl;
      }
    } catch (const std::exception &) {
      if (molTest->expectedResult != false) {
        throw;
      }
      return;
    }
    TEST_ASSERT(molTest->expectedResult == true);

    return;
  }

  void testMarvinAtrop(const MolTest *molTest) {
    BOOST_LOG(rdInfoLog) << "testing marvin atropisomers" << std::endl;

    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/MarvinParse/test_data/" + molTest->fileName;
    v2::MarvinParser::MrvParserParams params;
    params.removeHs = false;
    params.sanitize = molTest->sanitizeFlag;
    try {
      auto mol = MolFromMrvFile(fName, params);
      // mol  test
      TEST_ASSERT(mol != nullptr);
      RDKit::Chirality::removeNonExplicit3DChirality(*mol);

      TEST_ASSERT(mol->getNumAtoms() == molTest->atomCount)
      TEST_ASSERT(mol->getNumBonds() == molTest->bondCount)

      MolOps::Kekulize(*mol);
      if (molTest->reapplyMolBlockWedging) {
        RDKit::Chirality::reapplyMolBlockWedging(*mol);
      }

      {
        std::string expectedMrvName =
            fName + (molTest->sanitizeFlag ? "" : ".nosan") +
            (molTest->reapplyMolBlockWedging ? "" : ".noReapply") +
            ".expected.sdf";
        std::string outMolStr = MolToMolBlock(*mol, true, 0, true, true);

        generateNewExpectedFilesIfSoSpecified(
            fName + (molTest->sanitizeFlag ? "" : ".nosan") +
                (molTest->reapplyMolBlockWedging ? "" : ".noReapply") +
                ".NEW.sdf",
            outMolStr);

        TEST_ASSERT(GetExpectedValue(expectedMrvName) == outMolStr);
      }

      {
        std::string expectedMrvName =
            fName + (molTest->sanitizeFlag ? "" : ".nosan") +
            (molTest->reapplyMolBlockWedging ? "" : ".noReapply") +
            ".expected.mrv";

        std::string outMolStr = MolToMrvBlock(*mol, true, -1, true, false);

        generateNewExpectedFilesIfSoSpecified(
            fName + (molTest->sanitizeFlag ? "" : ".nosan") +
                (molTest->reapplyMolBlockWedging ? "" : ".noReapply") +
                ".NEW.mrv",
            outMolStr);

        TEST_ASSERT(GetExpectedValue(expectedMrvName) == outMolStr);
      }

      BOOST_LOG(rdInfoLog) << "done" << std::endl;
    } catch (const std::exception &) {
      if (molTest->expectedResult != false) {
        throw;
      }

      return;
    }
    TEST_ASSERT(molTest->expectedResult == true);

    return;
  }

  void testMolFiles(const MolTest *molFileTest) {
    BOOST_LOG(rdInfoLog) << "testing marvin writing" << std::endl;
    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase + "/Code/GraphMol/MarvinParse/test_data/" +
                        molFileTest->fileName;

    try {
      std::unique_ptr<RWMol> mol(
          MolFileToMol(fName, molFileTest->sanitizeFlag, false, false));

      if (molFileTest->reapplyMolBlockWedging) {
        Chirality::reapplyMolBlockWedging(*mol);
      }

      TEST_ASSERT(mol != nullptr);
      TEST_ASSERT(mol->getNumAtoms() == molFileTest->atomCount)
      TEST_ASSERT(mol->getNumBonds() == molFileTest->bondCount)

      {
        std::string expectedMrvName =
            fName + (molFileTest->sanitizeFlag ? "" : ".nosan") +
            ".expected.sdf";
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
          outMolStr = MolToMrvBlock(*mol, true, -1, true, false);
        } catch (const RDKit::KekulizeException &) {
          outMolStr = "";
        }
        if (outMolStr == "") {
          // try without kekule'ing
          outMolStr = MolToMrvBlock(*mol, true, -1, false, false);
        }

        generateNewExpectedFilesIfSoSpecified(
            fName + (molFileTest->sanitizeFlag ? "" : ".nosan") + ".NEW.mrv",
            outMolStr);

        TEST_ASSERT(GetExpectedValue(expectedMrvName) == outMolStr);
      }

      {
        std::string expectedSmiName =
            fName + (molFileTest->sanitizeFlag ? "" : ".nosan") +
            ".expected.cxsmi";

        SmilesWriteParams ps;
        ps.cleanStereo = molFileTest->sanitizeFlag;
        ps.canonical = true;
        unsigned int flags = SmilesWrite::CXSmilesFields::CX_COORDS |
                             SmilesWrite::CXSmilesFields::CX_MOLFILE_VALUES |
                             SmilesWrite::CXSmilesFields::CX_ATOM_PROPS |
                             SmilesWrite::CXSmilesFields::CX_BOND_CFG |
                             SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO |
                             SmilesWrite::CXSmilesFields::CX_SGROUPS |
                             SmilesWrite::CXSmilesFields::CX_POLYMER;

        auto restoreDir = RestoreBondDirOptionTrue;
        if (!molFileTest->reapplyMolBlockWedging) {
          restoreDir = RestoreBondDirOptionClear;
        }

        std::string smilesOut = MolToCXSmiles(*mol, ps, flags, restoreDir);

        generateNewExpectedFilesIfSoSpecified(
            fName + (molFileTest->sanitizeFlag ? "" : ".nosan") + ".NEW.cxsmi",
            smilesOut);

        TEST_ASSERT(GetExpectedValue(expectedSmiName) == smilesOut);
      }

      BOOST_LOG(rdInfoLog) << "done" << std::endl;
    } catch (const std::exception &) {
      if (molFileTest->expectedResult != false) {
        throw;
      }
      return;
    }

    TEST_ASSERT(molFileTest->expectedResult == true);

    return;
  }

  void testRxn(const RxnTest *rxnTest) {
    BOOST_LOG(rdInfoLog) << "testing RXN file to Marvin" << std::endl;

    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/MarvinParse/test_data/" + rxnTest->fileName;

    try {
      std::unique_ptr<ChemicalReaction> rxn(
          RxnFileToChemicalReaction(fName, false, false, false));

      // check for errors

      unsigned int nWarn = 0, nError = 0;

      TEST_ASSERT(rxn != nullptr);

      TEST_ASSERT(rxn->getNumReactantTemplates() == rxnTest->reactantCount);
      TEST_ASSERT(rxn->getNumProductTemplates() == rxnTest->productCount);
      TEST_ASSERT(rxn->getNumAgentTemplates() == rxnTest->agentCount);
      rxn->initReactantMatchers(true);

      if (rxn->getNumReactantTemplates() > 0 &&
          rxn->getNumProductTemplates() > 0) {
        TEST_ASSERT(rxn->validate(nWarn, nError, true));
      } else {
        nWarn = 0;
        nError = 0;
      }

      TEST_ASSERT(nWarn == rxnTest->warnings);
      TEST_ASSERT(nError == rxnTest->errors);

      // // make sure the Rxn is kekule'ed

      for (auto mol : rxn->getReactants()) {
        auto rwMol = (RWMol *)mol.get();
        if (rwMol->needsUpdatePropertyCache()) {
          rwMol->updatePropertyCache(false);
        }
        MolOps::Kekulize(*rwMol);
      }
      for (auto mol : rxn->getAgents()) {
        auto rwMol = (RWMol *)mol.get();
        if (rwMol->needsUpdatePropertyCache()) {
          rwMol->updatePropertyCache(false);
        }
        MolOps::Kekulize(*rwMol);
      }
      for (auto mol : rxn->getProducts()) {
        auto rwMol = (RWMol *)mol.get();
        if (rwMol->needsUpdatePropertyCache()) {
          rwMol->updatePropertyCache(false);
        }
        MolOps::Kekulize(*rwMol);
      }

      {
        std::string outMolStr = ChemicalReactionToMrvBlock(*rxn);

        generateNewExpectedFilesIfSoSpecified(fName + ".NEW.mrv", outMolStr);

        std::string expectedRxnName = fName + ".expected.mrv";

        TEST_ASSERT(GetExpectedValue(expectedRxnName) == outMolStr);
      }
      BOOST_LOG(rdInfoLog) << "done" << std::endl;
    } catch (const std::exception &) {
      if (rxnTest->expectedResult != false) {
        throw;
      }
      return;
    }

    TEST_ASSERT(rxnTest->expectedResult == true);

    return;
  }

  void testPrecision() {
    BOOST_LOG(rdInfoLog) << "testing marvin writing" << std::endl;
    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase +
        "/Code/GraphMol/MarvinParse/test_data/precision.mrv.expected.mrv";

    std::unique_ptr<ROMol> mol{v2::MarvinParser::MolFromMrvFile(fName)};

    TEST_ASSERT(mol != nullptr);
    unsigned numAtoms = mol->getNumAtoms();
    TEST_ASSERT(numAtoms == 7);

    MrvWriterParams params{true, true, true, 15};
    std::string outMolStr = MolToMrvBlock(*mol, params, -1);
    std::unique_ptr<ROMol> readMol{
        RDKit::v2::MarvinParser::MolFromMrvBlock(outMolStr)};
    TEST_ASSERT(numAtoms == readMol->getNumAtoms());
    const Conformer &conformer = mol->getConformer();
    const Conformer &readConformer = readMol->getConformer();
    for (size_t i = 0; i < numAtoms; i++) {
      TEST_ASSERT(std::abs(conformer.getAtomPos(i).x -
                           readConformer.getAtomPos(i).x) < 1e-15);
      TEST_ASSERT(std::abs(conformer.getAtomPos(i).y -
                           readConformer.getAtomPos(i).y) < 1e-15);
      TEST_ASSERT(std::abs(conformer.getAtomPos(i).z -
                           readConformer.getAtomPos(i).z) < 1e-15);
    }
    BOOST_LOG(rdInfoLog) << "done" << std::endl;
  }

  void testMarvinRgroupError() {
    BOOST_LOG(rdInfoLog) << "testing marvin Rgroup error" << std::endl;
    const auto molb = R"CTAB( 
  -INDIGO-06252516172D

  2  1  0  0  0  0  0  0  0  0999 V2000
    9.9900   -7.1250    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
    8.9900   -7.1250    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  RGP  1   1   1
M  END
    )CTAB";
    std::unique_ptr<ROMol> mol = v2::FileParsers::MolFromMolBlock(molb);

    SmilesWriteParams ps;
    ps.canonical = false;

    std::string smilesOut = MolToSmiles(*mol, ps);
    std::string cxSmilesOut = MolToCXSmiles(*mol, ps);

    auto outMolStr = MolToMrvBlock(*mol, true, -1, true, false);

    std::unique_ptr<ROMol> mol2 =
        RDKit::v2::MarvinParser::MolFromMrvBlock(outMolStr);

    std::string smilesOut2 = MolToSmiles(*mol2, ps);
    std::string cxSmilesOut2 = MolToCXSmiles(*mol2, ps);

    TEST_ASSERT(smilesOut == smilesOut2);
    TEST_ASSERT(cxSmilesOut == cxSmilesOut2);
  }

  void testMarvinRgroupError2() {
    BOOST_LOG(rdInfoLog)
        << "testing marvin Rgroup error starting with MRV block" << std::endl;
    const auto mrvb =
        R"CTAB(<cml xmlns="http://www.chemaxon.com" version="ChemAxon file format v20.20.0, generated by vunknown" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.chemaxon.com http://www.chemaxon.com/marvin/schema/mrvSchema_20_20_0.xsd">
    <MDocument>
        <MChemicalStruct>
            <molecule molID="m1">
                <atomArray>
                    <atom id="a1" elementType="R" x2="-0.666834670524129" y2="-0.38500847027952034" rgroupRef="1"/>
                    <atom id="a2" elementType="O" x2="0.666834670524129" y2="0.38500847027952023" lonePair="2"/>
                </atomArray>
                <bondArray>
                    <bond id="b1" atomRefs2="a1 a2" order="1"/>
                </bondArray>
            </molecule>
        </MChemicalStruct>
    </MDocument>
</cml>)CTAB";

    std::unique_ptr<ROMol> mol = v2::MarvinParser::MolFromMrvBlock(mrvb);
    std::string tempStr;
    auto firstAtom = mol->getAtomWithIdx(0);

    TEST_ASSERT(
        firstAtom->getProp<std::string>(common_properties::dummyLabel) == "R1");
    TEST_ASSERT(firstAtom->getIsotope() == 1);
    int iLabel(0);
    firstAtom->getProp(common_properties::_MolFileRLabel, iLabel);
    TEST_ASSERT(iLabel == 1);
  }

 public:
  void RunTests() {
    UseLegacyStereoPerceptionFixture useLegacy(false);
    printf("Using new chirality perception\n");

    auto sanitizeOff = false;
    auto reapplyWedgesOn = true;
    auto sanitizeOn = true;
    auto reapplyWedgesOff = false;

    // rxn test returning a single mol

    if (testToRun == "" || testToRun == "rxnMolTests") {
      std::list<RxnTest> rxnMolTests{
          RxnTest("rxnStereoMarkedCrossed.mrv", true, 1, 0, 1, 2, 0),
      };

      for (auto rxnMolTest : rxnMolTests) {
        BOOST_LOG(rdInfoLog) << "Test: " << rxnMolTest.fileName << std::endl;

        printf("Test\n\n %s\n\n", rxnMolTest.fileName.c_str());
        testMarvinRxnMols(&rxnMolTest);
      }
    }

    if (testToRun == "" || testToRun == "testMarvinRgroupError") {
      testMarvinRgroupError();
    }
    if (testToRun == "" || testToRun == "testMarvinRgroupError2") {
      testMarvinRgroupError2();
    }

    // the molecule tests - starting with molfiles/sdf
    if (testToRun == "" || testToRun == "sdfTests") {
      std::list<MolTest> sdfTests{
          MolTest("145323811.mol", true, 172, 176, true, false),
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
          MolTest("vendor839.sdf", true, 25, 26),
          MolTest("ProblemShort.mol", true, 993, 992),
          MolTest("badWedgeError.sdf", true, 12, 13),
          MolTest("CrossedDoubleBondWithChiralNbr2.sdf", true, 10, 9),
          MolTest("CrossedDoubleBondWithChiralNbr.sdf", true, 10, 9),
          MolTest("SimpleWiggleDoubleBond.sdf", true, 6, 5),
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

    // now the RXN reactions

    if (testToRun == "" || testToRun == "rxnFileTests") {
      std::list<RxnTest> rxnFileTests{
          RxnTest("BadRxn.rxn", true, 2, 0, 1, 3, 0),
      };

      for (auto rxnFileTest : rxnFileTests) {
        printf("Test\n\n %s\n\n", rxnFileTest.fileName.c_str());
        testRxn(&rxnFileTest);
      }
    }

    if (testToRun == "" || testToRun == "molFileTests") {
      std::list<MolTest> molFileTests{
          MolTest("FalseChiral.mrv", true, 4, 3, sanitizeOff,
                  reapplyWedgesOff),  // not sanitized, wedges NOT reapplied
          MolTest("FalseChiral.mrv", true, 4, 3, sanitizeOn,
                  reapplyWedgesOff),  // sanitized, wedges NOT reapplied

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

          MolTest("DataSgroupMissingUnitsDisplayed.mrv", true, 15, 16),
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
                  false),  // not sanitized
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

    if (testToRun == "" || testToRun == "scsrFileTests") {
      std::list<ScsrMolTest> scsrFileTests{
          ScsrMolTest("153944501_original_structure.mol", true, 859, 1025, 128),
          ScsrMolTest("rnaTest.mol", true, 22, 24, 4)};

      for (auto &scsrFileTest : scsrFileTests) {
        BOOST_LOG(rdInfoLog) << "Test: " << scsrFileTest.fileName << std::endl;

        printf("Test\n\n %s\n\n", scsrFileTest.fileName.c_str());
        testMarvinFromScsr(&scsrFileTest);
      }
    }

    if (testToRun == "" || testToRun == "chiral3dFileTests") {
      std::list<MolTest> chiral3dFileTests{
          MolTest("Cubane.mrv", true, 16, 20),
      };

      for (auto &molFileTest : chiral3dFileTests) {
        BOOST_LOG(rdInfoLog) << "Test: " << molFileTest.fileName << std::endl;

        printf("Test\n\n %s\n\n", molFileTest.fileName.c_str());
        testMarvin3dChiral(&molFileTest);
      }
    }

    // atropisomer tests
    if (testToRun == "" || testToRun == "atropisomerTests") {
      std::vector<MolTest> atropisomerTests{
          // first the tests with sanitize off,
          // and reapplyMolBlockWedging on
          MolTest("FalseAtropisomer.mrv", true, 17, 18, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("AtropEnhancedStereo.mrv", true, 16, 17, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("AtropManyChirals.mrv", true, 20, 20, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("AtropManyChiralsEnhanced.mrv", true, 20, 20, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("AtropManyChiralsEnhanced2.mrv", true, 20, 20, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("AtropManyChiralsEnhanced3.mrv", true, 20, 20, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("AtropManyChiralsEnhanced4.mrv", true, 20, 20, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("BMS-986142_3d_chiral.mrv", true, 72, 77, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("BMS-986142_3d.mrv", true, 72, 77, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("BMS-986142_atrop1.mrv", true, 42, 47, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("BMS-986142_atrop2.mrv", true, 42, 47, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("BMS-986142_atrop3.mrv", true, 42, 47, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("BMS-986142_atrop4.mrv", true, 42, 47, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("BMS-986142_atrop5.mrv", true, 42, 47, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("BMS-986142_atrop6.mrv", true, 42, 47, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("BMS-986142_atrop7.mrv", true, 42, 47, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("BMS-986142_atrop8.mrv", true, 42, 47, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("BMS-986142_atropBad2.mrv", true, 42, 47, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("JDQ443_3d.mrv", true, 66, 72, sanitizeOff, reapplyWedgesOn),
          MolTest("JDQ443_atrop1.mrv", true, 38, 44, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("JDQ443_atrop2.mrv", true, 38, 44, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("JDQ443_atrop3.mrv", true, 38, 44, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("JDQ443_atropBad1.mrv", true, 38, 44, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("RP-6306_atrop1.mrv", true, 24, 26, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("RP-6306_atrop2.mrv", true, 24, 26, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("RP-6306_atrop3.mrv", true, 24, 26, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("RP-6306_atrop4.mrv", true, 24, 26, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("RP-6306_atrop5.mrv", true, 24, 26, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("RP-6306_atropBad1.mrv", true, 24, 26, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("RP-6306_atropBad2.mrv", true, 24, 26, sanitizeOff,
                  reapplyWedgesOn),
          // note the rp-6306_3d.mrv is backwards from the 2D versions
          // the 2D version were based on images from drug hunter
          // the 3D version came from PUBCHEM
          MolTest("RP-6306_3d.mrv", true, 44, 46, sanitizeOff, reapplyWedgesOn),
          MolTest("Sotorasib_atrop1.mrv", true, 41, 45, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("Sotorasib_atrop2.mrv", true, 41, 45, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("Sotorasib_atrop3.mrv", true, 41, 45, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("Sotorasib_atrop4.mrv", true, 41, 45, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("Sotorasib_atrop5.mrv", true, 41, 45, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("Sotorasib_atropBad1.mrv", true, 41, 45, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("Sotorasib_atropBad2.mrv", true, 41, 45, sanitizeOff,
                  reapplyWedgesOn),
          // note the sotorasib_3d.mrv is backwards from the 2D versions
          // the 2D version were based on images from drug hunter
          // the 3D version came from PUBCHEM
          MolTest("Sotorasib_3d.mrv", true, 71, 75, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("ZM374979_atrop1.mrv", true, 45, 49, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("ZM374979_atrop2.mrv", true, 45, 49, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("ZM374979_atrop3.mrv", true, 45, 49, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("ZM374979_atropBad1.mrv", true, 45, 49, sanitizeOff,
                  reapplyWedgesOn),
          // note the mrtx1719_3d.mrv is backwards from the 2D versions
          // the 2D version were based on images from drug hunter
          // the 3D version came from PUBCHEM
          MolTest("mrtx1719_3d.mrv", true, 51, 55, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("mrtx1719_atrop1.mrv", true, 33, 37, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("mrtx1719_atrop2.mrv", true, 33, 37, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("mrtx1719_atrop3.mrv", true, 33, 37, sanitizeOff,
                  reapplyWedgesOn),
          MolTest("mrtx1719_atropBad1.mrv", true, 33, 37, sanitizeOff,
                  reapplyWedgesOn),

          // now the tests with sanitizeOff on,
          // and reapplyMolBlockWedging off

          MolTest("FalseAtropisomer.mrv", true, 17, 18, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("AtropEnhancedStereo.mrv", true, 16, 17, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("AtropManyChirals.mrv", true, 20, 20, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("AtropManyChiralsEnhanced.mrv", true, 20, 20, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("AtropManyChiralsEnhanced2.mrv", true, 20, 20, true, false),
          MolTest("AtropManyChiralsEnhanced3.mrv", true, 20, 20, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("AtropManyChiralsEnhanced4.mrv", true, 20, 20, true, false),
          MolTest("BMS-986142_3d_chiral.mrv", true, 72, 77, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("BMS-986142_3d.mrv", true, 72, 77, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("BMS-986142_atrop1.mrv", true, 42, 47, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("BMS-986142_atrop2.mrv", true, 42, 47, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("BMS-986142_atrop3.mrv", true, 42, 47, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("BMS-986142_atrop4.mrv", true, 42, 47, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("BMS-986142_atrop5.mrv", true, 42, 47, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("BMS-986142_atrop6.mrv", true, 42, 47, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("BMS-986142_atrop7.mrv", true, 42, 47, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("BMS-986142_atrop8.mrv", true, 42, 47, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("BMS-986142_atropBad2.mrv", true, 42, 47, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("JDQ443_3d.mrv", true, 66, 72, sanitizeOn, reapplyWedgesOff),
          MolTest("JDQ443_atrop1.mrv", true, 38, 44, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("JDQ443_atrop2.mrv", true, 38, 44, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("JDQ443_atrop3.mrv", true, 38, 44, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("JDQ443_atropBad1.mrv", true, 38, 44, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("RP-6306_atrop1.mrv", true, 24, 26, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("RP-6306_atrop2.mrv", true, 24, 26, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("RP-6306_atrop3.mrv", true, 24, 26, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("RP-6306_atrop4.mrv", true, 24, 26, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("RP-6306_atrop5.mrv", true, 24, 26, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("RP-6306_atropBad1.mrv", true, 24, 26, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("RP-6306_atropBad2.mrv", true, 24, 26, sanitizeOn,
                  reapplyWedgesOff),
          // note the rp-6306_3d.mrv is backwards from the 2D versions
          // the 2D version were based on images from drug hunter
          // the 3D version came from PUBCHEM
          MolTest("RP-6306_3d.mrv", true, 44, 46, sanitizeOn, reapplyWedgesOff),
          MolTest("Sotorasib_atrop1.mrv", true, 41, 45, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("Sotorasib_atrop2.mrv", true, 41, 45, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("Sotorasib_atrop3.mrv", true, 41, 45, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("Sotorasib_atrop4.mrv", true, 41, 45, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("Sotorasib_atrop5.mrv", true, 41, 45, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("Sotorasib_atropBad1.mrv", true, 41, 45, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("Sotorasib_atropBad2.mrv", true, 41, 45, sanitizeOn,
                  reapplyWedgesOff),
          // note the sotorasib_3d.mrv is backwards from the 2D versions
          // the 2D version were based on images from drug hunter
          // the 3D version came from PUBCHEM
          MolTest("Sotorasib_3d.mrv", true, 71, 75, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("ZM374979_atrop1.mrv", true, 45, 49, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("ZM374979_atrop2.mrv", true, 45, 49, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("ZM374979_atrop3.mrv", true, 45, 49, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("ZM374979_atropBad1.mrv", true, 45, 49, sanitizeOn,
                  reapplyWedgesOff),
          // note the mrtx1719_3d.mrv is backwards from the 2D versions
          // the 2D version were based on images from drug hunter
          // the 3D version came from PUBCHEM
          MolTest("mrtx1719_3d.mrv", true, 51, 55, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("mrtx1719_atrop1.mrv", true, 33, 37, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("mrtx1719_atrop2.mrv", true, 33, 37, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("mrtx1719_atrop3.mrv", true, 33, 37, sanitizeOn,
                  reapplyWedgesOff),
          MolTest("mrtx1719_atropBad1.mrv", true, 33, 37, sanitizeOn,
                  reapplyWedgesOff),
      };

      for (auto &atropisomerTest : atropisomerTests) {
        BOOST_LOG(rdInfoLog)
            << "Test: " << atropisomerTest.fileName << std::endl;

        printf("Test\n\n %s\n\n", atropisomerTest.fileName.c_str());
        testMarvinAtrop(&atropisomerTest);
      }
    }
    // now the reactions

    if (testToRun == "" || testToRun == "mrvRxnFileTests") {
      std::list<RxnTest> mrvRxnFileTests{
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

      for (auto rxnFileTest : mrvRxnFileTests) {
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

          SmilesTest(
              "Vendor839",
              R"(CCOc1cc(/C=C2\C(=O)OC(C)(C(C)(C)C)OC2=O)cc(Br)c1O |(2.481,0.6126,;1.7562,0.2127,;1.7438,-0.6123,;1.0313,-1.0125,;0.3188,-0.5874,;-0.3937,-0.9873,;-1.0979,-0.5708,;-1.1062,0.25,;-0.4061,0.6544,;0.2812,0.25,;-0.4061,1.4793,;-1.1062,1.8875,;-1.6937,2.475,;-0.5186,2.475,;0.0645,3.0623,;0.0645,1.8999,;-1.1062,3.0623,;-1.7935,1.4793,;-1.7935,0.6544,;-2.481,0.25,;-0.4061,-1.8122,;0.2936,-2.2373,;0.2812,-3.0623,;1.0189,-1.8375,;1.7186,-2.2622,)|)",
              true, 25, 26),
          SmilesTest("DoubleBondChain",
                     R"(CC1=C(\C=C\C(C)=C\C=C\C(C)=C/C(O)=O)C(C)(C)CCC1)", true,
                     22, 22),
          // this does NOT work - still working on a solution for then

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

      for (auto &smiTest : smiTests) {
        printf("Test\n\n %s\n\n", smiTest.name.c_str());
        testSmilesToMarvin(&smiTest);
      }
    }
    if (testToRun == "" || testToRun == "precision") {
      testPrecision();
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

  mrvTests.RunTests();

  return 0;
}
