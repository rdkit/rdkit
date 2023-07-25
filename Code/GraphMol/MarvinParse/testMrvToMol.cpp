//
//  Copyright (C) 2002-2021 Greg Landrum and other RDKit contributors
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
#include "MarvinParser.h"

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <string>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <filesystem>

using namespace RDKit;

enum LoadAs { LoadAsMolOrRxn, LoadAsMol, LoadAsRxn };

class MolOrRxnTest {
public:
    std::string fileName;
    bool expectedResult;
    LoadAs loadAs;

    MolOrRxnTest(std::string fileNameInit, bool expectedResultInit,
                 LoadAs loadAsInit)
            : fileName(fileNameInit),
              expectedResult(expectedResultInit),
              loadAs(loadAsInit){};

    virtual bool isRxnTest() const = 0;
};

class MolTest : public MolOrRxnTest {
public:
    unsigned int atomCount;
    unsigned int bondCount;

    MolTest(std::string fileNameInit, bool expectedResultInit, LoadAs loadAsInit,
            int atomCountInit, int bondCountInit)
            : MolOrRxnTest(fileNameInit, expectedResultInit, loadAsInit),
              atomCount(atomCountInit),
              bondCount(bondCountInit){};

    bool isRxnTest() const override { return false; }
};

class RxnTest : public MolOrRxnTest {
public:
    unsigned int reactantCount;
    unsigned int agentCount;
    unsigned int productCount;
    unsigned int warnings;
    unsigned int errors;

    RxnTest(std::string fileNameInit, bool expectedResultInit, LoadAs loadAsInit,
            int reactantCountInit, int agentCountInit, int productCountInit,
            int warnInit, int errorInit)
            : MolOrRxnTest(fileNameInit, expectedResultInit, loadAsInit),
              reactantCount(reactantCountInit),
              agentCount(agentCountInit),
              productCount(productCountInit),
              warnings(warnInit),
              errors(errorInit){};

    bool isRxnTest() const override { return true; }
};

class SmilesTest {
public:
    std::string name;
    std::string smiles;
    bool expectedResult;
    unsigned int atomCount;
    unsigned int bondCount;

    SmilesTest(std::string nameInit, std::string smilesInit,
               bool expectedResultInit, int atomCountInit, int bondCountInit)
            : name(nameInit),
              smiles(smilesInit),
              expectedResult(expectedResultInit),
              atomCount(atomCountInit),
              bondCount(bondCountInit){};

    bool isRxnTest() const { return false; }
};

void *GetMolOrReaction(const MolOrRxnTest *molOrRxnTest, bool &isReaction) {
  std::string rdbase = getenv("RDBASE");
  std::string fName =
          rdbase + "/Code/GraphMol/MarvinParse/test_data/" + molOrRxnTest->fileName;

  bool sanitize;
  int tryCount;
  for (tryCount = 0, sanitize = true; tryCount < 2;
       ++tryCount,
               sanitize =
                       false)  // try with sanitize on - if it fails try with sanitize off
  {
    try {
      switch (molOrRxnTest->loadAs) {
        case LoadAsMolOrRxn:
          return MrvFileParser(fName, isReaction, sanitize, false);

        case LoadAsMol:
          isReaction = false;
          return (void *)MrvMolFileParser(fName, sanitize, false);

        case LoadAsRxn:
          isReaction = true;
          return (void *)MrvRxnFileParser(fName, sanitize, false);
      }
    } catch (const std::exception &e) {
      std::cerr << e.what() << '\n';
    }
  }

  throw BadFileException("Could not parse the MRV block");
}

void testSmilesToMarvin(const SmilesTest *smilesTest) {
  BOOST_LOG(rdInfoLog) << "testing smiles to marin " << std::endl;
  std::string rdbase = getenv("RDBASE");
  std::string fName =
          rdbase + "/Code/GraphMol/MarvinParse/test_data/" + smilesTest->name;

  class LocalVars  // protect against mem leak on error
  {
  public:
      RWMol *smilesMol;

      LocalVars(){};

      ~LocalVars() { delete smilesMol; }
  } localVars;

  try {
    SmilesParserParams smilesParserParams;
    smilesParserParams.sanitize = true;

    localVars.smilesMol = SmilesToMol(smilesTest->smiles, smilesParserParams);
    reapplyMolBlockWedging(*localVars.smilesMol);

    TEST_ASSERT(localVars.smilesMol->getNumAtoms() == smilesTest->atomCount);
    TEST_ASSERT(localVars.smilesMol->getNumBonds() == smilesTest->bondCount);

    // test round trip back to smiles
    {
      std::string expectedMrvName = fName + ".expected.smi";

      SmilesWriteParams ps;
      ps.canonical = false;

      std::string smilesOut = MolToSmiles(*localVars.smilesMol, ps);

      // code to generate the expected files

      // {
      //   std::ofstream out;
      //   out.open(fName + ".NEW.smi");
      //   out << smilesOut;
      // }
      std::stringstream expectedMolStr;
      std::ifstream in;
      in.open(expectedMrvName);
      expectedMolStr << in.rdbuf();
      std::string expectedStr = expectedMolStr.str();

      TEST_ASSERT(expectedStr == smilesOut);
    }

    {
      std::string expectedMrvName = fName + ".expected.sdf";
      std::string outMolStr = "";
      try {
        outMolStr = MolToMolBlock(*localVars.smilesMol, true, 0, true, true);
      } catch (const RDKit::KekulizeException &e) {
        outMolStr = "";
      } catch (...) {
        throw;  // re-trhow the error if not a kekule error
      }
      if (outMolStr == "") {
        outMolStr = MolToMolBlock(*localVars.smilesMol, true, 0, false,
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

      TEST_ASSERT(expectedStr == outMolStr);
    }
    {
      std::string mrvBlock;
      std::string expectedMrvName = fName + ".expected.mrv";
      std::string outMolStr = "";
      try {
        outMolStr = MolToMrvBlock(*localVars.smilesMol, true, -1, true, false);
      } catch (const RDKit::KekulizeException &e) {
        outMolStr = "";
      } catch (...) {
        throw;  // re-trhow the error if not a kekule error
      }
      if (outMolStr == "") {
        outMolStr = MolToMrvBlock(*localVars.smilesMol, true, -1, false,
                                  false);  // try without kekule'ing
      }

      // code to generate the expected files
      // {
      //   std::ofstream out;
      //   out.open(fName + ".NEW.mrv");
      //   out << outMolStr;
      // }
      std::stringstream expectedMolStr;
      std::ifstream in;
      in.open(expectedMrvName);
      expectedMolStr << in.rdbuf();
      std::string expectedStr = expectedMolStr.str();

      TEST_ASSERT(expectedStr == outMolStr);
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

void testMarvin(const MolOrRxnTest *molOrRxnTest) {
  BOOST_LOG(rdInfoLog) << "testing marvin parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName =
          rdbase + "/Code/GraphMol/MarvinParse/test_data/" + molOrRxnTest->fileName;

  class LocalVars  // protext against mem leak on error
  {
  public:
      void *molOrRxn = nullptr;
      bool isReaction = false;

      LocalVars(){};

      ~LocalVars() {
        // printf ("deleting molOrRxn\n");
        if (isReaction) {
          delete (ChemicalReaction *)molOrRxn;
        } else {
          delete (RWMol *)molOrRxn;
        }
      }
  } localVars;

  try {
    localVars.molOrRxn = GetMolOrReaction(molOrRxnTest, localVars.isReaction);

    if (localVars.isReaction != molOrRxnTest->isRxnTest()) {
      // printf("Wrong type of MRV file\n");
      TEST_ASSERT(molOrRxnTest->expectedResult == false);
      // printf("Expected failure!\n");
      return;
    }

    if (localVars.isReaction) {
      // reaction test

      auto rxn = (ChemicalReaction *)localVars.molOrRxn;
      auto rxnTest = (RxnTest *)molOrRxnTest;

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

        //  code to create the expected files for new or changed tests

        // {
        //   std::ofstream out;
        //   out.open(fName + ".NEW.rxn");
        //   out << outMolStr;
        // }

        std::string expectedRxnName = fName + ".expected.rxn";

        std::stringstream expectedMolStr;
        std::ifstream in;
        in.open(expectedRxnName);
        expectedMolStr << in.rdbuf();
        std::string expectedStr = expectedMolStr.str();

        TEST_ASSERT(expectedStr == outMolStr);
      }

      {
        std::string outMolStr = ChemicalReactionToMrvBlock(*rxn, false);

        //  code to create the expected files for new or changed tests

        //{
        //   std::ofstream out;
        //   out.open(fName + ".NEW.mrv");
        //   out << outMolStr;
        // }

        std::string expectedRxnName = fName + ".expected.mrv";
        std::stringstream expectedMolStr;
        std::ifstream in;
        in.open(expectedRxnName);
        expectedMolStr << in.rdbuf();
        std::string expectedStr = expectedMolStr.str();

        TEST_ASSERT(expectedStr == outMolStr);
      }
      BOOST_LOG(rdInfoLog) << "done" << std::endl;
    } else {
      // mol  test

      auto mol = (RWMol *)localVars.molOrRxn;
      reapplyMolBlockWedging(*mol);

      auto molTest = (MolTest *)molOrRxnTest;
      TEST_ASSERT(mol != nullptr);

      TEST_ASSERT(mol->getNumAtoms() == molTest->atomCount)
      TEST_ASSERT(mol->getNumBonds() == molTest->bondCount)

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

        //{
        //    std::ofstream out;
        //   out.open(fName + ".NEW.sdf");
        // out << outMolStr;
        // }

        std::stringstream expectedMolStr;
        std::ifstream in;
        in.open(expectedMrvName);
        expectedMolStr << in.rdbuf();
        std::string expectedStr = expectedMolStr.str();
        
        TEST_ASSERT(expectedStr == outMolStr);
      }

      {
        std::string expectedMrvName = fName + ".expected.mrv";

        std::string outMolStr = "";
        try {
          outMolStr = MolToMrvBlock(*mol, true, -1, true, false);
        } catch (const RDKit::KekulizeException &e) {
          outMolStr = "";
        } catch (...) {
          throw;  // re-trhow the error if not a kekule error
        }
        if (outMolStr == "") {
          outMolStr = MolToMrvBlock(*mol, true, -1, false,
                                    false);  // try without kekule'ing
        }
        // code to create the expected files for new or changed tests

        // {
        //   std::ofstream out;
        //   out.open(fName + ".NEW.mrv");
        //   out << outMolStr;
        // }

        std::stringstream expectedMolStr;
        std::ifstream in;
        in.open(expectedMrvName);
        expectedMolStr << in.rdbuf();
        std::string expectedStr = expectedMolStr.str();

        TEST_ASSERT(expectedStr == outMolStr);
      }

      BOOST_LOG(rdInfoLog) << "done" << std::endl;
    }
  } catch (const std::exception &e) {
    if (molOrRxnTest->expectedResult != false) {
      throw;
    }
    return;
  }

  TEST_ASSERT(molOrRxnTest->expectedResult == true);

  return;
}

void testMolFiles(const MolTest *molFileTest) {
  BOOST_LOG(rdInfoLog) << "testing marvin writing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName =
          rdbase + "/Code/GraphMol/MarvinParse/test_data/" + molFileTest->fileName;

  class LocalVars  // protext against mem leak on error
  {
  public:
      RWMol *mol = nullptr;

      LocalVars(){};

      ~LocalVars() { delete mol; }
  } localVars;

  try {
    localVars.mol = MolFileToMol(fName, true, false, false);
    reapplyMolBlockWedging(*localVars.mol);

    TEST_ASSERT(localVars.mol != nullptr);
    TEST_ASSERT(localVars.mol->getNumAtoms() == molFileTest->atomCount)
    TEST_ASSERT(localVars.mol->getNumBonds() == molFileTest->bondCount)

    {
      std::string expectedMrvName = fName + ".expected.sdf";
      std::string outMolStr = "";
      try {
        outMolStr = MolToMolBlock(*localVars.mol, true, 0, true, true);
      } catch (const RDKit::KekulizeException &e) {
        outMolStr = "";
      } catch (...) {
        throw;  // re-trhow the error if not a kekule error
      }
      if (outMolStr == "") {
        outMolStr = MolToMolBlock(*localVars.mol, true, 0, false,
                                  true);  // try without kekule'ing
      }

      // code to create the expected files for new or changed tests

      {
        std::ofstream out;
        out.open(fName + ".NEW.sdf");
        out << outMolStr;
      }

      std::stringstream expectedMolStr;
      std::ifstream in;
      in.open(expectedMrvName);
      expectedMolStr << in.rdbuf();
      std::string expectedStr = expectedMolStr.str();
      if (expectedStr != outMolStr) {
        std::cerr << "ERROR: " << fName << " failed" << std::endl;
      }
      TEST_ASSERT(expectedStr == outMolStr);
    }

    {
      std::string expectedMrvName = fName + ".expected.mrv";

      std::string outMolStr = "";
      try {
        outMolStr = MolToMrvBlock(*localVars.mol, true, -1, true, false);
      } catch (const RDKit::KekulizeException &e) {
        outMolStr = "";
      } catch (...) {
        throw;  // re-trhow the error if not a kekule error
      }
      if (outMolStr == "") {
        // try without kekule'ing
        outMolStr = MolToMrvBlock(*localVars.mol, true, -1, false, false);
      }
      // code to create the expected files for new or changed tests

      // {
      //   std::ofstream out;
      //   out.open(fName + ".NEW.mrv");
      //   out << outMolStr;
      // }

      std::stringstream expectedMolStr;
      std::ifstream in;
      in.open(expectedMrvName);
      expectedMolStr << in.rdbuf();
      std::string expectedStr = expectedMolStr.str();

      TEST_ASSERT(expectedStr == outMolStr);
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
  // the molecule tests - starting with molfiles/sdf

  std::list<MolTest> sdfTests{
          MolTest("ProblemShort.mol", true, LoadAsMolOrRxn, 993, 992),
          MolTest("lostStereoAnd.sdf", true, LoadAsMolOrRxn, 6, 5),
          MolTest("DoubleBondChain.sdf", true, LoadAsMolOrRxn, 22, 22),
          MolTest("UnitsError.sdf", true, LoadAsMolOrRxn, 17, 18),
          MolTest("StarAtom.sdf", true, LoadAsMolOrRxn, 17, 16),
          MolTest("nonProprietary.mol", true, LoadAsMolOrRxn, 17, 16),
          MolTest("ChiralTest.sdf", true, LoadAsMolOrRxn, 8, 7),
          MolTest("TestBond1.mol", true, LoadAsMolOrRxn, 10, 10),
          MolTest("Sgroup_MUL_ParentInMiddle.sdf", true, LoadAsMolOrRxn, 17, 16)};

  for (auto sdfTest : sdfTests) {
    BOOST_LOG(rdInfoLog) << "Test: " << sdfTest.fileName << std::endl;

    printf("Test\n\n %s\n\n", sdfTest.fileName.c_str());
    testMolFiles(&sdfTest);
  }

  std::list<MolTest> molFileTests{
          MolTest("lostStereoAnd.mrv", true, LoadAsMolOrRxn, 6, 5),
          MolTest("DoubleBondChain.mrv", true, LoadAsMolOrRxn, 22, 22),
          MolTest("WigglyAndCrossed.mrv", true, LoadAsMolOrRxn, 8, 7),
          MolTest("BondTypes.mrv", true, LoadAsMolOrRxn, 26, 25),
          MolTest("EmbeddedSGroupSUP_MUL.mrv", true, LoadAsMolOrRxn, 17, 17),
          MolTest("EmbeddedSgroupCOP_SUP.mrv", true, LoadAsMolOrRxn, 10, 10),
          MolTest("EmbeddedSgroupDAT_SUP.mrv", true, LoadAsMolOrRxn, 10, 10),
          MolTest("EmbeddedSgroupMULTICENTER_SUP.mrv", true, LoadAsMolOrRxn, 11,
                  10),
          MolTest("EmbeddedSgroupMUL_MUL.mrv", true, LoadAsMolOrRxn, 141, 140),
          MolTest("EmbeddedSgroupMUL_MUL2.mrv", true, LoadAsMolOrRxn, 23, 22),
          MolTest("EmbeddedSgroupMUL_SUP.mrv", true, LoadAsMolOrRxn, 129, 128),
          MolTest("EmbeddedSgroupSRU_SUP.mrv", true, LoadAsMolOrRxn, 10, 10),
          MolTest("EmbeddedSgroupSUPEXP_SUP.mrv", true, LoadAsMolOrRxn, 10, 10),
          MolTest("EmbeddedSgroupSUPEXP_SUP2.mrv", true, LoadAsMolOrRxn, 10, 10),
          MolTest("EmbeddedSgroupSUP_MULTICENTER.mrv", true, LoadAsMolOrRxn, 10, 8),
          MolTest("EmbeddedSgroupSUP_SUP.mrv", true, LoadAsMolOrRxn, 12, 11),
          MolTest("EmbeddedSgroupSUP_SUP2.mrv", true, LoadAsMolOrRxn, 12, 12),
          MolTest("RgroupBad.mrv", true, LoadAsMolOrRxn, 9, 9),
          MolTest("valenceLessThanDrawn.mrv", true, LoadAsMolOrRxn, 14, 14),
          MolTest("data_sgroup_no_fieldname.mrv", true, LoadAsMolOrRxn, 4, 3),
          MolTest("data_sgroup_empty_field_data.mrv", true, LoadAsMolOrRxn, 2, 1),
          MolTest("radical_value.mrv", true, LoadAsMolOrRxn, 3, 2),
          MolTest("emptyOneLineAtomList.mrv", true, LoadAsMolOrRxn, 0, 0),
          MolTest("mrvValence_value.mrv", true, LoadAsMolOrRxn, 3, 2),
          MolTest("ChiralTest2.mrv", true, LoadAsMolOrRxn, 46, 47),
          MolTest("ChiralTest.mrv", true, LoadAsMolOrRxn, 8, 7),
          MolTest("SnCl2.mrv", true, LoadAsMolOrRxn, 3, 2),
          MolTest("SnH2Cl2.mrv", true, LoadAsMolOrRxn, 3, 2),
          MolTest("marvin01.mrv", true, LoadAsMolOrRxn, 11, 11),
          MolTest("marvin01.mrv", true, LoadAsMol, 11, 11),
          MolTest("marvin01.mrv", false, LoadAsRxn, 11, 11)  // should fail
          ,
          MolTest("marvin02.mrv", true, LoadAsMolOrRxn, 9, 9),
          MolTest("marvin07.mrv", true, LoadAsMolOrRxn, 12, 11),
          MolTest("marvin10.mrv", true, LoadAsMolOrRxn, 10, 10),
          MolTest("marvin06.mrv", true, LoadAsMolOrRxn, 11, 11),
          MolTest("marvin12.mrv", true, LoadAsMolOrRxn, 31, 33),
          MolTest("EmptyMol.mrv", true, LoadAsMolOrRxn, 0, 0),
          MolTest("Sparse.mrv", true, LoadAsMolOrRxn, 0, 0),
          MolTest("Sparse2.mrv", true, LoadAsMolOrRxn, 0, 0),
          MolTest("Sparse3.mrv", true, LoadAsMolOrRxn, 0, 0),
          MolTest("MarvinNoCoords.mrv", true, LoadAsMolOrRxn, 6, 6),
          MolTest("aspirin.mrv", true, LoadAsMolOrRxn, 13, 13),
          MolTest("MarvinStereoGroupsZeros.mrv", true, LoadAsMolOrRxn, 8, 8),
          MolTest("MarvinStereoGroupsAbs.mrv", true, LoadAsMolOrRxn, 8, 8),
          MolTest("triphenylphosphine.mrv", true, LoadAsMolOrRxn, 19, 21),
          MolTest("MarvinOldSuperGroupTest.mrv", true, LoadAsMolOrRxn, 89, 93),
          MolTest("RadicalTests.mrv", true, LoadAsMolOrRxn, 8, 7),
          MolTest("AnyBond.mrv", true, LoadAsMolOrRxn, 4, 3),
          MolTest("cisBenzene.mrv", true, LoadAsMolOrRxn, 6, 6),
          MolTest("DativeBond.mrv", true, LoadAsMolOrRxn, 6, 5),
          MolTest("MultipleSgroup.mrv", true, LoadAsMolOrRxn, 123, 122),
          MolTest("SgroupExpanded.mrv", true, LoadAsMolOrRxn, 5, 4),
          MolTest("SgroupMultAttach.mrv", true, LoadAsMolOrRxn, 44, 45),
          MolTest("MarvinMissingX2.mrv", true, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinMissingY2.mrv", true, LoadAsMolOrRxn, 12, 11),
          MolTest("DataSgroup.mrv", true, LoadAsMolOrRxn, 7, 6),
          MolTest("MulticenterSgroup.mrv", true, LoadAsMolOrRxn, 17, 16),
          MolTest("GenericSgroup.mrv", true, LoadAsMolOrRxn, 13, 13),
          MolTest("MonomerSgroup.mrv", true, LoadAsMolOrRxn, 4, 3),
          MolTest("modification_sgroup.mrv", true, LoadAsMolOrRxn, 54, 40),
          MolTest("copolymer_sgroup.mrv", true, LoadAsMolOrRxn, 19, 18),
          MolTest("MultipleSgroupParentInMiddleOfAtomBlock.mrv", true,
                  LoadAsMolOrRxn, 23, 22),
          MolTest("EmbeddedSgroups.mrv", false, LoadAsMolOrRxn, 14, 14),
          MolTest("marvin03.mrv", false, LoadAsMolOrRxn, 31, 33),
          MolTest("MarvinBadMissingMolID.mrv", false, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinBadMissingAtomID.mrv", false, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinBadX2.mrv", false, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinBadY2.mrv", false, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinBadStereoGroupsAbs.mrv", false, LoadAsMolOrRxn, 8, 8),
          MolTest("MarvinBadElementType.mrv", false, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinBadMissingBondID.mrv", false, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinBadMissingBondAtomRefs", false, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinBadMissingBondOrder.mrv", false, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinBadMissingSruMolID.mrv", false, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinBadMissingSruID.mrv", false, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinBadMissingSruRole.mrv", false, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinBadMissingSruAtomRef.mrv", false, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinBadMissingSruTitle.mrv", false, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinBadSruAtomRef.mrv", false, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinBadSruID.mrv", false, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinBadSruRole.mrv", false, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinBadSruAtomRef.mrv", false, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinBadSruAtomRef.mrv", false, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinBadSruConnect.mrv", false, LoadAsMolOrRxn, 12, 11),
          MolTest("MarvinBadSupAttachAtom.mrv", false, LoadAsMolOrRxn, 9, 9),
          MolTest("MarvinBadSupAttachBond.mrv", false, LoadAsMolOrRxn, 9, 9),
          MolTest("MarvinBadSupAttachOrder.mrv", false, LoadAsMolOrRxn, 9, 9),
          MolTest("MarvinBadSupAttachAtom.mrv", false, LoadAsMolOrRxn, 9, 9),
          MolTest("MarvinBadSupAttachAtom.mrv", false, LoadAsMolOrRxn, 9, 9),
          MolTest("MarvinBadSupAttachAtom.mrv", false, LoadAsMolOrRxn, 9, 9),
          MolTest("MarvinBadSupMissingAttachBond.mrv", false, LoadAsMolOrRxn, 9, 9),
          MolTest("MarvinBadSupMissingAttachOrder.mrv", false, LoadAsMolOrRxn, 9,
                  9)};

  for (auto molFileTest : molFileTests) {
    BOOST_LOG(rdInfoLog) << "Test: " << molFileTest.fileName << std::endl;

    printf("Test\n\n %s\n\n", molFileTest.fileName.c_str());
    testMarvin(&molFileTest);
  }

  // now the reactions

  std::list<RxnTest> rxnFileTests{

          RxnTest("AlexRxn.mrv", true, LoadAsMolOrRxn, 1, 0, 1, 2, 0),
          RxnTest("BadReactionSign.mrv", true, LoadAsMolOrRxn, 2, 0, 1, 3, 0),
          RxnTest("bondArray_node.mrv", true, LoadAsMolOrRxn, 2, 4, 1, 3, 0),
          RxnTest("marvin03.mrv", true, LoadAsMolOrRxn, 1, 1, 1, 2, 0),
          RxnTest("marvin03.mrv", true, LoadAsRxn, 1, 1, 1, 2, 0),
          RxnTest("marvin03.mrv", false, LoadAsMol, 1, 1, 1, 2,
                  0),  // should fail
          RxnTest("marvin04.mrv", true, LoadAsMolOrRxn, 2, 1, 2, 4, 0),
          RxnTest("marvin08.mrv", true, LoadAsMolOrRxn, 2, 3, 2, 4, 0),
          RxnTest("marvin09.mrv", true, LoadAsMolOrRxn, 2, 3, 2, 4, 0),
          RxnTest("marvin11.mrv", true, LoadAsMolOrRxn, 2, 0, 1, 0, 0),
          RxnTest("marvin05.mrv", true, LoadAsMolOrRxn, 2, 1, 1, 3, 0),
          RxnTest("EmptyRxn.mrv", true, LoadAsMolOrRxn, 0, 0, 0, 0, 0),
          RxnTest("RxnNoCoords.mrv", true, LoadAsMolOrRxn, 2, 0, 1, 3, 0),
          RxnTest("mrvValenceZero.mrv", true, LoadAsMolOrRxn, 3, 0, 1, 4, 0),
          RxnTest("condition_coordinates_mpoint.mrv", true, LoadAsMolOrRxn, 1, 0, 1,
                  0, 0),
          RxnTest("marvin01.mrv", false, LoadAsMolOrRxn, 2, 1, 1, 3, 0),
          RxnTest("aspirineSynthesisWithAttributes.mrv", true, LoadAsMolOrRxn, 2, 0,
                  1, 3, 0)};

  for (auto rxnFileTest : rxnFileTests) {
    printf("Test\n\n %s\n\n", rxnFileTest.fileName.c_str());
    testMarvin(&rxnFileTest);
  }

  // now smiles tests

  std::list<SmilesTest> smiTests{

          SmilesTest("DoubleBondChain",
                     R"(CC1=C(\C=C\C(C)=C\C=C\C(C)=C/C(O)=O)C(C)(C)CCC1)", true, 22,
                     22),
          // SmilesTest(
          //     "Macrocycle2",
          //     R"(CC1OC(=O)CC(O)CC(O)CC(O)CCC(O)C(O)CC2(O)CC(O)C(C(CC(O[C@@H]3O[C@H](C)[C@@H](O)[C@H](N)[C@@H]3O)\C=C\C=C\C=C\C=C\CC\C=C\C=C\C(C)C(O)C1C)O2)C(O)=O
          //     |t:42,44,46,48,52,54|)",
          //     true, 65, 67),
          SmilesTest("Na_Mg_Al_OH",
                     "[OH-].[OH-].[OH-].[O--].[Na+].[Mg++].[Al+3].[Si].OC([O-])=O",
                     true, 12, 3),
          SmilesTest("Pb", "[Pb]", true, 1, 0),
          SmilesTest("O_Mg_Si", "[O].[Mg].[Si]", true, 3, 0),
          SmilesTest("SquiggleBond", "CN1N=C(SC1=NC(C)=O)S(N)(=O)=O |c:2|", true,
                     14, 14),
          SmilesTest(
                  "BigMacrocycle",
                  "C[C@@H]1CCCCCCCCC(=O)OCCN[C@H](C)CCCCCCCCC(=O)OCCN[C@H](C)CCCCCCCCC(=O)OCCN1",
                  true, 48, 48),
          SmilesTest("Smiles1", "N[C@@H]([O-])c1cc[13c]cc1", true, 9, 9)};

  for (auto smiTest : smiTests) {
    printf("Test\n\n %s\n\n", smiTest.name.c_str());
    // RDDepict::preferCoordGen = true;
    testSmilesToMarvin(&smiTest);
  }
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  RDLog::InitLogs();
  BOOST_LOG(rdInfoLog) << " ---- Running with POSIX locale ----- " << std::endl;

  RunTests();  // run with C locale

  return 0;
}
