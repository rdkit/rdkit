#include "RDGeneral/test.h"
#include "catch.hpp"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/DetermineBonds/DetermineBonds.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/AddHs.cpp>
#include <iostream>
#include <fstream>
#include <GraphMol/Resonance.h>

using namespace RDKit;

TEST_CASE("Determine Connectivity") {
  SECTION("Van der Waals") {
    unsigned int numTests = 39;
    for (unsigned int i = 0; i < numTests; i++) {
      std::string rdbase = getenv("RDBASE");
      std::string fName =
          rdbase + "/Code/GraphMol/DetermineBonds/test_data/connectivity/" +
          "test" + std::to_string(i) + ".xyz";
      std::unique_ptr<RWMol> mol(XYZFileToMol(fName));
      REQUIRE(mol);
      std::string smiles = mol->getProp<std::string>("_FileComments");
      std::unique_ptr<RWMol> orig(SmilesToMol(smiles));
      REQUIRE(orig);

      determineConnectivity(*mol, false);
      MolOps::removeAllHs(*mol, false);

      auto numAtoms = mol->getNumAtoms();

      REQUIRE(orig->getNumAtoms() == numAtoms);
      for (unsigned int i = 0; i < numAtoms; i++) {
        for (unsigned int j = i + 1; j < numAtoms; j++) {
          const auto origBond = orig->getBondBetweenAtoms(i, j);
          const auto molBond = mol->getBondBetweenAtoms(i, j);
          if (origBond) {
            CHECK(molBond);
          } else {
            CHECK(!molBond);
          }
        }
      }
    }
  }  // SECTION

  SECTION("Hueckel") {
    unsigned int numTests = 39;
    for (unsigned int i = 0; i < numTests; i++) {
      std::string rdbase = getenv("RDBASE");
      std::string fName =
          rdbase + "/Code/GraphMol/DetermineBonds/test_data/connectivity/" +
          "test" + std::to_string(i) + ".xyz";
      std::unique_ptr<RWMol> mol(XYZFileToMol(fName));
      REQUIRE(mol);
      std::string smiles = mol->getProp<std::string>("_FileComments");
      std::unique_ptr<RWMol> orig(SmilesToMol(smiles));
      REQUIRE(orig);
      int charge = MolOps::getFormalCharge(*orig);

      determineConnectivity(*mol, true, charge);
      MolOps::removeAllHs(*mol, false);

      auto numAtoms = mol->getNumAtoms();

      REQUIRE(orig->getNumAtoms() == numAtoms);
      for (unsigned int i = 0; i < numAtoms; i++) {
        for (unsigned int j = i + 1; j < numAtoms; j++) {
          const auto origBond = orig->getBondBetweenAtoms(i, j);
          const auto molBond = mol->getBondBetweenAtoms(i, j);
          if (origBond) {
            CHECK(molBond);
          } else {
            CHECK(!molBond);
          }
        }
      }
    }
  }  // SECTION

  SECTION("DetermineBondOrdering using charged fragments") {
    unsigned int numTests = 38;
    for (unsigned int i = 0; i < numTests; i++) {
      std::string rdbase = getenv("RDBASE");
      std::string fName =
          rdbase +
          "/Code/GraphMol/DetermineBonds/test_data/charged_fragments/" +
          "test" + std::to_string(i) + ".xyz";
      std::unique_ptr<RWMol> mol(XYZFileToMol(fName));
      REQUIRE(mol);
      std::string smiles = mol->getProp<std::string>("_FileComments");
      std::unique_ptr<RWMol> orig(SmilesToMol(smiles));
      REQUIRE(orig);
      SmilesWriteParams params = {false, false, true, false, false, false, -1};
      std::string canonSmiles = MolToSmiles(*orig, params);
      int charge = MolOps::getFormalCharge(*orig);

      determineBonds(*mol, false, charge);
      MolOps::removeAllHs(*mol, false);

      auto numAtoms = mol->getNumAtoms();
      REQUIRE(orig->getNumAtoms() == numAtoms);

      ResonanceMolSupplier resMolSuppl(
          *mol, ResonanceMolSupplier::UNCONSTRAINED_CATIONS |
                    ResonanceMolSupplier::UNCONSTRAINED_ANIONS);
      bool valid = false;
      for (unsigned int i = 0; i < resMolSuppl.length(); i++) {
        std::unique_ptr<ROMol> firstResMol(resMolSuppl[i]);
        std::unique_ptr<RWMol> resMol(new RWMol(*firstResMol));
        MolOps::setAromaticity(*resMol);

        std::string molSmiles = MolToSmiles(*resMol, params);
        if (molSmiles == canonSmiles) {
          CHECK(true);
          valid = true;
          break;
        }
      }
      if (!valid) {
        CHECK(false);
      }
    }
  }  // SECTION

  SECTION("DetermineBondOrdering using radicals") {
    unsigned int numTests = 10;
    for (unsigned int i = 0; i < numTests; i++) {
      std::string rdbase = getenv("RDBASE");
      std::string fName = rdbase +
                          "/Code/GraphMol/DetermineBonds/test_data/radicals/" +
                          "test" + std::to_string(i) + ".xyz";
      std::unique_ptr<RWMol> mol(XYZFileToMol(fName));
      REQUIRE(mol);
      std::string smiles = mol->getProp<std::string>("_FileComments");
      std::unique_ptr<RWMol> orig(SmilesToMol(smiles));
      REQUIRE(orig);
      SmilesWriteParams params = {false, false, true, false, false, false, -1};
      std::string canonSmiles = MolToSmiles(*orig, params);
      int charge = MolOps::getFormalCharge(*orig);

      determineBonds(*mol, false, charge, 1.3, false);
      MolOps::removeAllHs(*mol, false);

      auto numAtoms = mol->getNumAtoms();
      REQUIRE(orig->getNumAtoms() == numAtoms);

      ResonanceMolSupplier resMolSuppl(
          *mol, ResonanceMolSupplier::UNCONSTRAINED_CATIONS |
                    ResonanceMolSupplier::UNCONSTRAINED_ANIONS);
      bool valid = false;
      for (unsigned int i = 0; i < resMolSuppl.length(); i++) {
        std::unique_ptr<ROMol> firstResMol(resMolSuppl[i]);
        std::unique_ptr<RWMol> resMol(new RWMol(*firstResMol));
        MolOps::setAromaticity(*resMol);

        std::string molSmiles = MolToSmiles(*resMol, params);
        if (molSmiles == canonSmiles) {
          CHECK(true);
          valid = true;
          break;
        }
      }
      if (!valid) {
        CHECK(false);
      }
    }
  }  // SECTION
}

//        std::string smiles[] = {
//            "C[C-](c1ccccc1)C",
//            "C[C-](C)c1ccccc1",
//            "C=C([O-])CC",
//            "C=C([NH3+])CC",
//            "CC(=O)[O-]",
//            "C[N+](=O)[O-]",
//            "CS(CC)(=O)=O",
//            "CS([O-])(=O)=O",
//            "C=C(C)CC",
//            "CC(C)CC",
//            "C=C(N)CC",
//            "C=C(C)C=C",
//            "C#CC=C",
//            "c1ccccc1",
//            "c1ccccc1c1ccccc1",
//            "[NH3+]CS([O-])(=O)=O",
//            "CC(NC)=O",
//            "[O-]c1ccccc1",
//            "O=C(C=C1)C=CC1=CCC([O-])=O",
//            "C#CC#C",
//            "Cc1ccc(cc1)C1C=CC2C(C=CC2(C#N)C#N)=CC=1",
//            "C[NH+]=C([O-])CC[NH+]=C([O-])C",
//            "C[NH+]=CC=C([O-])C",
//            "[C+](C)(C)CC[C-](C)(C)",
//            "O=C(C=C1)C=CC1=CCC([O-])=O",
//            "O=C([CH-]C=CC(C([O-])=O)=O)[O-]",
//            "[O-]c1ccccc1",
//            "CNC(C(C)=[NH+][CH-]CC(O)=O)=O",
//            "[CH2][CH2][CH]=[CH][CH2]",
//            "Cc1ccc(cc1)C1C=CC2C(C=CC2(C#N)C#N)=CC=1",
//            "CC1C=CC2C(C=CC2(C)C)=CC=1",
//            "CC1=CC=C(C=CC2)C2C=C1",
//            "CC1=CC=C(C2=CC=CC=C2)C=C1",
//            "C1(CC2=CC=CC=C2)=CC=CC=C1",
//            "[O-]c1ccccc1[O-]",
//            "C[N+](=O)[O-]",
//            "N#CC(C#N)=CC=C1C=CC=CC(=C1)c1ccc(cc1)[N+](=O)[O-]",
//            "CNC([O-])=C([NH+]=C/CC(O)=O)C",
//            "Cc1cn(C2CC(O)C(COP(=O)([O-])OP(=O)([O-])OC3OC(C)C([NH3+])C(O)C3O)O2)c(=O)[nH]c1=O"
//        };
